"""Forecast GTA VI's launch-day Twitch viewers hour by hour (make.fun bounty).

For every case this program is given (/input/<id>.in, or $MAKEFUN_INPUT_DIR) it
writes /output/<id>.out: for each of hours 2..24 of the game's release day, five
quantiles (0.05, 0.25, 0.5, 0.75, 0.95) of that hour's Twitch viewers.

Model (written from scratch; Python 3.13 standard library only, no network,
deterministic):

  1. Every usable comparable launch contributes a log-ratio series
     lr_h = ln(viewers(h) / viewers(1)) for h = 2..24 from its hourly data.
     A comparable with no hour-1 print (or a print below 1) is skipped, and so
     is GTA VI's own pre-release category history: that is hype data, not a
     launch-day shape (the game's category had a few hundred viewers an hour
     until the premiere weeks, against a launch-day print in the hundreds of
     thousands).

  2. Launches are pooled by log(p1) quartile, which cleanly separates the
     wildly different day shapes of tiny-p1 launches (pre-unlock noise, ratios
     that explode a hundredfold) from big-p1 launches (peak near hour 1, ratio
     near 1 all day).  For the case's own quartile the per-hour MEDIAN of the
     pooled log ratios is the forecast median ratio r_h (clamped to
     [0.05, 60] so one weird launch cannot blow the forecast up).

  3. The per-hour medians are smoothed across hours with a Gaussian kernel
     (bandwidth 2 h) so the forecast keeps a smooth diurnal curve and the
     per-hour estimates borrow strength from their neighbours.

  4. The log-spread s_h is the pooled interquartile range of the log ratios at
     that hour over 1.349 (the normal-equivalent width), floored so the band is
     never absurdly narrow, times a calibration constant (1.4).  The constant
     was chosen so the 5%-95% band covers about 90% of scored hours on
     leave-one-launch-out pseudo-cases built from the comparables.

  5. q_tau = p1 * r_h * exp(z_tau * s_h) with z = [-1.645, -0.674, 0, 0.674,
     1.645], rounded half-up to 2 decimals.  Quantiles are clamped to
     (0, 1e9] and kept non-decreasing.
"""

import json
import math
import os
import statistics as st
from datetime import datetime
from pathlib import Path

LEVELS = [0.05, 0.25, 0.5, 0.75, 0.95]
Z = [-1.645, -0.674, 0.0, 0.674, 1.645]
HOURS = list(range(2, 25))

# Model hyper-parameters (calibrated on the comparables, leave-one-out).
N_BUCKETS = 4          # log-p1 quartile buckets
SMOOTH_BW = 2.0        # Gaussian smoothing bandwidth over hours (h)
SPREAD_MULT = 1.4      # spread calibration constant
R_MIN = 0.05           # ratio clamp floor
R_MAX = 60.0           # ratio clamp ceiling
MIN_SCALE = 0.15       # log-spread floor per hour
FALLBACK_SCALE = 0.5   # log-spread when a bucket has too few launches

# GTA VI's own category history is pre-launch hype, not a launch-day shape.
SKIP_NAME_TAGS = ("grand theft auto vi category",)

# Fallback diurnal prior (median ratio to hour 1): dips overnight, rises as the
# US wakes, peaks in the late UTC evening, falls towards midnight.  Used only
# when a case carries no usable comparables at all.
DIURNAL_PRIOR = {
    2: 0.87, 3: 0.75, 4: 0.73, 5: 0.75, 6: 0.79, 7: 0.80, 8: 0.90,
    9: 0.88, 10: 0.91, 11: 0.93, 12: 0.98, 13: 1.06, 14: 1.18, 15: 1.52,
    16: 1.65, 17: 1.77, 18: 1.95, 19: 2.40, 20: 2.37, 21: 2.33, 22: 2.30,
    23: 2.18, 24: 2.04,
}
PRIOR_SCALE = 0.8      # log-spread of the fallback band


def parse_iso(s: str) -> float:
    """ISO 8601 UTC instant to a Unix timestamp."""
    return datetime.fromisoformat(s.replace("Z", "+00:00")).timestamp()


def round2(x: float):
    """Round half-up to 2 decimals; a whole number comes back as an int."""
    scaled = x * 100.0
    whole = math.floor(scaled)
    if scaled - whole >= 0.5:
        whole += 1
    value = whole / 100.0
    return int(value) if abs(value - round(value)) < 1e-9 else value


def comparable_series(comp):
    """(p1, {h: ln(v_h / p1)}) for one comparable, or (None, None) if unusable."""
    name = (comp.get("name") or "").lower()
    if any(tag in name for tag in SKIP_NAME_TAGS):
        return None, None
    day_start = parse_iso(comp["day_start"])
    by_hour = {}
    p1 = None
    for pt in comp.get("hourly") or []:
        v = pt.get("viewers")
        if v is None:
            continue
        k = (parse_iso(pt["from"]) - day_start) / 3600.0
        h = int(math.floor(k)) + 1
        if 1 <= h <= 24:
            if h not in by_hour:
                by_hour[h] = v
        if abs(k) < 0.5:
            p1 = v
    if p1 is None:
        p1 = by_hour.get(1)
    if p1 is None or p1 < 1:
        return None, None
    p1 = max(1.0, float(p1))
    lr = {}
    for h in HOURS:
        v = by_hour.get(h)
        if v is not None and v > 0:
            lr[h] = math.log(max(1.0, float(v)) / p1)
    if not lr:
        return None, None
    return p1, lr


def bucket_edges(log_p1s):
    """Quantile edges over the pool's log p1 values."""
    s = sorted(log_p1s)
    edges = []
    for i in range(1, N_BUCKETS):
        edges.append(s[int(round(i * (len(s) - 1) / N_BUCKETS))])
    return edges


def bucket_of(lp, edges):
    for i, e in enumerate(edges):
        if lp <= e:
            return i
    return len(edges)


def gauss_smooth(vals, bw):
    """Gaussian-smooth a {hour: value} map over hours 2..24."""
    if bw <= 0:
        return dict(vals)
    out = {}
    for h in HOURS:
        num, den = 0.0, 0.0
        for h2 in HOURS:
            if h2 in vals:
                w = math.exp(-0.5 * ((h - h2) / bw) ** 2)
                num += w * vals[h2]
                den += w
        if den > 0:
            out[h] = num / den
    return out


def spread_scale(vals):
    """Robust log-spread: normal-equivalent IQR, floored."""
    if len(vals) < 3:
        return FALLBACK_SCALE
    s = sorted(vals)
    lo = s[max(0, int(0.25 * (len(s) - 1)))]
    hi = s[min(len(s) - 1, int(0.75 * (len(s) - 1)))]
    return max(MIN_SCALE, (hi - lo) / 1.349)


def fit_curve(pool, case_p1):
    """Per-hour median log-ratio and log-spread for the case's p1 quartile.

    pool: list of (p1, lr) comparable series.  Falls back to the whole-pool
    curve when the case's own quartile is empty.
    """
    lp = math.log(max(1.0, case_p1))
    edges = bucket_edges([math.log(p) for p, _ in pool])
    b = bucket_of(lp, edges)
    meds, scales = {}, {}
    for h in HOURS:
        vals = [lr[h] for p, lr in pool
                if bucket_of(math.log(max(1.0, p)), edges) == b and h in lr]
        if not vals:
            vals = [lr[h] for _, lr in pool if h in lr]
        if not vals:
            continue
        meds[h] = st.median(vals)
        scales[h] = spread_scale(vals)
    meds = gauss_smooth(meds, SMOOTH_BW)
    return meds, scales


def forecast_case(case) -> dict:
    first = next((k for k in case["known"] if k["hour"] == 1), None)
    if first is None:
        raise ValueError("no hour-1 print")
    p1 = max(1.0, float(first["viewers"]))

    pool = []
    for comp in case.get("comparables") or []:
        sp1, lr = comparable_series(comp)
        if sp1 is not None:
            pool.append((sp1, lr))

    if pool:
        meds, scales = fit_curve(pool, p1)
    else:
        meds, scales = {}, {}

    rows = []
    for h in HOURS:
        if h in meds:
            r = min(R_MAX, max(R_MIN, math.exp(meds[h])))
            s = scales.get(h, FALLBACK_SCALE) * SPREAD_MULT
        else:
            r = DIURNAL_PRIOR.get(h, 1.0)
            s = PRIOR_SCALE
        q = [round2(min(1e9, p1 * r * math.exp(z * s))) for z in Z]
        for i in range(1, 5):  # keep the five quantiles non-decreasing after rounding
            if q[i] < q[i - 1]:
                q[i] = q[i - 1]
        rows.append({"hour": h, "q": q})
    return {"hours": rows}


def main() -> None:
    in_dir = Path(os.environ.get("MAKEFUN_INPUT_DIR", "/input"))
    out_dir = Path(os.environ.get("MAKEFUN_OUTPUT_DIR", "/output"))
    out_dir.mkdir(parents=True, exist_ok=True)
    for path in sorted(in_dir.glob("*.in")):
        case = json.loads(path.read_text(encoding="utf-8"))
        forecast = forecast_case(case)
        (out_dir / f"{path.stem}.out").write_text(
            json.dumps(forecast, separators=(",", ":")) + "\n", encoding="utf-8"
        )


if __name__ == "__main__":
    main()