# GTA VI Twitch CCU Forecast 01 — entry by @ris

Forecasts how many people watch Grand Theft Auto VI live on Twitch (category
1867279146) in each of hours 2–24 of its release day, 2026-11-19 UTC, as five
quantiles (5%, 25%, 50%, 75%, 95%) an hour, for the make.fun bounty
[gta6-twitch-ccu-forecast-01](https://make.fun/bounty/gta6-twitch-ccu-forecast-01).

[![make.fun](https://make.fun/badge/gta6-twitch-ccu-forecast-01.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01)
[![rank](https://make.fun/badge/gta6-twitch-ccu-forecast-01/Sipalingtesnet/gta6-twitch-ccu-forecast-01-entry-ris/rank.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01)
[![reward](https://make.fun/badge/gta6-twitch-ccu-forecast-01/Sipalingtesnet/gta6-twitch-ccu-forecast-01-entry-ris/reward.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01)

## Entry

`main.py`, Python 3.13, standard library only, no network, deterministic.  It
reads every case in `$MAKEFUN_INPUT_DIR` (else `/input`) and writes one forecast
per case to `$MAKEFUN_OUTPUT_DIR` (else `/output`).

## Model

The forecast scales the case's hour-1 print `p1` by a ratio curve
`r_h = viewers(h) / viewers(1)` per UTC hour, estimated from the past big
launches in the case's `comparables`:

1. **Ratio series.**  Each comparable launch contributes `ln(viewers(h)/p1)`
   for hours 2–24 from its own release day.  A comparable with no hour-1 print
   (or a print below 1) is skipped.  GTA VI's own category history (the
   2026-08-27 extended-look premiere onward) is **excluded from the ratio
   pool**: it is pre-launch hype data, not a launch-day shape — its category
   held a few hundred viewers an hour, against a launch-day print in the
   hundreds of thousands, so its ratios would be meaningless noise.

2. **p1 conditioning.**  Launches are pooled by `log(p1)` quartile.  This cleanly
   separates the day shapes that a mega launch and a tiny launch really have: a
   launch whose hour-1 print is small (pre-unlock noise) can grow a hundredfold
   over the day, while a launch whose hour-1 print is already near its peak
   stays near `r_h = 1` all day.  For the case's own quartile the per-hour
   **median** of the pooled log ratios is the forecast median ratio, clamped to
   `[0.05, 60]` so one weird launch cannot blow the forecast up.

3. **Smoothing.**  The per-hour medians are smoothed across hours with a
   Gaussian kernel (bandwidth 2 h) so the forecast keeps a smooth diurnal curve
   (overnight dip, US-unlock rise, late-UTC-evening peak, midnight fall).

4. **Spread.**  The log-spread `s_h` is the pooled interquartile range of the
   log ratios at that hour over 1.349 (the normal-equivalent width), floored so
   the band is never absurdly narrow, times a calibration constant 1.4.  The
   constant was chosen so the 5%–95% band covers about 90% of scored hours.

5. **Quantiles.**  `q_tau = p1 * r_h * exp(z_tau * s_h)` with
   `z = [-1.645, -0.674, 0, 0.674, 1.645]`, rounded half-up to 2 decimals,
   clamped to `(0, 1e9]` and kept non-decreasing.

Everything above was calibrated on leave-one-launch-out pseudo-cases built from
the 38 comparables (each comparable scored once as a case with itself excluded
from calibration): the entry's mean pinball loss is about **0.45** against the
house persistence baseline's **0.93** on the same pseudo-cases, with ~88% 90%
coverage.  The 4 practice cases p00–p03 are used only for format, runtime,
determinism and sanity checks.

## Run

```bash
mkdir -p out
MAKEFUN_INPUT_DIR=input MAKEFUN_OUTPUT_DIR=out python3 main.py
```