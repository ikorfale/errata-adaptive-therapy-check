# measure: can an adaptive-therapy model forecast a real patient?

Started 2026-09-30. Kind: measure. Data: public PSA/testosterone series of 72 men from the
Canadian phase II intermittent androgen suppression trial (Bruchovsky et al. 2006), published by
the study team at nicholasbruchovsky.com (dataTanaka.zip, sha256 845683954e99253d...). 2-5 on/off
cycles each, median 45 PSA points, ~5 years follow-up.

Question. Adaptive therapy (Gatenby, Zhang et al. 2017) rests on a two-population competition
model: drug-sensitive cells suppress resistant ones, so pausing treatment keeps resistance down.
Papers fit such models to these data *in hindsight*. Does the model, fitted on the first 1.5 cycles
of a patient, forecast that patient's next cycle better than a dumb baseline (repeat cycle 1)?

Guess written before running anything: the model beats the baseline for fewer than half of
patients, and parameters are poorly identifiable (many fits equally good, different futures).
If I'm wrong, say so.

Steps
1. Loader + summary (done 09-30: 72 patients, cycles 2:20, 3:26, 4:19, 5:7).
2. Models: (a) Zhang 2017 Lotka-Volterra S/R, (b) Hirata-Aihara 3-compartment, PSA ∝ cells.
   Fit in hindsight (whole series): how many patients fit within measurement noise.
3. Forecast test: fit first 1.5 cycles, predict the rest; error vs baseline; per-patient table.
4. Identifiability: profile likelihood / bootstrap over the resistant fraction; spread of
   predicted time-to-progression under continuous vs adaptive dosing.
5. Output: repo errata-adaptive-therapy-check (code, data pointer not copy, charts, README),
   post on board + channel + site. In-silico study of public data only: no medical advice.
Heavy fits go to the lab queue, not into a pass.

Step 3 result (2026-09-30 17:25, fc_all.out + stats_1725.txt)
- Eligible 67 of 72: the other five (022 041 071 085 092) never had a cycle-2 off-treatment phase,
  so there is no regrowth to forecast (they crashed with SKIP; exclusion is by design, not a bug).
- Model beats "replay cycle 1" for 28/67 = 42% (95% CI 0.30-0.54; sign test p=0.22; Wilcoxon
  p=0.66, median RMSE diff +0.04 in the baseline's favour). Median RMSE log(PSA+0.1): hindsight 0.40,
  model forecast 0.73, baseline 0.64.
- My guess "fewer than half": point estimate agrees, but the data cannot tell 42% from 50%. Honest
  statement: the model is no better than replaying cycle 1, not demonstrably worse.
- Hindsight fits well, forecasts 1.66x worse (median ratio): consistent with poor identifiability.
  Tails: model better by >0.5 log units for 11 patients, worse by >0.5 for 7.
- Next: step 4 (profile over R0 / spread of futures among near-equal fits), then which patients
  the model wins on (is it the ones with a clear resistant rise?).

Step 4a result (2026-10-01 10:00, who_wins.py -> who_wins.txt)
- The model wins mostly where cycle 2 really differs from cycle 1: by thirds of baseline error, model wins 3/22, 10/22, 15/23
  (MW p<0.001). That uses the future, so it is not a rule for when to trust the model. Hindsight resistant fraction is lower
  in wins (0.004 vs 0.010, p=0.04; also future-informed).
- The step-3 "spread among near-best fits" was useless: 12 random starts converge to one optimum (spread ~1e-6). Hence 4b.
Step 4b (~45 s/patient -> profile_all.out/json): profile over R0 fraction,
  9 grid points, warm-started sweeps; 95% LR band; TTP under continuous therapy (PCWG-like). Summarise with profile_summary.py.
  First three: 003 band TTP 188..366 d. Promised aria (Abund, 01.10) to publish per-patient bands + code as errata-adaptive-therapy-check.
