"""Task #78. Zhang adaptive abiraterone arm: does a fitted two-population model predict the next off-period better
than 'same as last' (54.0 d MAE) and 'last x 0.72' (39.2 d, LOPO ratio)? Same 36 next cycles as cycles.py.
Causal: at the stop that opens off-span j, fit on every visit up to and including that stop day; predict the day the
model PSA first reaches the patient's restart level (mean PSA at their earlier restart visits); off_hat = that - stop.
Model (daily Euler, log space): g = 1 - (S+R)/K; lS += rS g - dD D; lR += rR g; PSA = S + R. D = abi of last visit."""
import csv, sys, numpy as np
from collections import defaultdict
from scipy.optimize import least_squares
FLOOR, MAXOFF = 0.05, 1500
import os; THR = os.environ.get('THR', 'mean')   # restart level: mean / min of earlier restart PSAs, or day-0 PSA
P = defaultdict(list)
for r in csv.DictReader(open('../data/zhang/zhang_long.csv')):
    if r['arm'] == 'Adaptive': P[r['pid']].append((float(r['day']), float(r['psa']), int(float(r['abi']))))
def sim(th, dose, n):
    rS, rR, dD, K, S0, R0 = np.exp(th); lS, lR = np.log(S0), np.log(R0); out = np.empty(n)
    for t in range(n):
        S, R = np.exp(lS), np.exp(lR); out[t] = S + R; g = 1 - (S + R) / K
        lS = max(lS + rS * g - dD * dose[t], -30); lR = max(lR + rR * g, -30)
    return out
def fit(days, psa, dose):
    n = int(days[-1]) + 1; y = np.log(np.maximum(psa, FLOOR)); di = days.astype(int)
    res = lambda th: np.log(np.maximum(sim(th, dose, n)[di], FLOOR)) - y
    lo = np.log([1e-3, 1e-4, 1e-3, psa.max(), 1e-4, 1e-5]); hi = np.log([0.3, 0.1, 1.0, 100 * psa.max() + 1, 10 * psa[0] + 1, psa[0] + 1])
    best = None
    for rS, rR, dD, Kx, fr in ((0.03, 0.005, 0.1, 3, 0.99), (0.08, 0.01, 0.3, 10, 0.9), (0.01, 0.002, 0.05, 2, 0.999),
                               (0.05, 0.02, 0.2, 30, 0.5)):
        x0 = np.clip(np.log([rS, rR, dD, Kx * psa.max(), fr * psa[0] + 1e-4, (1 - fr) * psa[0] + 1e-5]), lo + 1e-6, hi - 1e-6)
        r = least_squares(res, x0, bounds=(lo, hi), max_nfev=300)
        if best is None or r.cost < best.cost: best = r
    return best
rows = []
for pid, v in sorted(P.items()):
    v.sort(); d = np.array([x[0] for x in v]); p = np.array([x[1] for x in v]); a = np.array([x[2] for x in v])
    ch = [i for i in range(1, len(a)) if a[i] != a[i - 1]]             # visit index where abi changes
    stops = [i for i in ch if a[i] == 0]; starts = [i for i in ch if a[i] == 1]
    # complete cycles as in cycles.py: on-span closed by a stop, off-span closed by a restart
    offs = [(s, next(t for t in starts if t > s)) for s in stops if any(t > s for t in starts)]
    for j in range(1, len(offs)):
        s, t = offs[j]; prev_restart = [p[t2] for (_, t2) in offs[:j]]; thr = {'mean': float(np.mean(prev_restart)), 'min': float(np.min(prev_restart)), 'base': float(p[0])}[THR]
        dose = np.zeros(int(d[s]) + MAXOFF + 1)
        for i in range(len(d)):
            if d[i] <= d[s]: dose[int(d[i]):] = a[i]
        dose[int(d[s]):] = 0.0
        r = fit(d[:s + 1], p[:s + 1], dose)
        traj = sim(r.x, dose, len(dose))[int(d[s]):]
        hit = np.nonzero(traj >= thr)[0]; off_hat = int(hit[0]) if len(hit) else MAXOFF
        off = d[t] - d[s]; last = d[offs[j - 1][1]] - d[offs[j - 1][0]]
        rows.append((pid, j, off, off_hat, last, thr, float(np.sqrt(2 * r.cost / (s + 1)))))
        print(f"{pid} cycle {j} off {off:.0f} model {off_hat} last {last:.0f} thr {thr:.2f} fit_rmse_log {rows[-1][-1]:.3f}", flush=True)
R = np.array([(r[2], r[3], r[4]) for r in rows])
print(f"n {len(R)} | MAE model {np.mean(abs(R[:,1]-R[:,0])):.1f} d | same-as-last {np.mean(abs(R[:,2]-R[:,0])):.1f} d | "
      f"median |err| model {np.median(abs(R[:,1]-R[:,0])):.0f} last {np.median(abs(R[:,2]-R[:,0])):.0f} | "
      f"model closer in {int(sum(abs(R[:,1]-R[:,0]) < abs(R[:,2]-R[:,0])))}/{len(R)} | capped at {MAXOFF}: {int(sum(R[:,1]>=MAXOFF))}")
