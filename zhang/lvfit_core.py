"""Task #78. Zhang adaptive abiraterone arm: does a fitted two-population model predict the next off-period better
than 'same as last' (54.0 d MAE) and 'last x 0.72' (39.2 d, LOPO ratio)? Same 36 next cycles as cycles.py.
Causal: at the stop that opens off-span j, fit on every visit up to and including that stop day; predict the day the
model PSA first reaches the patient's restart level (mean PSA at their earlier restart visits); off_hat = that - stop.
Model (daily Euler, log space): g = 1 - (S+R)/K; lS += rS g - dD D; lR += rR g; PSA = S + R. D = abi of last visit."""
import csv, sys, numpy as np
from collections import defaultdict
from scipy.optimize import least_squares
FLOOR, MAXOFF = 0.05, 1500
   # restart level: mean / min of earlier restart PSAs, or day-0 PSA
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
