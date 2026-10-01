"""Two-population Lotka-Volterra model of adaptive therapy, in the form used by Zhang et al. 2017
and Strobl et al. 2021 (Cancer Research) for this same trial:
  dS/dt = rS*S*(1-(S+R)/K)*(1-dD*D(t)) - dT*S     (drug-sensitive)
  dR/dt = rR*R*(1-(S+R)/K)            - dT*R     (resistant)
  PSA   = S + R                                   (PSA units)
D(t) = 1 while the patient is on androgen suppression, 0 off. Integrated in log space, 1-day steps.
Residuals are in log(PSA + 0.1), so a 0.3 residual is roughly a 35% miss."""
import numpy as np
from scipy.optimize import least_squares

EPS = 0.1
NAMES = ['rS', 'rR', 'dT', 'dD', 'K', 'S0', 'R0']
LO = np.log([1e-4, 1e-5, 1e-5, 1.0, 1.0, 1e-2, 1e-4])
HI = np.log([0.5, 0.5, 0.2, 20., 2000., 2000., 500.])

def series(rows):
    """Observed points and a daily on/off schedule."""
    T = int(rows[-1]['t']) + 1
    D = np.zeros(T)
    for a, b in zip(rows, rows[1:] + [None]):
        end = int(b['t']) if b else T
        D[int(a['t']):end] = a['on'] or 0
    obs = [(int(r['t']), r['psa']) for r in rows if r['psa'] is not None]
    t = np.array([o[0] for o in obs]); y = np.array([o[1] for o in obs])
    return t, y, D

def simulate(p, D, T):
    rS, rR, dT, dD, K, S0, R0 = np.exp(p)
    lS, lR = np.log(S0), np.log(R0)
    out = np.empty(T)
    for i in range(T):
        S, R = np.exp(lS), np.exp(lR)
        out[i] = S + R
        g = 1 - (S + R) / K
        lS += rS * g * (1 - dD * D[i]) - dT
        lR += rR * g - dT
        lS = max(lS, -40.); lR = max(lR, -40.)
    return out

def resid(p, t, y, D):
    T = int(t.max()) + 1
    m = simulate(p, D, T)
    return np.log(m[t] + EPS) - np.log(y + EPS)

def fit(t, y, D, starts=12, seed=0):
    rng = np.random.default_rng(seed)
    best = None; sols = []
    for k in range(starts):
        p0 = LO + rng.random(len(LO)) * (HI - LO)
        p0[5] = np.log(max(y[0], 0.2))  # start near the first PSA
        p0 = np.clip(p0, LO + 1e-6, HI - 1e-6)
        try:
            s = least_squares(resid, p0, bounds=(LO, HI), args=(t, y, D), max_nfev=400)
        except Exception:
            continue
        sols.append(s)
        if best is None or s.cost < best.cost: best = s
    return best, sols

def rmse(r): return float(np.sqrt(np.mean(np.square(r))))
