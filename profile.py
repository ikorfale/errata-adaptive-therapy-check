"""Step 4b: identifiability by profile over the resistant start R0.
For each patient, training data = points before the cycle-2 stop (as in step 3). R0 is fixed on a grid
(fraction f of the first PSA: 1e-4 .. 0.3), the other six parameters are refit (4 random starts + a warm start from the neighbour, swept up and down).
A grid point is 'equally good' if its training cost is within the 95% profile-likelihood band:
  (cost_f - cost_best) / sigma^2 <= 1.92, sigma^2 = 2*cost_best/n_train  (cost = 0.5*sum r^2).
For every grid point record: training cost, forecast RMSE on the held-out cycle-2 regrowth, and two
futures from the cycle-2 stop: PSA after 2 years of CONTINUOUS suppression, and days until PSA under
continuous suppression first rises 25% and 2 ng/ml above its running nadir (PCWG-like; cap 1500). Spread across equally good
points = how much the future is not pinned down by the data.
Usage: profile.py [ids...]; writes profile_<tag>.json, prints one line per patient."""
import sys, json, time
import numpy as np
from scipy.optimize import least_squares
from load import patients
from model import series, simulate, resid, LO, HI, EPS, rmse
from forecast import phase_start
FR = np.logspace(-4, np.log10(0.3), 9)

def fit_fixed(t, y, D, lR0, starts, rng, warm=None):
    lo, hi = LO[:6], HI[:6]
    f = lambda q: resid(np.r_[q, lR0], t, y, D)
    best = None
    for k in range(starts + (warm is not None)):
        if warm is not None and k == 0: q0 = warm.copy()
        else: q0 = lo + rng.random(6) * (hi - lo); q0[5] = np.log(max(y[0], 0.2))
        q0 = np.clip(q0, lo + 1e-6, hi - 1e-6)
        try: s = least_squares(f, q0, bounds=(lo, hi), max_nfev=300)
        except Exception: continue
        if best is None or s.cost < best.cost: best = s
    return np.r_[best.x, lR0], best.cost

def one(pid, rows, starts=4):
    rng = np.random.default_rng(1)
    t, y, D = series(rows)
    s2 = int(phase_start(rows, 2, 0)); e2 = phase_start(rows, 3, 1); e2 = e2 if e2 is not None else t.max() + 1
    tr = t < s2; te = (t >= s2) & (t < e2)
    H = 1500; Dc = np.r_[D[:s2], np.ones(H)]; Df = np.r_[D, np.zeros(max(0, s2 + H - len(D)))]
    nadir = y[tr].min()
    fits = {}
    for order in (range(len(FR)), reversed(range(len(FR)))):  # two warm-started sweeps, keep the better fit
        warm = None
        for i in order:
            lR0 = np.clip(np.log(FR[i] * max(y[0], 0.2)), LO[6], HI[6])
            p, c = fit_fixed(t[tr], y[tr], D, lR0, starts, rng, warm)
            if i not in fits or c < fits[i][1]: fits[i] = (p, c)
            warm = fits[i][0][:6]
    pts = []
    for i, f in enumerate(FR):
        p, c = fits[i]
        m = simulate(p, Df, s2 + H)
        fc = rmse(np.log(m[t[te]] + EPS) - np.log(y[te] + EPS))
        mc = simulate(p, Dc, s2 + H)[s2:]
        run = np.minimum.accumulate(mc)  # PCWG-like progression: 25% and >=2 ng/ml above the running nadir
        over = np.nonzero((mc >= 1.25 * run) & (mc - run >= 2))[0]
        pts.append(dict(f=float(f), cost=float(c), fc=fc, psa2y=float(mc[730]), ttp=int(over[0]) if len(over) else H))
    cb = min(q['cost'] for q in pts); s2v = 2 * cb / tr.sum()
    for q in pts: q['ok'] = bool((q['cost'] - cb) / max(s2v, 1e-9) <= 1.92); q['ok20'] = bool(q['cost'] <= 1.2 * cb)
    ok = [q for q in pts if q['ok']]
    return dict(pid=pid, n_train=int(tr.sum()), pts=pts, n_ok=len(ok),
                ttp_lo=min(q['ttp'] for q in ok), ttp_hi=max(q['ttp'] for q in ok),
                fc_lo=min(q['fc'] for q in ok), fc_hi=max(q['fc'] for q in ok),
                ttp20=[q['ttp'] for q in pts if q['ok20']])

if __name__ == '__main__':
    P = patients(); ids = sys.argv[1:] or sorted(P); out = []
    for pid in ids:
        t0 = time.time()
        try: r = one(pid, P[pid])
        except Exception as e: print(pid, 'SKIP', type(e).__name__, e, flush=True); continue
        out.append(r)
        print('%s ok %d/9  ttp %d..%d days  forecast rmse %.2f..%.2f  %.0fs' % (pid, r['n_ok'], r['ttp_lo'], r['ttp_hi'], r['fc_lo'], r['fc_hi'], time.time() - t0), flush=True)
        json.dump(out, open('profile_%s.json' % ('all' if not sys.argv[1:] else '_'.join(ids)), 'w'), indent=1)
