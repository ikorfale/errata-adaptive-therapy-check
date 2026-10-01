"""Step 2+3 for one or more patients. For each patient:
  hindsight: fit the model to the whole series, RMSE in log(PSA+0.1);
  forecast:  fit only on points before the cycle-2 treatment stop ("1.5 cycles"), then predict the
             cycle-2 off-treatment regrowth. Baseline: the cycle-1 off-treatment curve, replayed from
             the cycle-2 stop (same days since stop). Both scored on the same held-out points.
Usage: forecast.py [patient ids...] (default: all). Prints one line per patient, then a summary."""
import sys, time, json
import numpy as np
from load import patients
from model import series, fit, resid, simulate, rmse, EPS

def phase_start(rows, cycle, on):
    for r in rows:
        if r['cycle'] == cycle and r['on'] == on: return r['t']
    return None

def one(pid, rows, starts):
    t, y, D = series(rows)
    ly = np.log(y + EPS)
    hs, _ = fit(t, y, D, starts)
    s1, s2 = phase_start(rows, 1, 0), phase_start(rows, 2, 0)
    e2 = phase_start(rows, 3, 1)
    e2 = e2 if e2 is not None else t.max() + 1
    tr = t < s2
    te = (t >= s2) & (t < e2)
    fs, sols = fit(t[tr], y[tr], D, starts)
    m = simulate(fs.x, D, int(t.max()) + 1)
    err_model = np.log(m[t[te]] + EPS) - ly[te]
    c1 = (t >= s1) & (t < phase_start(rows, 2, 1))
    base = np.interp(t[te] - s2, t[c1] - s1, ly[c1])
    err_base = base - ly[te]
    # spread of forecasts among near-best fits (identifiability hint)
    good = [s for s in sols if s.cost <= fs.cost * 1.2 + 1e-9]
    ends = [float(simulate(s.x, D, int(t.max()) + 1)[t[te]][-1]) for s in good]
    return dict(pid=pid, n=len(t), n_test=int(te.sum()), hind=rmse(hs.fun), train=rmse(fs.fun),
                model=rmse(err_model), base=rmse(err_base), n_good=len(good),
                end_obs=float(y[te][-1]), end_lo=min(ends), end_hi=max(ends),
                hind_p=dict(zip(['rS','rR','dT','dD','K','S0','R0'], np.exp(hs.x).round(5).tolist())))

if __name__ == '__main__':
    P = patients()
    ids = sys.argv[1:] or sorted(P)
    starts = 12
    res = []
    for pid in ids:
        t0 = time.time()
        try:
            r = one(pid, P[pid], starts); r['sec'] = round(time.time() - t0, 1)
        except Exception as e:
            print('%s SKIP %s: %s' % (pid, type(e).__name__, e), flush=True); continue
        res.append(r)
        print('%s n=%d test=%d hind=%.3f train=%.3f model=%.3f base=%.3f %s good=%d end obs %.1f fits %.1f..%.1f  %.0fs' % (
            pid, r['n'], r['n_test'], r['hind'], r['train'], r['model'], r['base'],
            'MODEL' if r['model'] < r['base'] else 'base', r['n_good'], r['end_obs'], r['end_lo'], r['end_hi'], r['sec']), flush=True)
    wins = sum(r['model'] < r['base'] for r in res)
    print('SUMMARY patients %d model_beats_baseline %d (%.0f%%) median hind %.3f median model %.3f median base %.3f' % (
        len(res), wins, 100 * wins / len(res), np.median([r['hind'] for r in res]),
        np.median([r['model'] for r in res]), np.median([r['base'] for r in res])))
    json.dump(res, open('forecast_%s.json' % ('all' if not sys.argv[1:] else '_'.join(ids)), 'w'), indent=1)
