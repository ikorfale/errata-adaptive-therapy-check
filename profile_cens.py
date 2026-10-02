"""Sensitivity check (errata, 2026-10-02): the same profile fit as profile.py, but PSA values recorded at an assay
floor (exactly 0.1 or 0.02 ng/ml; Cancer 2007 doi:10.1002/cncr.22464 gives <0.02 as the lower reportable limit, and
0.1 piles up in 1996-97 rows) are treated as left-censored: a model value below the floor costs nothing.
Writes profile_cens.json. Compare live sets / eligible patients with profile_all.json."""
import json, sys, time
import numpy as np
import model, profile
from load import patients

def resid_cens(p, t, y, D):
    m = model.simulate(p, D, int(t.max()) + 1)[t]
    r = np.log(m + model.EPS) - np.log(y + model.EPS)
    floor = np.isclose(y, 0.1) | np.isclose(y, 0.02)
    return np.where(floor & (r < 0), 0.0, r)

profile.resid = resid_cens
if __name__ == '__main__':
    P = patients(); out = []
    for pid in sorted(P):
        t0 = time.time()
        try: r = profile.one(pid, P[pid])
        except Exception as e: print(pid, 'SKIP', type(e).__name__, e, flush=True); continue
        out.append(r)
        print('%s ok %d/9  %.0fs' % (pid, r['n_ok'], time.time() - t0), flush=True)
        json.dump(out, open('profile_cens.json', 'w'), indent=1)
