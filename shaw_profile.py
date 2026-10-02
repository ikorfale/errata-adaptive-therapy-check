"""Holdout cohort: run the same profile fit (profile.py, unchanged) on the 18 Shaw et al. patients of the
same public archive (dataTanaka.zip). Ids are the 4-digit file numbers. Writes profile_shaw.json."""
import glob, json, os, sys, time
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here); os.chdir(here)
import load, profile
D = os.path.join(here, 'data/bruchovsky/dataTanaka/Shaw_et_al')
P = {os.path.basename(p)[7:-4]: load.load(p) for p in sorted(glob.glob(D + '/patient*.txt'))}
out = []
for pid in sorted(P):
    t0 = time.time()
    try: r = profile.one(pid, P[pid])
    except Exception as e: print(pid, 'SKIP', type(e).__name__, e, flush=True); continue
    out.append(r)
    print('%s ok %d/9  ttp %d..%d days  forecast rmse %.2f..%.2f  %.0fs' % (pid, r['n_ok'], r['ttp_lo'], r['ttp_hi'], r['fc_lo'], r['fc_hi'], time.time() - t0), flush=True)
    json.dump(out, open('profile_shaw.json', 'w'), indent=1)
print('done', len(out), 'of', len(P))
