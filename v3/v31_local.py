"""v3.1 second refinement (errata, 2026-10-02): re-score an OK entry against a grid refined around its own mean burden
(targets mb-0.01..mb+0.01 every 0.0002, integral gains 2/5/15, proportional gains 5/10/15/20/30). Usage: v31_local.py PID SET"""
import json, sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..'))
import v31
prof = {r['pid']: r for r in json.load(open(os.path.join(HERE, '..', 'profile_all.json')))}
pid, m = sys.argv[1], int(sys.argv[2])
P = np.array([q['p'] for q in prof[pid]['pts'] if q['ok']])[m:m + 1]
entry = ('prop', 15., 0.88 / 0.75)
ttp, dose, mb = v31.run(P, [entry])
ts = np.round(np.arange(mb[0, 0] - 0.01, mb[0, 0] + 0.0101, 0.0002), 5)
grid = [('int', g, t) for g in (2., 5., 15.) for t in ts] + [('prop', g, t) for g in (5., 10., 15., 20., 30.) for t in ts]
ttp, dose, mb = v31.run(P, [entry] + grid)
b = v31.frontier(ttp, dose, mb, 0, list(range(1, len(grid) + 1)), 0)
print(pid, m, 'entry mb', round(float(mb[0, 0]), 4), 'dose', round(float(dose[0, 0]), 1), 'local score', round(b / dose[0, 0], 4) if b else None)
