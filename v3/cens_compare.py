"""Compare live band sets per patient: profile_all.json (PSA floors as exact) vs profile_cens.json (floors censored)."""
import json, os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import v31
def live(fn):
    out = {}
    for r in json.load(open(os.path.join(HERE, '..', fn))):
        P = np.array([q['p'] for q in r['pts'] if q['ok']])
        ttp, _, _ = v31.run(P, [('none', 0., 0.)])
        out[r['pid']] = (len(P), int((ttp[:, 0] < v31.H).sum()))
    return out
a, b = live('profile_all.json'), live('profile_cens.json')
el = lambda d: sorted(k for k, (n, l) in d.items() if l >= 2)
for k in sorted(set(a) | set(b)):
    if a.get(k) != b.get(k): print(k, 'exact', a.get(k), 'censored', b.get(k))
print('eligible exact', len(el(a)), 'censored', len(el(b)), 'gained', sorted(set(el(b)) - set(el(a))), 'lost', sorted(set(el(a)) - set(el(b))))
