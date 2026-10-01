"""CM-CANCER-Q01 v3.1 check (errata, 2026-10-01): aria's two holes in v3.
1. Frontier grid gap: add integral modulation with gains 2/5/15 and proportional setpoint rules
   D = clamp(g (x - target)) with targets every 0.002 up to the line; does agentcue's setpoint rule
   (CM-CANCER-103, D = clamp(15 (x - 1.173))) still beat the frontier anywhere?
2. Live sets: a set is live if the untreated tumour crosses the line within H; count patients by live sets."""
import json, sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
from load import patients
H, VISIT, LINE = 1825, 7, 1.2

def run(P, rules):
    """rules: list of (kind, gain, target); kind 'int' (D += g(x-t)), 'prop' (D = g(x-t)), 'none'. Clamped to [0,1]."""
    M, R = len(P), len(rules)
    rS, rR, dT, dD, K, S0, R0 = [np.repeat(np.exp(P[:, i]), R) for i in range(7)]
    lS, lR = np.log(S0), np.log(R0); N0 = S0 + R0
    kind = np.tile(np.array([r[0] for r in rules]), M); g = np.tile(np.array([r[1] for r in rules], float), M)
    tg = np.tile(np.array([r[2] for r in rules], float), M)
    D = np.where(kind == 'none', 0., 1.); alive = np.ones(M * R, bool); ttp = np.full(M * R, H)
    dose = np.zeros(M * R); area = np.zeros(M * R)
    for t in range(H):
        N = np.exp(lS) + np.exp(lR)
        if t % VISIT == 0:
            x = N / N0
            D = np.where(kind == 'int', np.clip(D + g * (x - tg), 0, 1), D)
            D = np.where(kind == 'prop', np.clip(g * (x - tg), 0, 1), D)
        prog = alive & (N > LINE * N0); ttp[prog] = t; alive &= ~prog
        dose += D * alive; area += (N / N0) * alive
        gg = 1 - N / K
        lS = np.maximum(lS + rS * gg * (1 - dD * D) - dT, -40.); lR = np.maximum(lR + rR * gg - dT, -40.)
    return ttp.reshape(M, R), dose.reshape(M, R), (area / np.maximum(ttp, 1)).reshape(M, R)

def frontier(ttp, dose, mb, m, cols, j):
    s = [dose[m, k] for k in cols if ttp[m, k] >= H and mb[m, k] <= mb[m, j] + 1e-9]
    return min(s) if s else None

if __name__ == '__main__':
    prof = {r['pid']: r for r in json.load(open(os.path.join(HERE, '..', 'profile_all.json')))}
    Pt = patients(); ids = sys.argv[1:] or sorted(prof)
    coarse = [('int', 2., t) for t in np.round(np.arange(0.10, 1.181, 0.02), 3)]
    fine_t = np.round(np.arange(1.000, 1.1995, 0.002), 3)
    fine = [('int', gg, t) for gg in (2., 5., 15.) for t in fine_t] + [('prop', gg, t) for gg in (5., 15., 30.) for t in fine_t]
    entry = [('prop', 15., 1.173)]
    rules = entry + coarse + fine + [('none', 0., 0.)]
    ci = list(range(1, 1 + len(coarse))); fi = list(range(1, len(rules) - 1))
    out = []
    for pid in ids:
        P = np.array([q['p'] for q in prof[pid]['pts'] if q['ok']])
        ttp, dose, mb = run(P, rules)
        live = [bool(v < H) for v in ttp[:, -1]]
        sets = []
        for m in range(len(P)):
            if ttp[m, 0] < H: sets.append(dict(live=live[m], status='FAIL', ttp=int(ttp[m, 0]))); continue
            a = frontier(ttp, dose, mb, m, ci, 0); b = frontier(ttp, dose, mb, m, fi, 0)
            sets.append(dict(live=live[m], status='OK', mb=round(float(mb[m, 0]), 4), dose=round(float(dose[m, 0]), 1),
                             v3=round(a / max(dose[m, 0], 1e-9), 3) if a else None, refined=round(b / max(dose[m, 0], 1e-9), 3) if b else None))
        out.append(dict(pid=pid, n=len(P), n_live=sum(live), sets=sets))
        print(pid, len(P), 'live', sum(live), sets, flush=True)
    json.dump(out, open('/tmp/v31.json' if sys.argv[1:] else os.path.join(HERE, 'v31_all.json'), 'w'), indent=0)
