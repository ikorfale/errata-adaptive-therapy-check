"""CM-CANCER-Q01 v3 prototype (errata, 2026-10-01): score dosing rules by drug used at matched burden,
on parameter sets fitted to real patients (Bruchovsky 2006 IAS trial, 67 patients, every fit inside the
patient's 95% profile band from errata-adaptive-therapy-check).
Model (Strobl/Zhang form, PSA units): dS = rS S g (1 - dD D) - dT S,  dR = rR R g - dT R,  g = 1 - (S+R)/K.
Start at each patient's fitted day-0 state (N0 = S0 + R0), daily steps, rules decide at 7-day visits from burden only.
Progression: burden > 1.2 N0 (same line as bench v2). Horizon H = 1825 d (about the trial's follow-up).
Score of a rule on one parameter set:
  - if it progresses before H: FAIL, ranked below every survivor by TTP;
  - else: frontier dose / its dose, where frontier dose = least cumulative dose of any weekly modulation
    (bench v2 family, gain 2, targets 0.10..1.18) that survives H at the same or lower mean burden. >1 = less drug.
A rule beats the bar for a patient only if score > 1.01 on EVERY parameter set in that patient's band."""
import json, sys, numpy as np
import os; HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
from load import patients

H, VISIT, LINE = 1825, 7, 1.2
FAMILY = 'both'
TARGETS = np.round(np.arange(0.10, 1.181, 0.02), 3)

def cycle1(rows):
    """Durations of the patient's own first on and off phases (days)."""
    off = next(r['t'] for r in rows if r['on'] == 0)
    nxt = next((r['t'] for r in rows if r['on'] == 1 and r['t'] > off), None)
    return off, (nxt - off) if nxt else 365

def run(P, rules, c1):
    """P: (M,7) log-params. rules: list of names. Vectorised over parameter sets x rules."""
    M, R = len(P), len(rules)
    rS, rR, dT, dD, K, S0, R0 = [np.repeat(np.exp(P[:, i]), R) for i in range(7)]
    lS, lR = np.log(S0), np.log(R0); N0 = S0 + R0
    kind = np.tile(np.array(rules, dtype=object), M)
    D = np.ones(M * R); on = np.ones(M * R, bool)
    alive = np.ones(M * R, bool); ttp = np.full(M * R, H); dose = np.zeros(M * R); area = np.zeros(M * R)
    tgt = np.array([float(k[4:]) if k.startswith('mod_') else 0 for k in kind])
    is_mod = np.array([k.startswith('mod_') for k in kind]); is_a50 = kind == 'adaptive50'; is_c1 = kind == 'cycle1'
    is_mtd = kind == 'mtd'; is_none = kind == 'none'
    is_hys = np.array([k.startswith('hys_') for k in kind])
    lo = np.array([float(k.split('_')[1]) if k.startswith('hys_') else 0 for k in kind])
    hi = np.array([float(k.split('_')[2]) if k.startswith('hys_') else 0 for k in kind])
    on_d, off_d = c1
    for t in range(H):
        S, Rr = np.exp(lS), np.exp(lR); N = S + Rr
        if t % VISIT == 0:
            x = N / N0
            D = np.where(is_mod, np.clip(D + 2.0 * (x - tgt), 0, 1), D)
            on = np.where(is_a50 & on & (x <= 0.5), False, np.where(is_a50 & ~on & (x >= 1.0), True, on))
            D = np.where(is_a50, on.astype(float), D)
            D = np.where(is_c1, float((t % (on_d + off_d)) < on_d), D)
            on = np.where(is_hys & on & (x <= lo), False, np.where(is_hys & ~on & (x >= hi), True, on))
            D = np.where(is_hys, on.astype(float), D)
            D = np.where(is_mtd, 1.0, np.where(is_none, 0.0, D))
        prog = alive & (N > LINE * N0)
        ttp[prog] = t; alive &= ~prog
        dose += D * alive; area += (N / N0) * alive
        g = 1 - N / K
        lS = np.maximum(lS + rS * g * (1 - dD * D) - dT, -40.); lR = np.maximum(lR + rR * g - dT, -40.)
    mb = area / np.maximum(ttp, 1)
    return ttp.reshape(M, R), dose.reshape(M, R), mb.reshape(M, R)

def score_patient(pid, rows, prof):
    band = [q for q in prof['pts'] if q['ok']]
    P = np.array([q['p'] for q in band])
    mods = ['mod_%.2f' % x for x in TARGETS]
    hys = ['hys_%.2f_%.2f' % (a, b) for a in np.arange(0.1, 1.01, 0.1) for b in np.arange(0.2, 1.16, 0.1) if b > a + 0.05]
    mods = mods + hys if FAMILY == 'both' else mods
    rules = ['mtd', 'adaptive50', 'cycle1'] + mods + ['none']
    ttp, dose, mb = run(P, rules, cycle1(rows))
    nm = len(mods); out = {}
    for j, name in enumerate(rules[:3] + ['mod_1.00', 'mod_0.50']):
        j = rules.index(name); sc = []
        for m in range(len(P)):
            if ttp[m, j] < H: sc.append(('FAIL', int(ttp[m, j]))); continue
            surv = [(mb[m, k], dose[m, k]) for k in range(3, 3 + nm) if ttp[m, k] >= H and mb[m, k] <= mb[m, j] + 1e-9]
            sc.append(('OK', round(min(d for _, d in surv) / max(dose[m, j], 1e-9), 3)) if surv else ('BELOW', None))
        out[name] = dict(scores=sc, mean_burden=[round(float(v), 3) for v in mb[:, j]],
                         dose=[round(float(v), 1) for v in dose[:, j]], ttp=[int(v) for v in ttp[:, j]])
    nmod = len(TARGETS)
    allsurv = bool((ttp[:, 3:3 + nmod] >= H).all())          # v2's problem: every containment rule survives -> TTP can't rank
    anysurv = bool((ttp[:, 3:] >= H).any(axis=1).all())
    fr = []
    for m in range(len(P)):  # does adding on/off containment move the bar? least surviving dose at each modulation's burden
        for k in range(3, 3 + nmod):
            if ttp[m, k] < H: continue
            a = min(dose[m, q] for q in range(3, 3 + nmod) if ttp[m, q] >= H and mb[m, q] <= mb[m, k] + 1e-9)
            b = min(dose[m, q] for q in range(3, ttp.shape[1] - 1) if ttp[m, q] >= H and mb[m, q] <= mb[m, k] + 1e-9)
            fr.append(round(float(a / max(b, 1e-9)), 3))
    mtd_fail = [int(v) for v in ttp[:, 0]]; none_ttp = [int(v) for v in ttp[:, -1]]
    return dict(pid=pid, n_band=len(P), frontier_gain=fr, every_mod_survives=allsurv, some_mod_survives_all_sets=anysurv, mtd_ttp=mtd_fail, none_ttp=none_ttp, rules=out)

if __name__ == '__main__':
    prof = {r['pid']: r for r in json.load(open(os.path.join(HERE, '..', 'profile_all.json')))}
    Pt = patients(); ids = sys.argv[1:] or sorted(prof)
    res = [score_patient(pid, Pt[pid], prof[pid]) for pid in ids]
    json.dump(res, open(os.path.join(HERE, 'v3_all.json') if not sys.argv[1:] else '/tmp/v3_test.json', 'w'), indent=0)
    for r in res:
        print(r['pid'], r['n_band'], 'MTD ttp', r['mtd_ttp'], 'allmod', r['every_mod_survives'],
              ' '.join('%s:%s' % (k, v['scores']) for k, v in r['rules'].items()))
