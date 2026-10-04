"""zenith-claude 72066 on patient 064: (1) how many DISTINCT trend on-off rules are in the 1830 widened pairs (identical
dose sequences on all live fits)? (2) score each pair as the minimum over the 7 visit phases (weekly visit on day
ph, ph+7, ...; drug on from day 0 until the first visit), frontier recomputed at the same phase. Same frontier family
and scoring as maxmin064.py (on-off top 1.18, setpoint gains as there). Phase 0 must reproduce maxmin064.out."""
import json, os, hashlib, numpy as np
import cm_cancer_q01_v3 as v3
from cm_cancer_q01_v3 import DATA, onoff, H, LINE, VISIT
from cm_cancer_q01_v31 import integral, proportional
from pd_entry import trend_onoff
def simulate(p, rule, ph=0, trace=None):
    rS, rR, dT, dD, K, S0, R0 = np.exp(p)
    lS, lR, N0, D, state, dose, area = np.log(S0), np.log(R0), S0 + R0, 1.0, {}, 0.0, 0.0
    for t in range(H):
        S, R = np.exp(lS), np.exp(lR); N = S + R
        if N > LINE * N0: return t, dose, area / max(t, 1)
        if t % VISIT == ph:
            D = rule(t, N / N0, state)
            if trace is not None: trace.append(D)
        dose += D; area += N / N0
        g = 1.0 - N / K
        lS = max(lS + rS * g * (1.0 - dD * D) - dT, -40.0); lR = max(lR + rR * g - dT, -40.0)
    return H, dose, area / H
prof = {r["pid"]: r for r in json.load(open(os.path.join(DATA, "errata_profile_all_ad66924.json")))}
none = lambda t, x, s: 0.0
fits = [np.array(q["p"]) for q in prof["064"]["pts"] if q["ok"]]
fits = [p for p in fits if simulate(p, none)[0] < H]
grid = np.round(np.arange(0.30, 1.51, 0.02), 2); pairs = [(lo, hi) for lo in grid for hi in grid if hi > lo + 0.01]
gf = np.round(np.arange(0.30, 1.19, 0.02), 2); sp = np.round(np.arange(0.30, 1.21, 0.01), 2)
famf = lambda: [onoff(lo, hi) for lo in gf for hi in gf if hi > lo + 0.01] + [integral(t, g) for g in (1.0, 2.0, 5.0, 10.0, 15.0, 30.0) for t in sp] + \
      [proportional(t, g) for g in (3.0, 5.0, 10.0, 15.0, 30.0, 60.0) for t in sp]
print("live fits", len(fits), "pairs", len(pairs), "VISIT", VISIT, flush=True)
S = np.zeros((len(pairs), len(fits), VISIT)); sig = [[None] * len(fits) for _ in pairs]
for ph in range(VISIT):
    for k, p in enumerate(fits):
        F = np.array([(mb, d) for ttp, d, mb in (simulate(p, r, ph) for r in famf()) if ttp >= H])
        for i, (lo, hi) in enumerate(pairs):
            tr = [] if ph == 0 else None; ttp, d, mb = simulate(p, trend_onoff(lo, hi), ph, tr)
            if ph == 0: sig[i][k] = hashlib.sha1(np.array(tr).tobytes()).hexdigest()
            if ttp < H: S[i, k, ph] = 0.0; continue
            ok = F[F[:, 0] <= mb + 1e-9]; S[i, k, ph] = ok[:, 1].min() / d if len(ok) else np.inf
        print("phase", ph, "fit", k, "done", flush=True)
np.save("phase064.npy", S)
print("distinct trend on-off rules (identical switch sequence on all live fits, phase 0):", len(set(map(tuple, sig))), "of", len(pairs))
for hi_max in (1.18, 1.50):
    idx = [i for i, (lo, hi) in enumerate(pairs) if hi <= hi_max + 1e-9]
    print(f"--- pairs with hi <= {hi_max}: {len(idx)}, distinct {len(set(tuple(sig[i]) for i in idx))}")
    m0 = S[idx, :, 0].min(1); i0 = idx[int(np.argmax(m0))]
    print(f"phase 0 max-min: lo {pairs[i0][0]:.2f} hi {pairs[i0][1]:.2f} {m0.max():.4f}")
    mp = S[idx].min((1, 2)); ip = idx[int(np.argmax(mp))]
    print(f"min over fits and phases, best: lo {pairs[ip][0]:.2f} hi {pairs[ip][1]:.2f} {mp.max():.4f}; the phase-0 winner scores {S[i0].min():.4f} (per phase {np.round(S[i0].min(0), 4).tolist()})")
    for k in range(len(fits)):
        b0 = idx[int(np.argmax(S[idx, k, 0]))]; bp = idx[int(np.argmax(S[idx, k, :].min(1)))]
        print(f"fit {k}: best at phase 0 lo {pairs[b0][0]:.2f} hi {pairs[b0][1]:.2f} {S[b0, k, 0]:.4f} -> its min over phases {S[b0, k].min():.4f} "
              f"(per phase {np.round(S[b0, k], 4).tolist()}); best by min over phases: lo {pairs[bp][0]:.2f} hi {pairs[bp][1]:.2f} {S[bp, k].min():.4f}")
