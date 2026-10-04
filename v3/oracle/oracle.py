"""Oracle headroom (zenith-claude 71591). For each live fit of the 14 eligible patients: the model is known, the schedule
is open loop (one dose in [0, 1] per week, 261 weeks), and must keep the tumour below the line for H days at mean burden
<= B. oracle_drug(B) = least total dose found; frontier_drug(B) = least dose among contained family rules with mean
burden <= B (on-off 0.02 grid, integral/proportional setpoints on a 0.01 grid, as in maxmin064.py).
oracle_ratio = frontier_drug / oracle_drug. The optimiser starts from the frontier rule's own weekly doses, so the
ratio is >= 1 by construction; the optimiser is local (SLSQP), so the oracle number is an upper bound on the true
optimum dose and the ratio a LOWER bound on the true headroom.  B = mean burden of the frontier rule at three levels:
the 25th, 50th and 75th percentile of mean burden among contained family rules on that fit."""
import json, os, sys, numpy as np
from scipy.optimize import minimize
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "v31cens"))
from cm_cancer_q01_v3 import DATA, simulate, onoff, H, VISIT, LINE
from cm_cancer_q01_v31 import integral, proportional
W = (H + VISIT - 1) // VISIT
def sim_batch(p, U):
    """U: (b, W) weekly doses. Same Euler steps as simulate(); returns daily dose sum, mean burden, max N/N0."""
    rS, rR, dT, dD, K, S0, R0 = np.exp(p); N0 = S0 + R0; b = U.shape[0]
    lS = np.full(b, np.log(S0)); lR = np.full(b, np.log(R0)); area = np.zeros(b); top = np.zeros(b)
    for t in range(H):
        N = np.exp(lS) + np.exp(lR); x = N / N0; area += x; top = np.maximum(top, x)
        D = U[:, t // VISIT]; g = 1.0 - N / K
        lS = np.maximum(lS + rS * g * (1.0 - dD * D) - dT, -40.0); lR = np.maximum(lR + rR * g - dT, -40.0)
    days = np.full(W, VISIT); days[-1] = H - VISIT * (W - 1)
    S, R = np.exp(lS), np.exp(lR)
    sim_batch.end = ((S + R) / N0, S / (S + R))                         # end state: burden and sensitive share at H
    return U @ days, area / H, top
def weekly(p, rule):
    rS, rR, dT, dD, K, S0, R0 = np.exp(p); N0 = S0 + R0
    lS, lR, state, u = np.log(S0), np.log(R0), {}, []
    for t in range(H):
        N = np.exp(lS) + np.exp(lR)
        if t % VISIT == 0: D = rule(t, N / N0, state); u.append(D)
        g = 1.0 - N / K
        lS = max(lS + rS * g * (1.0 - dD * D) - dT, -40.0); lR = max(lR + rR * g - dT, -40.0)
    return np.array(u)
def oracle(p, B, u0, line=LINE - 1e-3, term=None):
    """term=(xH, fH): the schedule must also end no worse than the frontier rule (burden at H <= xH and sensitive
    share at H >= fH), so it cannot coast on the end of the horizon."""
    days = np.full(W, float(VISIT)); days[-1] = H - VISIT * (W - 1)
    cache = {}; best = [None, None]                                     # best feasible (dose, u) seen, start included
    def get(u):
        k = u.tobytes()
        if k not in cache:
            u = np.clip(u, 0, 1); h = np.where(u < 0.5, 1e-4, -1e-4)       # step into the box, never clipped away
            d, mb, top = sim_batch(p, np.vstack([u, u + np.diag(h)])); xe, fe = sim_batch.end
            c = [B - mb, line - top] + ([term[0] + 1e-9 - xe, fe - term[1] + 1e-9] if term else [])
            f0 = np.array([d[0]] + [ci[0] for ci in c])
            if min(f0[1:]) >= 0 and (best[0] is None or d[0] < best[0]): best[0], best[1] = d[0], u.copy()
            J = np.vstack([days] + [(ci[1:] - ci[0]) / h for ci in c])
            cache.clear(); cache[k] = (f0, J)
        return cache[k]
    r = minimize(lambda u: get(u)[0][0], u0, jac=lambda u: get(u)[1][0], method="SLSQP", bounds=[(0, 1)] * W,
                 constraints=[{"type": "ineq", "fun": lambda u: get(u)[0][1:], "jac": lambda u: get(u)[1][1:]}],
                 options={"maxiter": 100, "ftol": 1e-6})
    get(np.asarray(r.x, float))
    return (None if best[0] is None else float(best[0])), best[1], r.nit
TERM = os.environ.get("TERM_END", "1") == "1"
if __name__ == "__main__":
    prof = {r["pid"]: r for r in json.load(open(os.path.join(DATA, "errata_profile_all_ad66924.json")))}
    none = lambda t, x, s: 0.0
    grid = np.round(np.arange(0.30, 1.19, 0.02), 2); sp = np.round(np.arange(0.30, 1.21, 0.01), 2)
    fam = [onoff(lo, hi) for lo in grid for hi in grid if hi > lo + 0.01] + \
          [integral(t, g) for g in (1.0, 2.0, 5.0, 10.0, 15.0, 30.0) for t in sp] + \
          [proportional(t, g) for g in (3.0, 5.0, 10.0, 15.0, 30.0, 60.0) for t in sp]
    pids = sys.argv[1:] or sorted(pid for pid, r in prof.items()
            if sum(1 for q in r["pts"] if q["ok"] and simulate(np.array(q["p"]), none)[0] < H) >= 2)
    skip = set(os.environ.get("SKIP", "").split(","))
    for pid in [x for x in pids if x not in skip]:
        fits = [np.array(q["p"]) for q in prof[pid]["pts"] if q["ok"]]
        for k, p in enumerate(fits):
            if simulate(p, none)[0] >= H: continue
            runs = [(simulate(p, r), r) for r in fam]; ok = [(mb, d, r) for (ttp, d, mb), r in runs if ttp >= H]
            mbs = np.array([o[0] for o in ok])
            if not ok:                                                     # no rule in the family survives the horizon
                print(f"{pid} fit {k} no feasible family rule (all progress before H) term {int(TERM)}", flush=True); continue
            for q in (25, 50, 75):
                B = float(np.percentile(mbs, q)); cand = [o for o in ok if o[0] <= B + 1e-9]
                fd = min(o[1] for o in cand); rule = [o[2] for o in cand if o[1] == fd][0]
                u0 = weekly(p, rule); sim_batch(p, u0[None, :]); term = (sim_batch.end[0][0], sim_batch.end[1][0])
                od, u, nit = oracle(p, B, u0, term=term if TERM else None)
                if od is None:                                             # start itself infeasible in sim_batch margins
                    d0, mb0, top0 = sim_batch(p, u0[None, :])
                    print(f"{pid} fit {k} q{q} B {B:.4f} frontier {fd:.1f} oracle None start_mb {mb0[0]:.6f} "
                          f"start_top {top0[0]:.6f} line {LINE} term {int(TERM)}", flush=True); continue
                ttp, vd, vmb = simulate(p, lambda t, x, s, u=u: float(u[t // VISIT]))     # independent replay check
                assert ttp >= H and vmb <= B + 1e-9 and abs(vd - od) < 1e-6, (ttp, vd, od, vmb, B)
                print(f"{pid} fit {k} q{q} B {B:.4f} frontier {fd:.1f} oracle {od if od is None else round(od, 1)} "
                      f"ratio {None if od is None else round(fd / od, 4)} iters {nit} term {int(TERM)}", flush=True)
