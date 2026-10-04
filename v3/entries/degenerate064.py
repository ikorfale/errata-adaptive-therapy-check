"""zenith-claude 72181: is the phase-robust best, trend on-off lo 0.98 / hi 1.18 (score exactly 1.0000), just a frontier
rule? Compare its weekly dose sequence with the plain on-off (0.98, 1.18), which is in the frontier family, per fit and phase."""
import numpy as np
from phase064 import simulate, fits, onoff, trend_onoff, VISIT, H
for k, p in enumerate(fits):
    for ph in range(VISIT):
        a, b = [], []; ra = simulate(p, trend_onoff(0.98, 1.18), ph, a); rb = simulate(p, onoff(0.98, 1.18), ph, b)
        print(f"fit {k} phase {ph}: trend ttp {ra[0]} dose {ra[1]:.1f} mb {ra[2]:.4f} | plain on-off ttp {rb[0]} dose {rb[1]:.1f} mb {rb[2]:.4f} | "
              f"identical weekly doses: {a == b} (weeks differing {sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))})")
