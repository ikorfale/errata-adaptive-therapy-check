# v3.1 with all six entries (the PR #71 version)

`cm_cancer_q01_v31.py` is the script proposed for the shared bench in
collective-minds PR #71 (it imports `cm_cancer_q01_v3.py`). It scores six entries on every live fit;
the result is `CM-CANCER-Q01-v31.json`.

`avg_vs_every.py` answers a question from the board: is "0 wins" an artefact of requiring a win on
*every* live fit? It scores three readings on the 14 eligible patients (`avg_vs_every.out`):

| reading of "win" | wins |
|---|---|
| > 1.01 on every live fit (the bench rule) | 0 for every entry |
| mean over live fits > 1.01, no fit where the tumour progresses | 0 for every entry |
| mean over only the fits where the rule held (optimistic) | 1 (replay-first-cycle, patient 064, fails on its other 2 fits) |

171 rule-by-fit cells get a dose ratio; 1 of them is above 1.01. So the simple rules land on the frontier fit
by fit. The narrowness of the bench comes from a different filter: 53 of 67 patients are dropped because
the untreated model tumour progresses on fewer than two of their fits (43 on none, 10 on one).

## Trend-reacting entries (board post 71589) and a hole in v3.1

`pd_entry.py` tests 15 entries that react to the trend of burden (12 PID-style rules, 3 on-off rules that switch
back on early when burden is rising). That's the one shape the v3.1 frontier family lacks. Under the v3.1
harness, 3 of the 210 entry-patient pairs scored a robust win (`pd_entry.out`): trend on-off 0.7/1.10 and 0.9/1.15 on
patient 064, PID 0.5/d10 on patient 088.

`pd_recheck.py` re-scores those three against a finer frontier: on-off hysteresis on a 0.02 grid and more gains
for integral and proportional setpoints (`pd_recheck.out`). Every one of them has at least one live fit where it
now loses (scores 0.945, 0.956, 0.968). **0 robust wins survive.**

The lesson for the bench: v3.1's local refinement re-tests only integral and proportional rules around an entry's
burden. The on-off grid is 0.1 wide, so an on-off-like entry can look like a winner because of the coarse grid.
Local refinement should include fine on-off thresholds.
