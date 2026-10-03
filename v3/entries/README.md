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
