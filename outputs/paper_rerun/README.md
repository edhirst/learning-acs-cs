# Paper rerun: run records

Records of the full rerun of the paper's experiments (25–26 Oct 2026), made with the code of commit 12796ba. `experiments/run_all.py` ran 72 jobs; four further S³×S¹ jobs were run by hand with the same code (`scheduler_manual.log`). All 76 jobs finished and none failed. The test suite (63 tests) passes on this code.

| Path | Contents |
|---|---|
| `runs/<job>/` | Per-run records: `result.json` (final metrics), `progress.jsonl` (training history), `done.json` or `certificate.json` (analysis and reconstruction jobs), evaluation grids (`*.npz`) |
| `logs/<job>.txt` | Console log of each job |
| `scheduler.log`, `scheduler_manual.log`, `scheduler_status.json` | Job order, durations and exit codes |
| `tables/*.tex`, `tables/numbers.tex` | LaTeX table bodies and number macros |
| `figures/*.pdf` | Figures |
| `numbers.txt`, `numbers.json` | Every quoted number with the run record it comes from (`experiments/paper_numbers.py`) |
| `code_state/README.txt` | Campaign log: code changes and reruns made while the campaign ran |

Not included:
- the model checkpoints (`*.pt`, about 225 MB), which the repository ignores;
- the frozen copy of the code used during the campaign. It is identical to commit 12796ba except for a diagnostic in `extra/ambient_supervised_hopf.py` that records the base's Hopf invariant.

To regenerate the tables, figures and number index from these records:

```
PAPER_OUT=outputs/paper_rerun .venv/bin/python experiments/make_figures.py
.venv/bin/python experiments/paper_numbers.py --runs outputs/paper_rerun/runs \
    --json outputs/paper_rerun/numbers.json > outputs/paper_rerun/numbers.txt
```

To rerun the campaign into a new folder (CPU, double precision): `PAPER_OUT=<folder> .venv/bin/python experiments/run_all.py --budget 9`. Summed over jobs this took about 114 hours of wall-clock time, most of it the four S⁵×S¹ supervised fits (about 19 h each).
