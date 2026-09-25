# Ambient ACS campaign: implementation and Lovelace handoff

The next experiment is **ACS only**. No CS job, hold or dependency is needed.
This implements the first phase of [the ambient ACS plan](s4xt2_ambient_acs_plan.md).
Stop after the first four pilots and review their results before spending the
larger multi-seed/long-run budget. Existence of an ACS does not guarantee that
Adam finds it in the chosen budget.

## What is implemented

- The same nine ambient inputs, 81 raw matrix outputs and `K = P M P` as before.
  Random-start, residual-only training; no oracle labels or analytic warm start.
- Generic `all2+trilin` features: 55 degree-at-most-two features plus 20
  three-factor products on S4×S1×S1. Coefficient-level representability test.
- A genuinely polynomial-only first phase; zero-output MLP activation preserves
  its field. Existing checkpoint interpretation and the old CLI remain available.
- Uniform, adaptive/replay, and uniform-plus-hard-tail objectives; independent
  training, replay, model-selection and final-audit random streams.
- Worst-residual-first checkpoint selection on common validation/challenge points;
  saved polynomial candidate; atomic phase-aware restart including Adam/scheduler,
  RNG, replay and actual counters. L-BFGS is chunked, gated and rolled back when
  validation RMS or maximum worsens. Completed outputs cannot be silently reused.
- YAML configurations and a separate checkpoint-hashed numerical auditor, including
  structured points, independent ascent, tangent spectra, corrected residuals,
  tangency, nonsymmetric SciPy square-root comparison, derivatives and chart overlap.
  It never calls finite sampling a rigorous certificate.

The new runner is selected by `python extra/campaign_unit.py --config ...`.
Legacy argument-only ACS/CS commands keep their original training loops. Use the
new `_v2` submission script, not the older `submit_S4xS1xS1_ACS.sh`.

## Get the files onto Lovelace

The changes are local until you commit and push them. The entire `scratch/`
directory remains gitignored. From this local checkout, inspect and commit the
intended changes, including the new files (ordinary `git commit -a` misses them):

```bash
git diff --check
git status --short
git add src/geometry/ambient.py extra/campaign_unit.py \
  extra/campaign_summary.py extra/campaign_recertify.py \
  extra/acs_training.py extra/acs_local_checks.py \
  extra/validate_ambient_acs.py extra/test_acs_training.py \
  configs/ambient/s4xs1xs1_acs docs/s4xt2_ambient_acs_plan.md \
  docs/s4xt2_research_plan.md docs/s4xt2_acs_hpc.md
git commit -m "Add reproducible ambient ACS pilot campaign"
git push origin no_overlaps
```

Copy the two ignored scripts separately, using your working SSH hostname/alias:

```bash
scp scratch/submit_S4xS1xS1_ACS_v2.sh scratch/submit_S4xS1xS1_ACS_audit.sh \
  ehirst@lovelace:/home/lovelace/proj/proj925/ehirst/nijenhuis/
```

On Lovelace:

```bash
cd /home/lovelace/proj/proj925/ehirst/nijenhuis/github
git branch --show-current
git status --short
test "$(git branch --show-current)" = no_overlaps && git pull --ff-only origin no_overlaps
test -f extra/acs_training.py
test -f configs/ambient/s4xs1xs1_acs/pilot_trilin_tail.yaml
cd ..
```

The branch should be `no_overlaps`; stop and resolve a different branch or local
conflicts rather than resetting your checkout. Pulling a branch into a different
checked-out branch is not the same as switching to that branch.

## 1. Submit a short smoke test

Both scripts request the previously working `serial`, `mem=10G`, one-thread setup,
with no `nodes`/`ppn` directive. They additionally request four hours. This limit
has **not** been verified against your live queue; inspect `qstat -Qf serial` or
your site's queue documentation if rejected. If changing the PBS wall limit,
also set `WALL_MINUTES` to no more than 85% of its minutes.

```bash
qsub -v RECIPE=smoke,SEED=7,WALL_MINUTES=10 submit_S4xS1xS1_ACS_v2.sh
qstat -u "$USER"
```

The smoke run is 100 polynomial + 100 refinement steps. It is expected to finish
quickly and is **not** expected to learn an accurate ACS. Confirm:

- The job log under `logs/S4T2_acs_smoke_*.log` ends with
  `completed_pending_independent_audit`.
- Its output directory contains `state.pt`, `best_poly.pt`, `best.pt`, `final.pt`,
  `manifest.json`, `progress.jsonl`, and `result.json`.
- There is no exception or missing module/config error.

The output directory is printed near the beginning of the job log. By default:

```text
outputs/s4xs1xs1_acs/pilots/<recipe>/s<seed>/<full-PBS-job-id>/
```

## 2. First bounded HPC pilots: seed 7 only

After the smoke succeeds, submit these four independent jobs. They do not depend
on each other and none should be held for a CS job:

```bash
qsub -v RECIPE=pilot_all2_uniform,SEED=7 submit_S4xS1xS1_ACS_v2.sh
qsub -v RECIPE=pilot_trilin_uniform,SEED=7 submit_S4xS1xS1_ACS_v2.sh
qsub -v RECIPE=pilot_trilin_adaptive,SEED=7 submit_S4xS1xS1_ACS_v2.sh
qsub -v RECIPE=pilot_trilin_tail,SEED=7 submit_S4xS1xS1_ACS_v2.sh
```

| Recipe | Features | Collocation | Loss |
| --- | --- | --- | --- |
| A0 `pilot_all2_uniform` | all degree ≤2 | uniform | mean squared residual |
| A1 `pilot_trilin_uniform` | degree ≤2 + three-factor | uniform | mean squared residual |
| A2 `pilot_trilin_adaptive` | same as A1 | ½ uniform, ¼ worst-pool, ¼ replay | mixed mean |
| A3 `pilot_trilin_tail` | same as A1 | same as A2 | uniform mean + ramped hard-tail mean |

Each allows 5,000 polynomial and 5,000 refinement Adam updates, plus at most 200
L-BFGS iterations **only if** a candidate's maximum found is at most 0.5.
Refinement is skipped if the polynomial alone already reaches the screening
targets. MLP widths are `[256,256,256,256]`. Learning rates are 0.005 polynomial,
0.001 refinement, cosine decay to 1e-5, in float64 on one CPU thread.

The three-hour application guard is per allocation, not accumulated over resumes.
In-progress checkpoints are saved every 250 updates or five minutes, and on
cooperative SIGTERM/SIGINT/application timeout. A hard SIGKILL can recover only
the last completed atomic checkpoint. A partial run exits **75**, has status
`incomplete`, and does not create `final.pt`. Exit zero means the scheduled
experiment completed; it says nothing by itself about ACS accuracy.

To resume a partial run, use the **same recipe, seed and exact existing output
directory** (replace the example job directory with the path printed in your log):

```bash
qsub -v RECIPE=pilot_trilin_tail,SEED=7,RESUME=1,OUT=/home/lovelace/proj/proj925/ehirst/nijenhuis/outputs/s4xs1xs1_acs/pilots/pilot_trilin_tail/s7/ACTUAL_JOB_ID \
  submit_S4xS1xS1_ACS_v2.sh
```

Changing the step budget or other numerical settings is a **continuation**, not
an exact resume. The CLI supports `--init-from` into a new directory, requiring
`poly_steps: 0`, matching hidden widths and features (or explicit
`--migrate-features` for a zero-coefficient basis extension). It restarts the
optimizer/scheduler and records the source checksum. We will configure the
seed-10 legacy comparisons after inspecting that checkpoint's actual metadata;
do not silently assume which ladder feature class it saved.

## 3. Bring back results; pause before the long campaign

With the `nijenhuis` conda environment active:

```bash
cd /home/lovelace/proj/proj925/ehirst/nijenhuis
python3 github/extra/campaign_summary.py outputs/s4xs1xs1_acs/pilots --all
```

Send that output and the four pilots' `result.json` and `progress.jsonl` files.
We need the **raw RMS and worst residual**, phase counters, loss curves and wall
times. A low mean with maximum near 1.4–1.8 is not success. No candidate from the
local smoke runs passed the target ACS accuracy gate.

After review, run the useful recipes on seeds 8 and 9, choose one or two based on
tails and time rather than mean alone, then prepare the 20k+20k runs on fresh
seeds 20–22. Those jobs are deliberately **not** submitted/configured blindly
now. Normalized features, matrix-free Gauss–Newton/LM, and other escalation
options remain deferred until the pilot evidence justifies them.

## Independent audit, when a promising candidate exists

Do not spend the final held-out audit repeatedly during tuning. Use `--quick`
with a separate diagnostic seed for preliminary checks. Only promote candidates
with raw RMS ≤1e-4 and maximum found ≤1e-2 to the full audit:

```bash
qsub -v CHECKPOINT=/absolute/path/to/candidate/final.pt \
  submit_S4xS1xS1_ACS_audit.sh
```

The full audit uses 100,000 fresh uniform points, structured pole/circle grids,
and three independent 1,024-start, 200-step adversarial searches with 32,768-point
initial pools each. It checks the corrected field separately (identity ≤1e-10,
tangency ≤1e-12), derivatives and independent matrix roots. The checkpoint
checksum must remain unchanged. The report records `numerically_validated`,
never a rigorous global certificate. A quick report **cannot** set this flag.
An audit killed by the scheduler produces no completed report; rerun with a
sufficient allocation. This audit is ACS-only; Nijenhuis/CS work remains deferred.

## Local verification

Verified on 2026-09-25 with Python 3.10 / PyTorch 2.5.1, CPU float64:

- Focused tests: **33 passed** (including the eight existing product tests).
- All three test files: **42 passed, ten pre-existing dtype failures**.
- Thirteen bounded actual runs completed: six cross-manifold controls,
  three target smoke seeds, and four full-width pilot wiring checks.
- S2 linear control: raw RMS 1.39e-4, maximum found 4.08e-4.
- Torus bilinear control: raw RMS 1.98e-3, maximum found 6.74e-3.
- Target smoke seeds 7/8/9: raw maxima 1.718 / 1.676 / 2.015, **not success**.
- A preliminary torus `all2` control stalled near 1.415; the result is retained
  locally. Simplifying to the existing `bilin` class solved the positive-control
  check, supporting the need for the target feature/sampling comparisons.
- Shell syntax and mocked module/conda/Python execution passed, including path
  quoting and propagation of exit 75. Live PBS submission is not locally tested.
- A quick independent audit rejected the target smoke, as expected; the exact
  test fixture passed numerical checks but correctly could not receive full-audit
  acceptance. No full-scale held-out audit has yet been spent on a learned field.

```bash
python -m pytest -q extra/test_acs_training.py extra/test_ambient_products.py
python -m extra.acs_local_checks --out scratch/acs_local_checks/NEW_RUN_NAME
bash -n scratch/submit_S4xS1xS1_ACS_v2.sh scratch/submit_S4xS1xS1_ACS_audit.sh
```

The existing `extra/test_all.py` had ten dtype mismatch failures before these
changes (17 tests passed across that file plus the eight product tests).
Those legacy chart-model failures are not evidence that this new runner works;
the focused suite and actual bounded runs provide that evidence. Conversely,
they must not be reported as new regressions or hidden under a blanket
“all repository tests pass.”
