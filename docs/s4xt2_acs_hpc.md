# Ambient ACS campaign: implementation and Lovelace handoff

**Current handoff:** the seed-7 lineage now passed its full numerical ACS audit.
Use [the production ACS recipe and audited CS smoke instructions](s4xt2_reproduce_and_cs.md)
for new runs. The historical pilot/follow-up instructions below are retained
for provenance, not the recommended default for new ACS runs.

The next experiment is **ACS only**. The four pilots are independent; an
optional summary job waits for all four using an `afterany` dependency.
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
  extra/validate_ambient_acs.py extra/test_acs_training.py extra/test_acs_pbs.py \
  configs/ambient/s4xs1xs1_acs docs/s4xt2_ambient_acs_plan.md \
  docs/s4xt2_research_plan.md docs/s4xt2_acs_hpc.md
git commit -m "Add reproducible ambient ACS pilot campaign"
git cat-file -e HEAD:extra/acs_training.py
git cat-file -e HEAD:configs/ambient/s4xs1xs1_acs/pilot_trilin_tail.yaml
git push origin s4t2
```

The branch is now `s4t2`, as selected for this campaign. Check
`git branch --show-current` before committing/pushing; do not accidentally push
the old `no_overlaps` branch. `git status --short` entries starting with `??`
are untracked and will not be included by `git commit -a` or `git push`.
The `git cat-file` checks above verify that the files actually exist in the
commit, not just in the working directory.

Copy the four ignored scripts separately, using your working SSH hostname/alias:

```bash
scp scratch/submit_S4xS1xS1_ACS_v2.sh scratch/submit_S4xS1xS1_ACS_audit.sh \
  scratch/submit_S4xS1xS1_ACS_pilots.sh scratch/submit_S4xS1xS1_ACS_summary.sh \
  ehirst@lovelace:/home/lovelace/proj/proj925/ehirst/nijenhuis/
```

On Lovelace:

```bash
cd /home/lovelace/proj/proj925/ehirst/nijenhuis/github
git branch --show-current
git status --short
test "$(git branch --show-current)" = s4t2 && git pull --ff-only origin s4t2
git cat-file -e HEAD:extra/acs_training.py
git cat-file -e HEAD:configs/ambient/s4xs1xs1_acs/pilot_trilin_tail.yaml
test -f extra/acs_training.py
test -f configs/ambient/s4xs1xs1_acs/pilot_trilin_tail.yaml
cd ..
```

The branch should be `s4t2`; stop and resolve a different branch or local
conflicts rather than resetting your checkout. Pulling a branch into a different
checked-out branch is not the same as switching to that branch.

### If you see “Missing extra/acs_training.py”

Being on `s4t2` does not establish that its latest commit contains the new files.
At the time this error was reported, this local checkout had the new runner and
configs as **untracked files**, with HEAD still at `1383f80`. A push of that
branch alone cannot publish them. Add/commit/push the explicit file list above,
then pull and verify `git cat-file` on Lovelace before resubmitting.

During the fix, the implementation was committed as `f1cf650` (`astra plan`);
that commit does contain the runner and configs. The follow-up changes to this
runbook, `extra/campaign_summary.py`, and `extra/test_acs_pbs.py` also need their
own commit/push. The four `scratch/` scripts still need copying separately.

If the file is committed and present remotely, the remaining possibility is
checkout discovery: the job searches the directory where `qsub` was called and
its `github/` child. Override that explicitly with
`-v REPO_DIR=/home/lovelace/proj/proj925/ehirst/nijenhuis/github` (combine it with
other `-v` variables using commas). The updated scripts print the searched paths
on failure. Reading a `#PBS -k oe` log from your home directory does not by itself
tell you where the job was submitted.

## 1. Submit a short smoke test

Both scripts request the previously working `serial`, `mem=10G`, one-thread setup,
with no `nodes`/`ppn` directive. They additionally request four hours. This limit
has **not** been verified against your live queue; inspect `qstat -Qf serial` or
your site's queue documentation if rejected. If changing the PBS wall limit,
also set `WALL_MINUTES` to no more than 85% of its minutes.

```bash
qsub -v REPO_DIR=/home/lovelace/proj/proj925/ehirst/nijenhuis/github,RECIPE=smoke,SEED=7,WALL_MINUTES=10 submit_S4xS1xS1_ACS_v2.sh
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

After the smoke succeeds, run this launcher on the **login node with bash**, not
with qsub. It submits all four pilots and their held summary, saves the returned
job IDs, and refuses to construct a dependency if any submission fails:

```bash
cd /home/lovelace/proj/proj925/ehirst/nijenhuis
bash submit_S4xS1xS1_ACS_pilots.sh
```

The default seed is 7; use `SEED=8 bash ...` for a later approved repetition.
Each invocation uses a new timestamped campaign root, preventing old/missing
results from being confused with this campaign. The launcher prints the exact
summary report path and writes `submission.tsv` alongside it. If submission
partly fails, already-submitted jobs remain queued; inspect that manifest and
`qstat` rather than blindly rerunning all four.

Alternatively, these manual commands submit only the pilots (no automatic
summary). They do not depend on each other and none should be held for a CS job:

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

## 3. Automatic held summary; pause before the long campaign

The launcher already submits the following equivalent scheduler dependency:

```bash
qsub -W "depend=afterany:JOB_ID_0:JOB_ID_1:JOB_ID_2:JOB_ID_3" \
  -v "REPO_DIR=/absolute/path/to/github,SUMMARY_ROOT=/absolute/campaign/root,PILOT_SEED=7" \
  submit_S4xS1xS1_ACS_summary.sh
```

Use real full job IDs, not these placeholders. `afterany` waits for all four jobs
to terminate even if they exit with errors or code 75; `afterok` would suppress
the summary in those cases. The dependency hold releases automatically: **do not
add `-h`**, which would introduce a separate user hold. The summary reports
missing results and includes incomplete statuses and actual phase counters.

If the pilots were submitted manually in the older layout, omit `PILOT_SEED`
and set `SUMMARY_ROOT` to the existing
`/home/lovelace/proj/proj925/ehirst/nijenhuis/outputs/s4xs1xs1_acs/pilots`.
`PILOT_SEED` enables exact expected-path checks for the launcher's new layout only.

PBS caveat: a job deleted before it ever runs may not satisfy `afterany`.
If you cancel such a predecessor, inspect/recreate the summary dependency rather
than assuming it will release. Submit the summary immediately after the pilots,
while their IDs are known to PBS. These semantics are documented in the
[PBS Professional reference guide](https://help.altair.com/2022.1.0/PBS%20Professional/PBSReferenceGuide2022.1.pdf).

The summary job requests one GB and ten minutes and writes
`SUMMARY_ROOT/summary.<summary-job-id>.txt`. This is just a report: it does not
submit CS or declare a mathematical certificate.

To summarize manually, with the `nijenhuis` conda environment active:

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

## Follow-up after the seed-7 quick audit

The reported campaign `pilots_20260925_181049_3240935` completed all four pilots.
The uniform `all2+trilin` recipe was the strongest on this seed: training
validation RMS 2.82e-4 and maximum found 6.78e-4. The five L-BFGS iterations were
rejected because both comparison RMS and maximum worsened; the previous weights
were restored. Do not force that finisher to continue or resume a completed run.

Independent quick audit `final.audit.1023437.ada.json` (seed 800007) reported:

| Diagnostic | Observed |
| --- | ---: |
| Fresh uniform raw RMS | 2.779e-4 |
| Raw maximum found | 6.783e-4 |
| Corrected identity maximum | 2.351e-15 |
| Corrected tangency maximum | 7.240e-16 |
| Independent matrix-root difference | 4.815e-15 |
| First derivative finite-difference error | 6.521e-11 |
| Second derivative finite-difference error | 5.304e-6 |
| Tangent condition number at tested worst points | 1.0074 |

The correction and derivative checks pass their thresholds. The chart-overlap
score near 1e-10 is the legacy overlap helper's reporting floor, below the audit
threshold of 1e-9. `screening_targets_passed` remains false because raw uniform
RMS exceeds 1e-4. `numerically_validated` also cannot be granted by a quick
audit. This is encouraging numerical evidence, not a uniform global certificate
or an integrability/CS result. Keep the original candidate and audit unchanged.

### Prepare the three-job follow-up

Publish the new tracked configuration/runner/summary/test/doc changes on `s4t2`,
then copy the **updated** ignored wrappers:

```bash
git add extra/acs_training.py extra/campaign_summary.py \
  extra/test_acs_training.py extra/test_acs_pbs.py \
  configs/ambient/s4xs1xs1_acs/polish_trilin_uniform.yaml \
  docs/s4xt2_acs_hpc.md docs/s4xt2_ambient_acs_plan.md
git commit -m "Add checked ACS polishing and reproducibility follow-up"
git push origin s4t2
scp scratch/submit_S4xS1xS1_ACS_v2.sh scratch/submit_S4xS1xS1_ACS_pilots.sh \
  scratch/submit_S4xS1xS1_ACS_summary.sh \
  ehirst@lovelace:/home/lovelace/proj/proj925/ehirst/nijenhuis/
```

On Lovelace, pull the commit and verify the new config before submitting:

```bash
cd /home/lovelace/proj/proj925/ehirst/nijenhuis/github
git pull --ff-only origin s4t2
git cat-file -e HEAD:configs/ambient/s4xs1xs1_acs/polish_trilin_uniform.yaml
cd ..
STAGE=followup \
SOURCE_CHECKPOINT="$PWD/outputs/s4xs1xs1_acs/pilots_20260925_181049_3240935/pilot_trilin_uniform/s7/final.pt" \
SOURCE_SHA256=a0cc96594b91b5dccb8be014bd80a971c3991dc56f0023c5e3c79b9fae185a86 \
bash submit_S4xS1xS1_ACS_pilots.sh
```

Run that launcher **once** on the login node, not via qsub. It verifies the source
checksum before any submission, and the Python loader verifies it again inside
the polishing job. It submits:

1. `polish_trilin_uniform/s107`: the existing seed-7 candidate, unchanged at
   initialization, then 10,000 additional refinement steps with uniform samples,
   no tail term and no L-BFGS. Adam starts at 1e-5 (the old schedule's terminal
   rate) and decays to 1e-6. Batch size rises from 192 to 768 to reduce stochastic
   variation. No polynomial/ladder stage is rerun. The optimizer is new; this is
   a continuation, not an exact resume. Seed 107 changes the sampling stream;
   **it is not an independent random-start ACS discovery**. Source checksum and
   source counters are retained in the manifest, results and saved checkpoints.
2. `pilot_trilin_uniform/s8`: the unchanged winning pilot recipe from random weights.
3. `pilot_trilin_uniform/s9`: the same independent experiment with seed 9.
4. A summary held with `afterany` on these three jobs, reading the submission
   manifest so missing continuation/repeat outputs are explicitly reported.

The output root is a new `outputs/s4xs1xs1_acs/followup_<timestamp>_<pid>/`.
Resource requests stay at four hours/10 GB with the three-hour application guard.
Polishing does more work per step than the old pilots; no runtime guarantee is
inferred from their four-minute timings. As before, an interrupted job exits 75
and can be resumed with the same config, seed and OUT, `RESUME=1`, and no INIT_FROM.

Send the held summary and the polishing `result.json`/`progress.jsonl`. If a new
candidate reaches raw RMS <=1e-4 and maximum found <=1e-2, select it for the full
held-out audit; do not consume that audit repeatedly while choosing settings.
If it stalls, inspect the trajectory before extending its budget or changing
the optimizer. No CS jobs or automatic acceptance/promotion are submitted.

Local follow-up verification: 48 focused tests pass, including a two-step run at
the full 4×256 width and batch 768, source-hash rejection, immutable source
weights, continuation/resume provenance, and the three-job held summary. The
combined suite has 57 passes and the same ten previously recorded legacy dtype
failures. Long polishing and repeat runs have not been executed locally.

## Local verification

Verified on 2026-09-25 with Python 3.10 / PyTorch 2.5.1, CPU float64:

- Initial focused tests: **33 passed** (including the eight existing product tests).
- Initial three test files: **42 passed, ten pre-existing dtype failures**.
- Submission-fix regression suite: **41 passed**, adding eight PBS/summary
  tests. These include failed or empty job IDs, explicit checkout paths,
  automatic afterany dependencies, and missing/incomplete pilot reports.
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
python -m pytest -q extra/test_acs_pbs.py
python -m extra.acs_local_checks --out scratch/acs_local_checks/NEW_RUN_NAME
bash -n scratch/submit_S4xS1xS1_ACS_v2.sh scratch/submit_S4xS1xS1_ACS_audit.sh \
  scratch/submit_S4xS1xS1_ACS_pilots.sh scratch/submit_S4xS1xS1_ACS_summary.sh
```

The existing `extra/test_all.py` had ten dtype mismatch failures before these
changes (17 tests passed across that file plus the eight product tests).
Those legacy chart-model failures are not evidence that this new runner works;
the focused suite and actual bounded runs provide that evidence. Conversely,
they must not be reported as new regressions or hidden under a blanket
“all repository tests pass.”
