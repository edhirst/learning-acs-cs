# Resumable S4 × S1 × S1 CS refinement

## Evidence and plan

The s507 precision campaign audited RMS 0.288767385 and maximum found
0.784116698. Its 150-minute blocks completed only 900 Adam, 48 L-BFGS and
10 Gauss–Newton pilot updates; the correction received five GN updates.
The rule requiring simultaneous RMS/tail improvement discarded a monitored
tradeoff of RMS 0.289626 / tail 0.600297 versus 0.286719 / 0.747038.
Those intermediate weights were not retained or independently audited.

This is a new campaign, not an in-place edit of those results. The implementation:

1. Validates the old campaign, four full audits, chosen checkpoint, ACS lineage
   and historical point-cloud hashes using a standard-library login preflight.
2. Rechecks the selected field with coordinate/finite-difference derivatives;
   builds common selection sets with fresh random samples/searches. Old audit
   points are now explicitly training data. Independent evidence for new fields
   comes from fresh uniform samples and fresh-start searches; deterministic
   structured grids are supplemental stress tests and can overlap earlier grids.
3. Runs three Adam branches: low tail weight 0.05, tail weight 0.1, and an
   identity-start [64,64] correction with weight 0.05 and frozen parent parameters
   (but NOT frozen spatial derivatives). All use the existing ambient geometry,
   round product metric, and unnormalised signed Nijenhuis tensor.
4. Chains four allocations per branch. Model, Adam moments, cumulative update
   count, learning-rate schedule, RNG states and replay points persist. Branches
   continue from latest validated iterates, never silently from selected extrema.
   The initial LR is 3e-6, decaying by cumulative update count to 1e-7 at 30,000
   updates. Each allocation allows up to 10,000 updates or 170 minutes of work,
   whichever comes first. A time-limited block is not convergence evidence.
5. Retains separate best RMS, tail and balanced nominations. They are **not** a
   complete Pareto frontier. They allow useful tradeoffs, not just joint gains.
   Small common monitoring (512 uniform + 128 challenge) occurs every 250
   updates; replay searches every 1,000. Full validation occurs at block ends.
   Full-validation extrema are also retained across allocation boundaries, so
   a later cheap-monitor winner cannot discard an earlier full-validation winner.
   Phase timings and runtime versions are logged. Direct chunked backward avoids functional parameter
   vector reconstruction in Adam. No GPU speedup is assumed.
6. Evaluates source, extrema and latest candidates on the same expanded selection
   set. Subject to unchanged ACS/tangency and source-relative conditioning gates,
   selects minimum max(RMS, maximum found), matching the two equal absolute
   targets. This need not improve both quantities relative to the old checkpoint.
7. Runs independent audits: 100,000 uniform, 20,480 structured, and three fresh
   searches with 16,384 pools, 512 starts and 200 steps, plus derivative checks.
   Summary passes only if RMS AND maximum found <= 0.01, ACS/tangency <= 1e-10,
   provenance, derivative and growth checks pass, and the full audit completes.

The squared mean loss corresponding to RMS 0.01 is 1e-4. No rescaling of the
reported tensor or relaxation of targets is used. Passing is finite sampled
numerical evidence, not a rigorous supremum bound, a proof of integrability, or
a guarantee of proximity to an exact complex structure. The existence of a
small-residual field in this search family is not guaranteed.

## Resources and failure handling

Default graph: preparation → three parallel chains of four training jobs →
selection → four parallel audits → summary (19 jobs). Each job requests serial,
10G, 4 hours, except the 30-minute summary: at most 72.5 requested single-thread
job-hours. Actual runtime can be shorter. Dependencies use `afterany`; Python
receipts gate downstream work, and the summary reports missing/invalid stages.
No automatic plateau stop based on a single discarded tradeoff is used.
Numerical/growth gate violations roll back to a validated state and stop that
branch for review. External interruption produces no success receipt; `last.pt`
is a recovery artifact, not authorization to silently skip a failed predecessor.

Do not change Python code/configs or pull the remote checkout while jobs run.
Do not delete old runs: their hashes and base models are required for provenance.
Do not pass these checkpoints to legacy diagnostic loaders; use this campaign's
precision-compatible loader. Old scripts and runners are unchanged.

## Submission

Commit/push the new Python files, YAML, tests, and this runbook on `s4t2`.
Both new scripts live in ignored `scratch/` and must be copied separately:

```bash
scp scratch/launch_S4xS1xS1_CS_refine.sh scratch/submit_S4xS1xS1_CS_refine.sh ehirst@lovelace:/home/lovelace/proj/proj925/ehirst/nijenhuis/
```

Then on the HPC:

```bash
(
set -euo pipefail
cd /home/lovelace/proj/proj925/ehirst/nijenhuis/github
git switch s4t2
git pull --ff-only origin s4t2
git cat-file -e HEAD:extra/cs_refine.py
cd ..
bash launch_S4xS1xS1_CS_refine.sh
qstat -u "$USER"
)
```

The launcher defaults to the reviewed `precision_s507_Y8VCEIUP` source. It prints
the unique output root and exact `cat .../comparison.txt` command. Submit once;
no additional `qsub`, holds, or manual handoffs are needed. If submission fails,
inspect `submission.tsv` and `qstat` first: already-submitted jobs are not cancelled.

## Verification

`extra/test_cs_refine.py` covers strict inputs, actual model/optimizer continuation,
RNG/schedule persistence, separate extrema, correction parent freezing,
fail-closed receipts, selection/audits, small-budget non-certification, and PBS
dependency mocks. Run alongside `extra/test_cs_precision.py` and the prior ACS/CS
focused suite. Synthetic tests do not predict convergence on the HPC checkpoint.

Local verification: 121 existing focused ACS/CS tests passed. The new runner's
20 tests passed after the final recovery/nomination changes, including exact
uninterrupted-versus-split weights, Adam moments and RNG, and a mid-batch
wall-budget interruption before the first scheduled observation. Python compile
and Bash syntax checks also passed. Local Torch emits its existing matrix-exp
batching performance warnings; the HPC runtime/version and phase timings are
recorded separately. No new HPC training has been run locally.
