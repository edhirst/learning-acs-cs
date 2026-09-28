# S4 × T2 precision campaign

Implementation plan (2026-09-28):

1. Preserve the existing ACS and CS runners. Add a separate precision campaign
   with strict configuration, hash-bound inputs, explicit finite-budget outcomes,
   and independently evaluable checkpoints.
2. Cross-check ambient Nijenhuis components against local coordinate AD and a
   central finite-difference sweep on common random/worst-point probes. The local
   chart is orthonormal at its centre, so tensor norms are directly comparable.
   Include a nonconstant integrable six-torus control and a nonintegrable control.
3. Compare effective-batch Adam, deterministic chunked L-BFGS and matrix-free
   damped Gauss–Newton from the same tail checkpoint, data and selection probes.
   Keep collocation and weights fixed inside line searches/linear solves.
4. Select only numerically admissible, jointly improving candidates. Continue
   the selected method in bounded blocks with explicit top-tail weights and
   separate training/selection adversarial searches. Persist intermediate best
   checkpoints. Each block intentionally restarts its optimizer; this is not an
   exact optimizer resume across blocks.
5. If improvement stalls, allow one zero-initialized anticommuting correction
   layer about the retained learned field. Retain all input derivatives through
   the frozen parent. Do not change homotopy class or introduce an analytic ACS.
   Persistent plateau or excessive norm growth ends automatic training.
6. Freeze the selected artifact and run independent uniform/structured and three
   adversarial audit jobs. A separate report checks hashes, tensor cross-checks,
   ACS/tangency and the agreed RMS/max-found thresholds. No finite test grants an
   integrability certificate or a rigorous supremum bound.
7. Submit dependency-linked jobs with `afterany` and application-level receipts:
   safe review stops propagate as no-ops, failures block descendants, and the
   final summary always reports missing or failed work. Pin the Git revision and
   source/config content; do not silently continue against changed code.

The target is RMS <= 0.01, maximum found <= 0.01, and ACS/tangency <= 1e-10 in
the fixed unit-sphere round product metric. Norm/conditioning guards are
operational safeguards, not geometric nonexistence statements.

## Implemented files and scope

- `extra/cs_precision_fields.py`: existing-checkpoint loader, frozen-parent
  correction head, coordinate and finite-difference tensor verification.
- `extra/cs_precision_optim.py`: fixed-weight, chunked objectives, CG and damped
  matrix-free Gauss–Newton. Rejected GN trials do not mutate parameters.
- `extra/cs_precision_inputs.py`: standard-library preflight for the completed
  s407 campaign, including both diagnostic/checkpoint hashes and the previous
  checkpoint. No Torch or training on the login node.
- `extra/cs_precision.py`: preparation, pilots, selection, continuation, audit
  parts, receipts and final summary.
- `configs/ambient/s4xs1xs1_cs/precision.yaml`: explicit budgets and controls.
- `extra/test_cs_precision.py`: real numerical tests and mocked PBS tests.
- Two ignored scripts in `scratch/`: launcher and compute-job wrapper.

The existing long-run implementation/configurations are unchanged. New Adam has
an independent configurable `lr_min=1e-7`; legacy schedules retain their old
semantics. New checkpoints use a `precision` architecture descriptor and must
be loaded through `cs_precision_fields.load_model`, not the legacy CLI loader.
All automatically scheduled audits support both original and correction models.

The initial source is:

```
outputs/s4xs1xs1_cs/long_replay_s407_F7mKwAFO/run
```

The learned ACS and its full audit remain unchanged. Preparation compares the
previous, average and tail CS fields on common random/structured/search points;
it verifies tensor components at shared random and saved worst-point locations.
A nonlinear pullback CS on T6 and a nonintegrable deformation are independent
verification controls only; neither supplies training labels or an initial field.

## Automatic decisions

The launcher submits 16 jobs with the default six rounds:

```
prepare -> {Adam, L-BFGS, GN} -> select -> round 1 -> ... -> round 6
                                                       |
                         {uniform audit, searches 1, 2, 3} -> summary
```

The three pilots use the same parent, initial collocation and random streams,
with 12.5% training witness weight and no additional top-tail penalty. They have
the same 150-minute working budget, plus time for final selection/reporting.
They also have optimizer-specific iteration caps (5,000 Adam, 500 L-BFGS, 100 GN),
so report actual elapsed time and accepted work rather than equating iterations.
All work has a 210-minute application deadline inside a four-hour PBS allocation.
GN is validated after every step, L-BFGS at least every five steps, and Adam
every 25 steps, so expensive steps do not delay best-checkpoint retention.
Training searches refresh every 100 Adam steps or ten deterministic steps;
tail weights are reconsidered every five deterministic steps. Unchanged weights
do not reset L-BFGS history. Preparation exercises the actual GN derivative
kernels on the cluster before releasing the pilots.

Each comparison re-evaluates candidates on the same expanded selection set.
Promotion requires at least 1% improvement in BOTH RMS and maximum found, ACS
and tangency <= 1e-10, and source-relative conditioning gates. The norm gate is
1.5 times the original tail field's sampled maximum J norm; the G-condition gate
is 2 times its maximum, evaluated on the same points. Rebasing does not reset
these gates. They are conservative operational choices, not universal bounds on
valid complex structures.

If a pilot qualifies, round 1 checks the selected optimizer from the average
parent. Subsequent rounds continue the retained winner with new training RNG.
The uniform-plus-top-10%-tail objective has tail weights 0.1, 0.2, 0.4, 0.8, 1,
1 across the rounds. Collocation and weights are fixed within an optimizer solve;
L-BFGS history resets when the objective changes. Selection searches never enter
training. Historical diagnostic witnesses are explicitly training data now.

No qualifying pilot, or two unproductive continuation blocks, triggers ONE
zero-output correction head with hidden widths [64,64], optimized by damped GN.
It retains the entire parent field and its spatial derivatives, and projects
new tangent generators onto the anticommuting directions at that parent. A
plateau after this attempt stops automatic training. Later queued rounds become
safe no-ops, and audits still inspect the retained candidate.

The chain does not automatically enlarge networks indefinitely, change the
metric, relax tolerances, change ACS homotopy class, or infer nonexistence. A
new spectral-feature architecture or another learned ACS class requires review.

Every block is a declared optimizer restart from a selected checkpoint, not an
exact continuation of the previous optimizer state. `best.pt` is saved during
validation. A successful block exports `final.pt` and a hash-bound receipt.
An interrupted/failed block DOES NOT allow the next job to consume a partial
checkpoint. Failed jobs leave a `failure_*.json` and logs for review; resubmitting
the whole launcher creates a NEW campaign and is not a recovery operation.

## Independent final evidence

After round 6 (possibly a no-op following a review stop), all audit jobs read the
same frozen candidate. One job evaluates 100,000 fresh uniform points and a
32-point grid on each circle at sphere poles and nearby points. Three separate
jobs each use 16,384 pool points and up to 512 diverse starts, with 200 ascent
steps. Pool points are retained too: a search cannot hide a larger pool residual.
Each audit cross-checks its eight worst points independently and compares field
conditioning with the original tail field on the same points.

The report requires all jobs/artifact hashes to validate, all four audits to
agree on the checkpoint, the full audit budget, passing independent tensor and
growth checks, RMS/max-found <= 0.01, and ACS/tangency <= 1e-10. Reduced test
budgets cannot grant acceptance. Even a passing report sets
`integrability_certified=false`; maximum found is not a global supremum bound.

## Publish from the Mac

```bash
cd /Users/ed/Documents/Projects/DG/Nijenhuis/github
git add extra/cs_precision.py extra/cs_precision_fields.py extra/cs_precision_optim.py extra/cs_precision_inputs.py extra/test_cs_precision.py configs/ambient/s4xs1xs1_cs/precision.yaml docs/s4xt2_cs_precision.md
git commit -m "Add autonomous audited CS precision campaign"
git push origin s4t2
```

After the push succeeds, copy the ignored scripts separately:

```bash
scp scratch/launch_S4xS1xS1_CS_precision.sh scratch/submit_S4xS1xS1_CS_precision.sh ehirst@lovelace:/home/lovelace/proj/proj925/ehirst/nijenhuis/
```

## Launch on Lovelace

Copy only these commands, without terminal prompt characters:

```bash
(
set -euo pipefail
cd /home/lovelace/proj/proj925/ehirst/nijenhuis/github
git switch s4t2
git pull --ff-only origin s4t2
git cat-file -e HEAD:extra/cs_precision.py
cd ..
bash launch_S4xS1xS1_CS_precision.sh
qstat -u "$USER"
)
```

Use **bash**, not qsub, for the launcher. It performs preflight and submits the
whole dependency graph. Each compute job requests serial/10 GB/four hours;
the summary requests 30 minutes. The maximum reserved total is 60.5 CPU-hours,
not a predicted runtime. Pilots/audits may run concurrently if the queue permits.
No additional held-job commands are needed. `afterany` ensures reporting can
run after failures; Python receipts enforce scientific/data dependencies.

The launcher prints an exact command of the form:

```bash
cat /home/lovelace/proj/proj925/ehirst/nijenhuis/outputs/s4xs1xs1_cs/precision_s507_XXXXXXXX/comparison.txt
```

Use the actual printed suffix. Send `comparison.txt` back when finished; it
includes the selected artifact path in the accompanying `comparison.json`.
`submission.tsv` records job IDs and `logs/` has per-job output. Do not pull or
edit Python code in the checkout while this campaign runs: revision/content
checks deliberately reject mixed-code work. The configuration is copied into
the new campaign directory at submission and must also remain unchanged.

If qsub fails part-way through submission, inspect `qstat` and `submission.tsv`
before taking further action. Already submitted jobs are NOT automatically
cancelled. A failed/missing receipt causes downstream computation to fail closed;
the summary reports the blockage instead of presenting partial work as success.

## Verification

The focused new/existing ACS/CS regression suite passed 121 tests locally.
Coverage includes coordinate/finite-difference agreement, detection of a wrong
tensor, recovery of a six-dimensional integrable control below 1e-6 on fresh
points, GN finite-difference and adjoint identities, actual updates with all
three optimizers, correction initialization/derivatives/serialization, full
small-budget campaign execution, immutable inputs, failed provenance, terminal
no-ops, all four audit jobs, and mocked 16-job PBS dependencies and failures.
Shell syntax checks and Python compilation also passed. The matrix-exponential
batching warnings on local PyTorch concern performance, not failed assertions.

These are synthetic/local verification fixtures, not training on the remote
learned field. The preparation job validates the actual HPC artifacts and Torch
kernels before releasing the pilots. No remote jobs were submitted by the code
author during implementation, and no success at the research target is implied.
