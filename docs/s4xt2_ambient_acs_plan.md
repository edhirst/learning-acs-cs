# Learn an ambient ACS on S4 × S1 × S1 before returning to CS

Date: 2026-09-25. Planning baseline: `1383f80` on `no_overlaps`.

This document is the active execution plan. It supersedes the earlier proposal
to start CS training from the analytic frame reference. The broader
[research assessment](s4xt2_research_plan.md) retains the literature review and
the deferred CS work. The initial implementation and local checks are now in
place; see the [execution record and HPC handoff](s4xt2_acs_hpc.md). The longer
experiments below remain a gated plan, not completed work. No HPC jobs or CS
training have been submitted by the implementation agent.

### Execution checkpoint (2026-09-25)

Implemented: generic three-factor features; true polynomial-only phase;
field-preserving MLP activation; residual-adaptive replay and tail loss;
phase/RNG-aware restart and continuation; safeguarded L-BFGS; strict YAML pilot
configs; independent numerical audit; and ignored PBS submission scripts.

Verified locally: the coefficient-level exact-ACS representation, phase and
restart invariants, legacy checkpoint compatibility, audit/derivative/root
checks, mocked PBS paths and error propagation, and bounded actual learning.
The focused suite has 33 passing tests. Across all three test files there are
42 passes and the same ten legacy chart-model dtype failures measured before
the patch. There are no new failures in those files.

Actual short runs: 2,000 polynomial updates each on S2 (`lin`) and S1×S1
(`bilin`) produced raw maxima about 4.08e-4 and 6.74e-3. The preliminary torus
`all2` control stalled near 1.415; its failed result was retained, and the
bilinear control tests the simpler existing feature family rather than
pretending arbitrary Adam runs must converge. Three target smoke seeds ran
100+100 updates with maxima 1.718, 1.676 and 2.015: **none is an accurate learned
ACS**. Four 4-update wiring runs exercised the actual 4×256 pilot MLP widths.
Short S3×S1, S6, S4×S2 and S4 regressions also ran without claiming convergence.

Next: user submits the Lovelace smoke, then the four 5k+5k seed-7 pilots, and
returns their results. Stop there before choosing long-run recipes. Source
checkpoint-specific seed-10 continuations, normalized features, GN/LM, full
independent candidate validation and CS remain deferred. Only the auditor's
quick path has been run locally; its full 100k/three-restart budget is reserved
for a genuinely promising learned candidate on the HPC.

### HPC evidence and next checkpoint

The user has now completed the seed-7 pilots. `all2+trilin` with uniform sampling
was best: raw RMS 2.82e-4, maximum found 6.78e-4. The five-step L-BFGS attempt was
rejected and rolled back, not a job failure. The quick independent audit agreed
on raw residuals and gave corrected identity error 2.35e-15, tangency error
7.24e-16, root agreement 4.81e-15 and passing first/second derivative checks.
The raw RMS gate is still unmet; the quick audit is not final validation.

The next bounded batch is a checksum-verified 10,000-step uniform continuation
at Adam 1e-5 -> 1e-6, batch 768, no L-BFGS, plus unchanged winning-recipe random
starts on seeds 8 and 9. The continuation uses sampling seed 107 and is explicitly
identified as descended from seed 7, not an independent discovery. The original
checkpoint is read-only and a new output directory is required. See the
[follow-up commands](s4xt2_acs_hpc.md#follow-up-after-the-seed-7-quick-audit).
Stop for these HPC results before spending the full audit budget or returning
to CS. No acceptance thresholds have been weakened.

## 1. Objective and fixed experimental rules

Learn a globally defined ACS on S4 × S1 × S1 using the same ambient formulation
as the other manifold ACS runs:

    p = (x,u,v) in R5 × R2 × R2,   |x| = |u| = |v| = 1,
    M_theta(p) in R^(9 × 9),
    K_theta(p) = P(p) M_theta(p) P(p),
    R_theta(p) = K_theta(p)^2 + P(p),
    r_theta(p) = ||R_theta(p)||_F.

Retain the full 81-component ambient output, round product metric, float64,
and factorwise uniform sampling/retraction. Train from random coefficients
using the ACS residual, with sampling and optimization improvements described
below. Do not use analytic field values as labels, initialize training with
the known reference, switch to a learned 6 × 6 frame representation, or
hard-code an exact ACS into the trainable model. Frames may be used offline
to measure tangent spectra and to construct test fixtures.

The ACS exists, and an explicit zero-residual field can be represented by the
proposed polynomial features. This makes the problem attainable within the
model, not just abstractly possible. It does not guarantee that a particular
nonconvex optimizer or finite step budget finds that field. Failed experiments
must lead to diagnosis and changes supported by evidence.

Polar correction remains part of the existing workflow, but must not hide a
poorly learned field. Report raw K and its corrected J separately, and meet
both accuracy criteria below. Defer all CS training and CS-specific model
changes until this ACS campaign succeeds.

## 2. Success criteria and promotion rules

Use the true Frobenius residual for validation. The current helper adds 1e-20
inside its square root, giving a reporting floor near 1e-10. That floor must
not make a roundoff-level acceptance test impossible. Optimize the squared
residual directly; any stabilized norm used for ascent must be labelled.

| Stage | Required evidence | Action |
| --- | --- | --- |
| Representation verified | A coefficient-level test reproduces the analytic ACS using the actual new feature map, with maximum identity error <= 1e-12 in float64 | Proceed to learning tests |
| Optimization progressing | Independent validation tails/maxima improve, not only the mean; gradients and parameters remain finite | Continue a bounded pilot |
| Admissible candidate | Raw maximum found <= 0.5 on fresh points and independent ascent searches | Save an immutable candidate and test polar correction |
| Accurate learned ACS | Raw RMS <= 1e-4 and maximum found <= 1e-2; corrected identity error <= 1e-10 and tangency error <= 1e-12 | Submit separate final validation |
| Primary objective achieved | One random-start model passes the final independent audit, reload, and derivative checks | Preserve the artifact; do not automatically start CS |
| Reproducible method | At least three independent random starts reach the accurate-ACS gate under a recorded recipe | Finish the campaign report and later revisit CS |

The 0.5 threshold provides margin below the sufficient condition r < 1 for
polar correction. It is not a required mathematical bound for every possible
correction method. These are proposed numerical targets under the fixed metric.
Do not lower them after seeing results without recording and explaining that
decision. Budget-limited partial progress is not completion.

Finite sampling and ascent do not bound the true supremum from above. Use
`numerically_validated`, with the validation budget attached, rather than
claiming a rigorous certificate. Existence on this manifold is known
independently; whether the particular learned approximation works everywhere
still needs numerical scrutiny.

## 3. Model and training changes to implement first

### 3.1 Supply the missing three-factor interactions

Add a generic triple-product feature family over every distinct factor triple.
For the present product it has 5 × 2 × 2 = 20 features x_i u_j v_k. Test
`all2+trilin`, combining these with the existing constant, linear, and quadratic
features. It has 75 raw features and 75 × 81 = 6,075 polynomial coefficients.
The feature ordering, factor blocks, normalization, and schema version must
be saved with every checkpoint.

Retain the old feature modes and their loading behavior. Products with fewer
than three factors should either skip this ladder stage in auto mode or reject
an explicit unsupported selection with a clear message. Implement the feature
generator generically, not by hard-coded S4 × S1 × S1 coordinate indices.

The prior assessment derives a raw polynomial M whose projection is an exact
ACS. Use its coefficients in a test of the real model class. Also verify the
parity obstruction: the reference circle-circle block is odd separately in
all three factors, while degree-two raw features cannot have that parity.
This distinguishes representation failure from optimization failure.

Use deterministic feature scaling based on sphere moments if needed. Record
it and apply it identically during reload. The raw quadratic features contain
linear dependencies on the spheres, since the squared coordinates in each
factor sum to 1. If these harm conditioning, remove a documented redundant
square from each factor or use a rank-revealing orthogonalization. Never
attempt to invert a singular feature covariance. Keep the unnormalized mode
as a compatibility baseline.

### 3.2 Make the polynomial phase genuinely polynomial-only

The current `CurriculumNet.forward` always adds its randomly initialized MLP,
even when only the polynomial optimizer is active. Add an explicit stage
mechanism that bypasses the MLP completely during polynomial training and
evaluation. Inactive parameters must not accumulate gradients.

For this target, go directly to selected feature classes in the controlled
experiments instead of spending another full campaign on every already-failed
ladder rung. Keep the generic ladder available for cross-manifold regression.

When refinement is useful, enable an MLP correction with zero initial output,
random hidden layers, and a zero final layer. This preserves the accepted
polynomial field exactly at the switch and still lets the correction learn.
Do not zero the entire ACS network: K = 0 is a stationary point of the squared
ACS residual. Test both field preservation and subsequent nonzero updates.

Save the best polynomial candidate independently. If it already meets the
final ACS target, an MLP is unnecessary. If enabling the MLP worsens the
validation maximum, retain or restore the polynomial candidate.

### 3.3 Train against the high-residual regions

Implement ACS-aware versions of residual-adaptive and persistent adversarial
sampling, using r or r² rather than Nijenhuis. A starting mixed batch of 192
points is 96 uniform, 48 selected from a fresh candidate pool, and 48 replay
points. The pool can initially contain 768 points. These are tunable settings,
not assumptions about optimality.

Return enough distinct ascent endpoints to fill the replay quota; the existing
one-maximum-per-eight-starts implementation does not do this. Deduplicate by
distance on the product, replenish if necessary, and explore perturbed
neighborhoods so the replay set does not collapse to one hotspot. Refresh part
of the set every 250 optimizer steps initially; benchmark the overhead. Stop
gradients through the sample-selection process when updating model weights.

Compare two objectives explicitly:

    L_mean = average(r² on the chosen collocation distribution),
    L_tail = average(r² on an independent uniform component)
             + lambda * average(largest q-fraction of r² on the hard component).

Start q at 0.1 and ramp lambda from 0 to 1 over the first 20% of a phase.
Specify rounding when the hard set is small. The uniform term always uses
uniform points; the second term deliberately emphasizes difficult regions.
Report uniform mean/RMS independently. Neither the mixed objective nor the
top-fraction term is an unbiased estimate of uniform mean loss or a bound
on the true maximum.

Use separate training, selection-validation, adversarial-audit, and final-test
random streams. A change in logging frequency must not change the training
samples. Final audit seeds stay unused until the final job.

### 3.4 Select and preserve candidates by the intended objective

Record RMS, quantiles, maximum found, fraction r >= 1, parameter/field norms,
and actual compute time. Rank eligible candidates first by whether the
admissibility gate passes, then by the audited maximum, with RMS as a
secondary criterion. Retain a few candidates when their tradeoffs differ.

Use a fixed held-out set and an accumulated validation challenge set to compare
checkpoints. Refresh the audit periodically and re-evaluate shortlisted models
on the same points. This avoids comparing noisy maxima from unrelated samples.
Keep the untouched final audit separate from this model-selection process.

At worst points, record singular values and eigenvalues of the tangent
restriction E K E^T. Ambient K has three normal zero modes by construction;
these must not be counted as tangent degeneracy. Residuals above 1 do not by
themselves prove singularity; inspect the tangent restriction.

### 3.5 Optimize and checkpoint robustly

Keep Adam and a scheduled learning rate as the baseline. Begin with the current
polynomial rate 5e-3 and refinement rate 1e-3, checking stability; compare a
factor-of-three reduction only when diagnostics justify it. Start fresh
continuations near 3e-4. These initial settings must be measured in pilots.

Use L-BFGS only after a promising phase, on a fixed collocation set per solve.
Use chunks of at most 50 internal iterations, initially 200 total, refreshing
samples only between solves. Reject nonfinite parameters, losses, or gradients;
rollback on worse validation tails even when the mean improves. Save weights
and optimizer state consistently after any rollback or stage switch.

Separate exact restart from continuation. Restart restores the same schedule,
phase, feature map, optimizer, sample state, and RNG streams. Continuation loads
the chosen weights into a new run and starts a new explicitly recorded
schedule. Reconstruct the saved feature family and do not rerun the ladder.
If adding cubic features to an old checkpoint, migrate coefficients by named
features, initialize new coefficients to zero, and test field preservation.

Checkpoint during the polynomial phase, ascent refreshes where practical,
refinement, and L-BFGS boundaries. Check wall-time in all phases, including
evaluation. Keep per-job wall time separate from cumulative experiment time.
The current long polynomial ladder has no inner wall guard or mid-attempt
checkpoint; extending it without fixing this is unsuitable for long HPC runs.

Use atomic checkpoint replacement, record successful optimizer steps separately
from attempts, and handle scheduler termination signals by requesting a save
at a safe boundary. Keep best, last, and incomplete artifacts distinct.

### 3.6 Keep routine ACS evaluation affordable

The current dense report also computes Nijenhuis, although ACS training only
requires K² + P. Add an ACS-only monitoring path. Frequent training diagnostics
should not differentiate the full field for Nijenhuis or run the entire chart
protocol. Run derivative and chart-consistency checks on finalists; computing
Nijenhuis then is optional descriptive information, not an ACS acceptance test.

This keeps the objective and field representation consistent with other ACS
runs while allowing substantially more ACS training and validation per hour.

## 4. Tests required before expensive runs

Tests of implementation and stochastic convergence experiments have different
roles. Algebraic and restart tests are deterministic release checks. Fixed-seed
learning experiments are recorded as integration/optimization benchmarks;
failure to converge must not be hidden by weakening assertions or changing
seeds. Test results below are requirements, not results already obtained.

| Test group | Tests and expected behavior |
| --- | --- |
| Existing geometry | Retain all eight product tests; verify each factor norm, projector symmetry/idempotence/rank, tangency, finite off-manifold derivatives, and poles/seams |
| Factor-generic features | Verify 20 triple terms and 75 combined features here; compare to hand-enumerated products; test four factors, factor permutations, and unsupported one/two-factor requests |
| Representation | Set known coefficients only inside the test fixture; reproduce the analytic field and raw identity to <= 1e-12; confirm the degree-two parity limitation |
| Feature conditioning | Check scaling and normalization round trips; detect expected quadratic dependencies rather than silently using a singular inverse |
| Polynomial isolation | Changing disabled MLP weights has no effect on K or loss; no MLP gradient accumulation; stage/reload metadata preserves the active forward path |
| Refinement transition | Enabling the zero-output correction preserves K; one or more valid optimizer steps can activate it; rejected refinement restores the saved candidate |
| Algebraic gradients | Compare parameter and tangent-direction derivatives against centered finite differences at nonstationary, nonsingular points; check full and optional skew outputs |
| Mean/tail objective | Use hand-calculated residual arrays to verify squared norms, weights, top-fraction rounding, zero-weight limit, and finite gradients |
| Sampling and replay | Exact requested counts; points remain on the product; no unwanted parameter gradients from ascent; persistence, deduplication, and replenishment work |
| Defect detection | Deliberately scale the reference to zero in a narrow region; show that targeted ascent can expose the known defect and that the acceptance path rejects it |
| Polar correction | Test well-conditioned admissible nonsymmetric tangent perturbations; compare to an independent Schur-based root; verify actual corrected identity, finiteness, tangency, and derivatives |
| Boundary/invalid cases | Candidates with r near or above the policy gate, tangent real eigenvalues, NaNs, infinities, or singular corrections are not silently accepted |
| Checkpoint loading | Round-trip every old/new model kind and phase; strict case/dimension checks; explicit legacy metadata handling; feature migration preserves K |
| Exact restart | A deterministic short CPU run matches an interrupted/restarted run in weights, learning rates, replay, and sample stream within a stated tolerance |
| Fresh continuation | Starts at the requested new learning rate, keeps the selected field, preserves the original directory, and skips completed initialization |
| Wall and rollback | Interrupt inside polynomial/refinement/L-BFGS; exercise a deliberately nonfinite update; keep true step counts and a loadable last accepted state |
| Evaluation isolation | Changing evaluation cadence leaves training RNG unchanged; chunk size does not change aggregate statistics; final audit cannot alter model parameters |
| CLI/config and scheduler wrapper | Validate choices and unknown keys; CLI overrides tracked defaults; wrong checkout/checkpoint fails early; shell quoting, exit statuses, paths with spaces, and missing PBS metadata behave predictably |

Use float64 and tolerances based on scale: roughly 1e-12 for small algebraic
identities and 1e-5 relative error for appropriately scaled finite-difference
gradients, with a step-size sweep and an absolute tolerance near zero. Do not
demand a relative error bound on an identically zero derivative.

Also run the existing broader test suite on the baseline and modified tree.
Classify any baseline failures explicitly; the acceptance requirement is no
new failures attributable to the ACS patch, plus passing all new relevant tests.

## 5. Short learning experiments and pilot design

### Local smoke and positive/negative controls

Run two short random-start target smokes to exercise both polynomial and MLP
phases, adaptive refresh, checkpoint/reload, and final reporting. Use a small
MLP (two layers of width 32) and 100–300 steps per phase. They test execution,
not whether that small budget solves the target.

Then run ambient ACS learning controls on S2 and T2 using known low-degree
feature modes. Require real loss reduction and an accepted learned result
under a modest fixed-seed budget calibrated against the baseline. Extend
regression checks to S3 × S1, S6, and S4 × S2 using at least one recorded seed
each, without changing their default feature modes. Known ACS are validation
fixtures, never supervised targets.

Use S4 as a topological negative control and the manufactured narrow defect as
a deterministic adversarial-control case. A numerical validator accepting an
S4 candidate is evidence of insufficient validation or a bug; reject the
validator, not the mathematical obstruction. Failure of one finite adversary
to find a defect does not invalidate that obstruction.

### Paired S4 × T2 pilots

Use full ambient outputs and seeds 7, 8, 9 for each recipe. Each pilot has up
to 5,000 polynomial steps plus 5,000 refinement steps, with at most 200 L-BFGS
iterations on promising candidates. Default target width is 256 × 4 to match
the existing runs; smoke sizes are not the production architecture. Record
stage time, sample counts, validation cost, and maximum memory.

| Recipe | Features | Sampling | Loss | What the comparison tests |
| --- | --- | --- | --- | --- |
| A0 | all2 | uniform | mean r² | Matched degree-two control with corrected stage semantics |
| A1 | all2+trilin | uniform | mean r² | Effect of missing representation capacity |
| A2 | all2+trilin | mixed uniform/pool/replay | mean on that mixture | Effect of adaptive sampling |
| A3 | all2+trilin | same as A2 | explicit uniform mean plus ramped hard-tail term | Effect of the loss beyond sampling |

There are 12 primary pilots. A0 differs from the old curriculum because its
polynomial stage is actually isolated; label it accordingly. Compare matched
step counts and actual compute time, since adversarial updates add cost.

Run two additional bounded continuations of the remote seed-10 checkpoint:
one uniform and one adaptive/tail, with the same saved architecture and feature
family, fresh schedules, and 10,000 additional refinement steps. Preserve the
original. These diagnose whether the existing weights can be repaired; they
do not replace the random-start success requirement. Missing remote weights
do not block the new random-start experiments.

For a cheap first screen, evaluate all four recipes on seed 7 before expanding
to seeds 8 and 9. Skip recipes only for an explicit fault or recorded lack of
promise; retain the matched control needed to interpret changes.

Pilot selection: prefer repeated admissibility, then lower audited maxima and
continued improvement per wall time. If none reaches 0.5, do not automatically
run all recipes longer. A longer run is justified by continuing tail reduction;
a flat or worsening maximum requires the diagnosis in Section 8 first.

## 6. Longer HPC campaign

Use the winning recipe(s) in resumable stages, with no CS job in the chain.
These are optimizer-step budgets, not predictions of wall time.

| Stage | Runs and budget | Promotion criterion |
| --- | --- | --- |
| HPC smoke | 1–2 jobs, 100–300 steps per phase and a tiny reload/audit | Environment, memory, checkpoints, and path handling pass |
| Controlled pilots | Up to 12 recipes/seeds, 5k polynomial + 5k refinement each; 2 seed-10 continuations if available | Rank recipes from shared independent audits |
| Longer polynomial learning | Best 2 recipes × 3 fresh seeds 20,21,22; up to 20k polynomial steps each | Lower tail/maxima; preserve every admissible checkpoint |
| Refinement | Best up to 3 eligible candidates; up to 20k MLP-refinement steps each | Raw accuracy target met or continued measured progress |
| Extended refinement | Only still-improving finalists; another 20k steps each, reaching a 60k total cap for a 20k-poly run | Pass accuracy target, or trigger a specific diagnosis |
| Independent final audit | Separate job for each selected candidate; no training | All final numerical criteria pass |

If one pilot recipe clearly dominates, use three rather than six fresh long
runs. If the polynomial model reaches the target, skip MLP refinement. Preserve
already-accurate polynomial results even if later refinement fails. Additional
L-BFGS work may total 500 iterations per finalist in bounded chunks, conditional
on validated improvement and the measured memory budget.

Validation schedule for initial long runs: 4,096 uniform points every 500
steps; an audit every 2,000 steps from 128 starts selected from a fresh larger
pool and 50 ascent steps; a stronger 512-start/150-step audit at stage boundaries.
Increase or reduce frequency based on a recorded timing benchmark, without
altering the accuracy gates. Save every 250 steps or at most five minutes
apart, whichever calls for a save sooner.

Retain the currently accepted Lovelace serial settings (10 GB, one numerical
thread, no unverified node request). Read queue wall-time limits before
choosing a PBS walltime. Measure 500 representative steps of each phase plus
one audit and compute memory/time estimates. Configure the application stop
before the scheduler deadline, reserving at least 15% of each allocation for
checkpointing and bounded reporting. A 20k-step phase may span multiple jobs.
Do not copy the current 4,000-minute application setting into an allocation
without checking its actual limit. Peak measured memory should leave a useful
margin below 10 GB; reduce batch/chunk/history sizes if necessary.

## 7. Submission package to implement and test

Keep short site scripts entirely inside gitignored `scratch/`, as requested.
Keep experiment definitions and the loader in tracked files. Use a small,
strictly validated `--config` layer on the same ambient trainer, preserving
its CLI; the legacy atlas YAML files are not appropriate for this pipeline.

Proposed tracked YAML presets under `configs/ambient/s4xs1xs1_acs/`:

- `smoke.yaml`;
- `pilot_all2_uniform.yaml`, `pilot_trilin_uniform.yaml`,
  `pilot_trilin_adaptive.yaml`, `pilot_trilin_tail.yaml`;
- `continue_seed10_uniform.yaml`, `continue_seed10_tail.yaml`;
- `long.yaml`, resolved from the chosen pilot recipe and recording that choice;
- `validate.yaml`, specifying independent audit budgets and acceptance thresholds.

Proposed local submission scripts:

- `scratch/submit_S4xS1xS1_ACS_smoke.sh`;
- `scratch/submit_S4xS1xS1_ACS_pilot.sh`;
- `scratch/submit_S4xS1xS1_ACS_long.sh`;
- `scratch/submit_S4xS1xS1_ACS_validate.sh`;
- an optional small batch-submit helper that checks each returned job ID.

The scripts take a preset, seed, unique run ID, and optional source run.
They locate either a checkout at the submission directory or its `github/`
child, verify the required revision/capabilities and conda environment, and
write under `outputs/s4xs1xs1_acs/<campaign>/<recipe>/s<seed>/<run_id>/`.
They record the resolved configuration and PBS job ID before training. Fail
on output collisions unless exact restart is explicitly requested. Validate
config paths and quote every path; report a missing source checkpoint before
consuming a training allocation where possible.

The handoff must include exact `scp` and `qsub` commands for these implemented
scripts, plus a serial submission loop if job arrays are unavailable. Do not
provide commands for nonexistent flags as if they work now. Test with mocked
PBS variables and then with the actual smoke job before the pilot batch.
The HPC environment is not accessible from this local review; remote timing
and queue validation remain steps of the planned handoff.

Use afterok only for successful execution dependencies such as smoke -> pilot
or training -> audit. Mathematical pass/fail is recorded separately. A final
audit job should return a documented nonzero status on validation failure;
an incomplete training job may be resumable but must not produce a success
marker. Long-run promotion requires a selected recipe/candidate manifest, not
just a predecessor's zero exit code. Do not chain or submit any CS job.

Save `resolved_config.yaml`, versioned model/state checkpoints, `progress.jsonl`,
`result.json`, selection audit data, adversarial endpoints, code/environment
metadata, and final `validation.json`. Keep raw and corrected metrics in
separate fields. A final artifact includes the frozen learned weights, feature
definition, correction recipe, hashes, and a minimal load/evaluate command.

## 8. Escalation if the first recipes still fail

Change one hypothesis at a time, using bounded paired experiments:

1. **Representation:** if the exact-coefficient test fails, repair that before
   training. If the polynomial solution is representable but optimization
   stalls, longer MLP runs are not the first diagnostic.
2. **Conditioning:** examine feature rank/scaling, raw-field norms, gradient
   norms, and tangent spectra. Compare normalized features and a smaller
   learning rate with the same samples and initialization distribution.
3. **Optimizer:** trial damped Gauss–Newton or Levenberg–Marquardt on the
   6,075-parameter polynomial residual using matrix-free Jacobian products.
   Test on small known residual problems and compare improvement per wall
   time. Do not allocate the full residual Jacobian/Hessian without a memory
   estimate. Reuse existing optimizer code only after relevant verification.
4. **Optional skew ambient restriction:** compare K = P(M-M^T)P/2. This is still
   an ambient field, and the known solution is in that family; reducing
   nonnormality could help. Label it as a restricted ACS ablation. It must not
   silently replace the full-matrix baseline or become a restriction on future
   CS training. Release the restriction for a full-model polish if useful.
5. **Defect concentration:** persist diverse hotspot neighborhoods; compare
   stronger tail weighting and a weak annealed derivative penalty, monitoring
   whether defects disappear or merely narrow. Avoid a default large
   Dirichlet weight, which can favor shrinking the field.
6. **Capacity:** if optimization and tails are demonstrably progressing but
   saturating, compare full cubic features or a wider smooth correction. Keep
   the same ambient input/output geometry and run held-out derivative checks.

Analytic coefficients and supervised fitting may diagnose a model or optimizer
inside tests, but cannot count as campaign success. Do not change to frame
learning or an analytic ACS export to satisfy the random-start requirement.
If bounded experiments plateau, preserve their data, state the diagnosis,
and choose the next measured intervention; ACS existence rules out interpreting
their failure as a topological obstruction on this manifold.

## 9. Final independent ACS audit

For each frozen candidate:

1. Evaluate raw and corrected fields on at least 100,000 fresh uniform points
   in chunks; report true RMS, mean, quantiles, maxima, and fractions exceeding
   0.01, 0.5, and 1 for the raw field.
2. Run three independent adversarial audits, each with 1,024 starts selected
   from its own 32,768-point pool and up to 200 ascent steps. Include old
   hotspots and their neighborhoods as additional tests, not substitutes for
   fresh starts. Audit both raw residual and corrected identity/conditioning.
3. Test Cartesian poles, antipodal combinations, circle-angle grids, and
   structured slices. Check periodic seams through ambient coordinates.
4. Verify polar correction against an independent matrix-root method at
   representative and worst-conditioned points, with raw and corrected
   tangency checks. Reject nonfinite or unstable corrections.
5. Check first derivatives by tangent finite differences and automatic
   differentiation; check bounded, finite second derivatives on representative
   and worst regions. Inspect growth across checkpoints. Do not require
   Nijenhuis to be small: this stage is ACS learning only.
6. Run the eight-chart consistency check, keeping it distinct from the ACS
   identity test. Repeat chunk sizes and reload in a fresh process. Hash the
   checkpoint before and after auditing to verify no parameter updates.
7. Compare the shortlist on the same final points and preserve a report
   identifying the training source, method, seed, metrics, costs, and all
   remaining limitations. If final-audit points enter subsequent training,
   retire them from final-test status and generate a new independent audit.

Completion requires the learned-ACS accuracy gate in Section 2, with stable
correction and derivatives. Aim for three successful independent starts to
show the recipe is repeatable, while distinguishing that evidence from the
first successful individual structure. Preserve Chern/topological analysis
and all CS decisions for the later return to integrability.

## 10. Implementation and delivery order

1. Implement generic triple features and polynomial-stage isolation; add
   representation, gradient, and compatibility tests.
2. Implement ACS sampling/objectives, correct replay counts, and test the
   manufactured defect. Add economical ACS-only monitoring.
3. Repair continuation, phase checkpointing, wall guards, candidate selection,
   rollback, and result semantics; verify interrupt/restart behavior.
4. Implement strict tracked configurations and the independent ACS validator;
   create and syntax/mock-test the ignored submission scripts.
5. Run deterministic tests, local learning controls, and target smokes; publish
   the commands, outcomes, baseline regressions, and unresolved failures.
6. Handoff the smoke/pilot submission package, benchmark on Lovelace, run the
   paired pilots, select recipes, and submit bounded longer runs.
7. Audit successful candidates independently, repeat with fresh seeds, and
   deliver the learned ACS artifacts and evidence report. Then revisit CS.

The literature basis remains the constructive existence result in
[Parton](https://arxiv.org/abs/math/0009156), adaptive sampling experiments in
[Wu et al.](https://arxiv.org/abs/2207.10289), and optimizer comparisons in
[Rathore et al.](https://proceedings.mlr.press/v235/rathore24a.html). The latter
two motivate controlled numerical experiments; their PDE results do not
establish convergence guarantees for this algebraic tensor-field search.
