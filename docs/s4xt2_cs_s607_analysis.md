# S4 × T2: checkpoint analysis and the next CS search

2026-09-30. This is a research plan supported by local numerical experiments,
not a successful CS result and not an implemented HPC campaign.

## Conclusion

Do not repeat the previous small-batch, decaying-learning-rate campaign at
larger wall time. Keep both the low-RMS and low-tail checkpoints. The next
implementation should combine an equivalent, faster residual evaluator;
better control of gradient sampling noise; and regularized, line-searched
field updates. Increasing optimizer strength without improving collocation
coverage is demonstrably unsafe for these checkpoints.

The agreed targets remain uniform RMS <= 0.01 AND maximum found <= 0.01,
with ACS/tangency <= 1e-10 in the original unit-round product metric. An RMS
of 0.01 means an unweighted uniform mean-square Nijenhuis loss of 1e-4.
No division by tensor component count, learned metric, or rescaled sphere
radius is permitted. Passing sampled tests would establish an approximately
integrable ACS at stated tolerances, not an exact complex structure theorem.

## 1. Evidence and reproducibility

The supplied model and report archives were safely extracted into the ignored
directory `scratch/cs_s607_analysis_0wj_or15/`. All seven supplied files have
hashes present in the corresponding campaign receipts. Original artifacts
were not edited.

Important CS checkpoint SHA256 values:

- `low_4/last.pt`:
  `5be5768908469aca3fe1891335c16d79a84c6fb572be97d378980b618089d201`.
- `correction_4/best_tail.pt`:
  `4091d0ff0faeb8999c2efcf35acbcf21914966ba7c22d4a0b889adc5abcc819f`.
- `correction_4/last.pt`:
  `8639f5ed361ee52d9e8cbfce01bc3be1e35b40500cefb7b173fa61bc66b78d90`.

The ACS and its full-audit hashes match the original s307 base. Local
experiments use float64, Torch 2.5.1, one Torch thread. HPC timing must be
rechecked under its Torch 2.9.1 build. Some local experiments ran concurrently;
timings are indicative microbenchmarks, not promised cluster speedups.

Existing full-audit results, on different audit samples:

| Campaign | Uniform RMS | Maximum found |
|---|---:|---:|
| s407 tail | 0.309967 | 0.934936 |
| s507 precision | 0.288767 | 0.784117 |
| s607 selected correction | 0.286677 | 0.628805 |

The last campaign used 43,227 updates. The RMS improvement over s507 was
only about 0.72%, although the maximum found improved by about 19.8%.
The latest audited RMS and maximum still need factors of about 28.7 and
62.9 reduction. Extrapolating the recent learning curve is not a credible
way to budget that reduction.

### Fresh local baseline

On the same 2,048 newly generated uniform points (seed 930607):

| Checkpoint | RMS | Maximum in these 2,048 points |
|---|---:|---:|
| Low-RMS checkpoint | 0.278681 | 0.469069 |
| Audited correction checkpoint | 0.286986 | 0.448921 |

These maxima are NOT substitutes for adversarial searches. These points
became development data during the experiments below, so cannot be reused
as independent final audit data.

On 256 of these points, the low checkpoint's squared tensor norm decomposed
as approximately 40.9% sphere–sphere input pairs, 48.2% sphere–circle pairs,
and 10.9% circle–circle pairs (both antisymmetric input orders included).
This is not a diagnosis of a single bad factor or seam. The remaining bulk
residual needs attention as well as the maximum.

### Gradient noise and rare influential points

Eight disjoint uniform batches of 128 gave mean pairwise gradient cosines
of -0.0168 (low checkpoint) and -0.0209 (correction head). The individual
low-model gradient norms were 3.88–6.10, versus 1.70 for their average.
This measures substantial disagreement in the sampled directions; it does
not prove that the population gradient vanishes or that Adam cannot work.
These diagnostics used uniform MSE, not the previous composite tail loss.

For 512 individual-point gradients of the low checkpoint:

- The six largest points contributed 49.3% of the sum of squared gradient
  norms; the largest 52 contributed 68.6%.
- Spearman correlation of gradient norm with J's Frobenius norm was 0.658.
- On this finite pool, unbiased importance sampling with probability
  proportional to a 25% uniform / 75% J-norm-squared mixture reduced the
  calculated single-sample gradient variance by a factor 1.85. The same
  mixture using residual magnitude gave only 1.14. An oracle using gradient
  norms gave 1.98.

This supports testing sensitivity-aware sampling, not claiming a universal
variance reduction. The reported concentration is not a theorem about the
whole manifold. Repeated clipping in the HPC logs alone was not sufficient
to diagnose this. Global clipping and Adam have interacting scale effects;
turning clipping off is not a justified standalone fix.

### Exact acceleration, without replacing the learned ACS

At fixed collocation points, cache the frozen reference field's value AND
its first coordinate derivatives. Its affine first-order germ has exactly
the same first jet there, hence exactly the same Nijenhuis tensor and
parameter derivatives. Caching values alone would be wrong.

Use a smooth global orthonormal tangent frame E (rows 6 × 9) internally:

    J_frame = E J_ambient E^T,
    J_ambient = E^T J_frame E.

The frame is geometry, not an analytic ACS target. The original network,
learned ACS, R9 input, ambient output, field, and metric are unchanged.
Six-dimensional matrix exponentials and six tangent derivatives replace
their larger ambient counterparts. The noncommuting-frame bracket terms
must be included. The signed antisymmetric tensor can be stored in 90
components with a sqrt(2) factor, preserving the full Frobenius norm.

Median 128-point forward/backward times (three repetitions, chunk 32):

| Field | Existing evaluator | Cached frame evaluator | Approximate ratio |
|---|---:|---:|---:|
| Low model, all parameters | 1.122 s | 0.280 s | 4.0× |
| Correction, parent frozen | 2.211 s | 0.146 s | 15.1× |
| Correction, parent unfrozen | 2.596 s | 0.422 s | 6.2× |

Tensor discrepancies were below 1.6e-14; relative full parameter-gradient
discrepancies were below 1.3e-13. A separate ambient first-jet cache also
matched the original gradients. The frozen-parent evaluator avoids
recomputing the parent twice through the current correction implementation.

Additional checks at points with J norm up to 21.56 passed: frame versus
ambient gradients agreed to relative error below 4.3e-13, parameter-gradient
finite differences agreed below 9e-9 at the best tested step, and live
adversarial-input finite differences converged below 4e-9 absolute error.
Independent coordinate tensor checks passed at these points. Four existing
focused control/gradient tests also passed; this was not a whole-repository
test run.

Critical restrictions:

- A first-jet cache is valid for residual/parameter derivatives at its fixed
  points. It is NOT valid for adversarial input gradients, which need second
  derivatives of the reference. Use the original live field for searches.
- If a parent is unfrozen, its first jet cannot remain cached. Cache only
  the unchanged ACS base, and differentiate the trainable parent normally.
- Never use a nearby cached point as an interpolated replacement for a new
  collocation point. Compute a fresh jet at each new point.
- Bind caches to point hashes, reference hashes, geometry/frame version,
  precision, and schema. Rebuilding a reference invalidates its cache.
- Preserve the original conditioning baseline across rebases. Internal
  tangent coordinates must not silently redefine the legacy ambient
  conjugation condition number.

### Optimizer and representation experiments

Three 60-step fixed-data L-BFGS trials used 1,024 uniform training points
plus 128 historical challenge points, a 0.05 witness-MSE weight, and strong
Wolfe line searches. Separate uniform and historical-challenge development
sets were not used by the optimizer. Training objectives decreased, but
generalization deteriorated. These are rejected experiments, not new
candidate structures.

An exploratory linearized, degree-one global-frame correction used 360
coefficients, 256 training points and 512 held-out points. Its analytic
linearization matched a central finite difference to relative error 3.9e-11.
The column-scaled Gram matrix had null directions and a broad spectrum.
Small damping improved training more than held-out error. None of the 18
damping/step combinations improved the held-out baseline RMS of 0.271237.
Thus a small coarse basis is not, by itself, a solution either.

A degree-two probe increased this to 1,980 coefficients. One trial reduced
training RMS from 0.277599 to 0.175512 while increasing held-out RMS from
0.271237 to 1.04608. None of its 18 damping/step combinations improved
held-out RMS. This is particularly important: increasing representation
capacity before addressing coverage can make this problem worse.

Hyperparameter comparisons on these development sets are selection evidence,
not independent validation. Completed pilot figures are recorded below.

| Trial | Updates | Held-out RMS before → after | Held-out uniform maximum before → after |
|---|---:|---:|---:|
| Low model, 1,024 uniform training points | 60 | 0.27694 → 0.38027 | 0.44906 → 3.70875 |
| Frozen-parent correction, same points | 60 | 0.28505 → 0.31177 | 0.44892 → 1.42384 |
| Unfrozen correction + parent, same points | 60 | 0.28505 → 0.38563 | 0.44892 → 3.66865 |
| Low model, 8,192 uniform training points | 20 | 0.27868 → 0.28038 | 0.46907 → 0.59429 |
| Same larger set, coefficient-H1 step penalty 0.1 | 20 | 0.27868 → 0.28022 | 0.46907 → 0.60373 |

The larger trials additionally used 256 historical training witnesses;
their 2,048-point uniform development set was disjoint from training. Their
historical challenge maximum improved from 0.54859 to about 0.467, but that
does not erase the worse uniform behavior. None is promoted. The 0.1 penalty
on frame-coefficient changes and their directional derivatives was not a
sufficient cure; it was a coordinate-H1 prototype, not the covariant
trust-region implementation proposed below. Different local trial sizes
and durations are explicit and should not be confused with matched HPC
efficiency benchmarks.

### Topological diagnostic

Numerically integrate the first Chern form of the +i eigenbundle projector
Pi = (P - iJ)/2 over two torus slices. The learned ACS base gives pairings
at roundoff on both 16 × 16 and 32 × 32 grids. The low CS checkpoint gives
pairings below 6e-11 on the coarser grid and below 3e-17 on the finer grid.
This is consistent with c1 pairing zero on the torus generator, not a
rigorous integer certificate or a classification of its full homotopy class.

Every finite sequence of smooth exponential conjugations is homotopic
through ACS to its base. More layers or a new CS seed cannot leave that
base's homotopy class. There is no evidence yet that this is the active
bottleneck; it is a necessary fallback investigation if well-resolved,
larger-capacity optimization genuinely stalls.

## 2. Recommended implementation, in priority order

### A. Build and verify the equivalent fast evaluator first

Add isolated, tested modules rather than rewriting successful ACS training:

1. `extra/cs_jet_cache.py`: content-addressed frozen jets and explicit cache
   invalidation; bounded/chunked construction and exact resume metadata.
2. A registered S4 × T2 orthonormal-frame helper, with direct ambient
   equivalence tests. Keep ambient exports and ambient independent audits.
3. `extra/cs_resolved.py`: a new runner using the fast evaluator and strict
   configuration. Existing campaign results remain readable and unchanged.
4. Safe checkpoint loaders, streaming optimizers, full optimizer/RNG state,
   accepted-step counts, closure/CG counts, and wall-time accounting.

New correction architectures must round-trip through the existing independent
ambient audit path. An internal 6 × 6 generator can export C_ambient=E^T C E
as an ordinary tangent conjugation stage, preserving a live global ambient
field. Record the architecture explicitly; do not make an audit depend on a
cache or a collection of sampled J values.

Do not silently reuse a first-order cache in the existing adversarial search.
Prove numerical equivalence before using the speedup to enlarge datasets.

### B. Reduce sampling noise and explicitly cover sensitive regions

Use the saved low-RMS model as the primary starting point and retain the
audited correction as the lower-tail branch.

- Start with a 32,768-point training pool, combining uniform draws and a
  carefully tested uniform-measure randomized low-discrepancy construction.
  Do not sample spherical polar angles uniformly: that is not uniform area.
  One explicit construction uses seven scrambled Sobol coordinates: transform
  five through the normal inverse CDF and normalize to S4; use the other two
  as uniform circle angles. Test moments, endpoint handling, independent
  scrambles and random product rotations. The Sobol index/scramble belongs
  in resumable state. Compare independent scrambles to estimate quadrature
  error before declaring an apparent stationary point.
- Add separate, diversified residual witnesses and sensitivity witnesses.
  J-norm-squared is a cheap initial sensitivity proxy, not the final target.
- For stochastic training use effective uniform batches of 1,024, increasing
  to 4,096 if two independent batch gradients remain poorly aligned.
  Keep chunks 32–64 for memory; accumulation provides the effective batch.
- For a uniform finite pool of M points and sampling probabilities q_i,
  weight each stochastic uniform-loss contribution by 1/(M q_i).
  A 25% uniform mixture bounds this weight by 4. Detach sampling probabilities
  and weights. Tail penalties remain a separately declared objective.
- Compare this with uniform sampling under the same compute budget. Do not
  infer that a 1.85× finite-pool variance result guarantees faster training.
- Refresh training pools and searches between optimization blocks, not
  inside line searches or Krylov solves. Keep validation/search RNG separate.

### C. Compare three controlled update strategies

Do not send all three into indefinite long runs. Use matched short pilots
on both retained checkpoints, with equal wall budgets and common evaluation.

**Pilot 1 — importance-sampled, larger-batch Adam.** Preserve a matched
uniform control. Try declared optimizer restarts at LR 3e-6 and 1e-5 rather
than continuing an already tiny decay to 1e-7. Choose rates by pilot evidence;
do not assume either is safe. Compare clipping policies while logging actual
field displacement, not just parameter-gradient norms.

**Pilot 2 — larger-collocation L-BFGS with a field-space trust penalty.**
Use 8,192–32,768 uniform points plus separate hard/sensitive points. Bound
the update to J and its first derivatives relative to the block's initial
field; refresh that reference after accepted blocks. This penalty controls
an update, not the existence of a nonparallel or nonorthogonal target.
Reset quasi-Newton history when the objective/data change. Reject updates
creating new validation peaks, even if the closure loss falls.

**Pilot 3 — regularized function-space Gauss–Newton.** Solve for a correction
to the current field using signed tensor components, not scalar norms:

    minimize ||N(J) + D N_J(delta J)||_W^2
             + lambda0 ||delta J||_L2^2
             + lambda1 ||nabla delta J||_L2^2,
    subject to J delta J + delta J J = 0.

Realize the accepted correction by an ACS-preserving conjugation, and
measure the actual nonlinear residual. Use trust-region agreement between
predicted and observed decrease to adjust damping. Evaluate on common
held-out data as well as the solve points.

Use diagonal plus low-rank/Nystrom or block preconditioning, a matrix-free
operator, and streaming accumulation. Never build a dense full network
Jacobian under the 10 GB queue limit. Start with the existing learned
features or modest correction heads; select new directions by their ability
to reduce held-out linearized residual, not arbitrary width doubling.

Initial numerical settings to test explicitly: PCG relative tolerance 0.1,
at most 30 iterations, rank-32 preconditioner refreshed when its measured
effectiveness decays; field-step damping screened at 10, 1, 0.1, 0.01 and
0.001. Use the full quadratic prediction
`||r||^2 - ||r + A delta||^2`, not just `-2 g dot delta`, in the gain ratio.
Reject ratios below 0.1, increase damping on rejection, and permit larger
steps only after good agreement and held-out stability. Prototype relative
step caps are 0.05 for J and 0.1 for its covariant derivative, with denominator
`1 + reference norm`, evaluated on collocation and selection probes. These
are operational caps to calibrate, not global bounds or existence constraints.

The coordinate choice can be further improved by normalizing a reference
with its J-invariant positive metric h=(I+J^T J)/2 before choosing a general
anticommuting generator, then converting back. This is a proposed ablation,
not a verified cure. It MUST allow general real generators, not only skew
ones; otherwise it would impose an unwanted orthogonality restriction.

A frozen reference is not intrinsically pointwise too restrictive. Every
tangent ACS variation D satisfies JD+DJ=0 and is generated by C=-DJ/2,
since [C,J]=D. The current head is therefore locally general at a point;
its finite feature representation and conditioning may still be limiting.
The pilots must distinguish that issue from mere freezing of the reference.

### D. Increase representation only when the diagnostics justify it

Use the following decision rules before automatically choosing an escalation:

| Observation | Next action |
|---|---|
| Linearized/training error falls, held-out error rises | Improve quadrature, sensitivity coverage and update regularization; do not merely widen the model |
| Resolved linearized error barely falls on training or held-out data | Test richer correction directions and conditioning; topology is a later investigation |
| Linear prediction improves, actual nonlinear error worsens | Reduce the field step/increase damping; inspect derivative checks |
| RMS improves but new peaks grow | Refresh diversified tail constraints and retain a separate low-tail lineage |
| Both targets pass on selection data | Freeze and run genuinely independent audits |

The current small-data experiments primarily exhibit the FIRST row, not
evidence for a hard architectural or topological floor.

If a small correction cannot fit even well-resolved training residuals,
compare trainable-parent and frozen-parent models under the SAME data and
trust controls. The large raw gradient of an unfrozen parent is not by
itself evidence that unfreezing will generalize.

Next test smooth spherical-harmonic × torus-Fourier correction features in
an increasing-degree ladder. The model remains a globally defined ambient
field; this does not reintroduce an atlas. Spectral truncation and smooth
partition-of-unity enrichment must be selected using held-out performance.
No unsmoothed nearest-neighbor/RBF patching or discontinuous per-chart J.

Possible further escalation is a Sobolev-accurate distillation/rebase of
the SAME learned field to a better-conditioned architecture. Require both
J and first-derivative agreement before evaluating Nijenhuis. A small error
in J alone is insufficient. Rebase hashes and original growth anchors persist.

If these approaches reach a resolved stationary floor across capacity and
seed tests, investigate additional learned ACS homotopy classes explicitly.
Measure class diagnostics and audit each new base. More exponential stages
on the current base are not a substitute. This is a separate research stage,
not an automatic promise that another class contains a CS.

## 3. Numerical gates, tests and autonomous HPC workflow

### Required tests before submission

- E E^T = I6, E^T E = P, correct noncommuting-frame bracket signs, random
  points, sphere poles, circle seams and highly conditioned saved points.
- Ambient/frame agreement for J, signed N, norm, all parameter gradients,
  and GN JVP/VJP products. Test ordinary, correction, and unfrozen models.
- Frozen-cache equivalence after parameter updates and at new cache points;
  invalidation tests for point/reference changes and parent unfreezing.
- Adversarial searches use live second derivatives. Check against finite
  differences; deliberately test that a stale first-jet shortcut is rejected.
- Nonconstant integrable T6 positive control and nonintegrable controls;
  known CS controls such as S3 × S1 for end-to-end optimizer tests. A control
  passing is not a convergence theorem for S4 × T2.
- Importance weights reproduce the uniform finite-pool expectation and
  gradients; no weighting of the separately reported uniform validation.
- Signed-component GN derivative checks; exact retraction; trust-region
  rejection restores parameters AND optimizer/RNG state.
- Resume versus uninterrupted equivalence, wall/signal interruption,
  corrupted/missing artifact detection, and relocated-path resolution by
  hash rather than editing hash-bound metadata.
- Verify closure loss, actual uniform RMS, maxima, and per-factor tensor
  energies independently. Do not treat a sum of regularizers as the target.
- PBS dependency and failure-injection tests; memory and runtime smoke on
  the actual single-thread serial queue before scheduling larger work.

### Data and promotion policy

Maintain three separate roles: training; development/selection; final audit.
Previously audited witnesses may become training data only in the next
declared campaign. Repeatedly inspected development data are not independent
final evidence. All candidates in a comparison use the same expanded probes.

Keep up to two nondominated lineages: lower RMS and lower maximum. Never
average their J matrices: that generally destroys J^2=-P. A temporary
RMS/tail tradeoff can be retained without satisfying simultaneous improvement
at every update. Still enforce algebra, tangency and physical growth guards.

At the current RMS near 0.28, focus the low-RMS lineage on reducing the bulk
residual while imposing a maximum-growth safeguard, rather than selecting
solely on a modest improvement in the worst found. Once RMS is near 0.03,
increase emphasis on the tail. Compare a fixed-inner-solve top-tail objective
with an augmented-Lagrangian witness constraint r(x)^2 <= 1e-4. Adapt those
weights only between solves, keep uniform RMS explicit, and require new
searches: satisfying the old witness constraints is not global control.

Suggested engineering promotion rule: at least 10% improvement in one
primary metric over a bounded pilot allocation, with the other no more than
5% worse, or a comparable documented joint gain. Keep the two-sided target
unchanged. These percentages are budget decisions, not mathematical laws;
calibrate them from matched pilot variability before fixing the campaign.

### Dependency graph

    preflight + HPC numerical/timing smoke
                    |
       matched bounded optimizer/sampling pilots
                    |
       compare on common fresh selection probes
                    |
        retain <= 2 admissible Pareto lineages
                    |
       resumable blocks + pool/search refresh
                    |
     progress? ---- no --> bounded representation probe --> review if stalled
         |
        yes
         |
       target margin reached --> freeze artifacts --> independent audit --> report

Use existing serial/10 GB/4 h allocations, one thread, a <=170-minute work
budget, and saved accepted state. Start with 90-minute method pilots; expand
only measured winners into resumable blocks. Submission dependencies can use
`afterany` with hash-bound success/review receipts, so a numerical review
terminal state propagates safely and a failure cannot feed a bad checkpoint.
No modification of code/config while a campaign is running.

An explicit first budget is eight 90-minute pilots: uniform Adam,
importance-sampled Adam, trust-controlled L-BFGS, and preconditioned GN, each
from the low-RMS and low-tail parent. That is at most 12 single-core training
hours, plus preparation/comparison. Within the Adam pilots, use a short
declared LR screen before the scored continuation, rather than silently
multiplying the number of long branches. Initially authorize at most two
170-minute continuation blocks per retained lineage (at most 11.3 further
training hours for two lineages) before the next explicit progress gate.
Budgets can be extended only on measured progress, not a forecast that 0.01
must eventually be reached.

Provisional milestones, NOT forecasts:

| Gate | Uniform RMS | Maximum found |
|---|---:|---:|
| First material improvement | 0.20 | 0.45 |
| Clear progress toward low residual | 0.10 | 0.20 |
| Precision phase | 0.03 | 0.05 |
| Margin before final audit | 0.007 | 0.008 |
| Agreed final sampled acceptance | 0.01 | 0.01 |

Do not automatically spend many more allocations when two comparable blocks
fail to improve either frontier substantially. Trigger the declared
representation/sampling diagnostic, then request review if it also fails.
No plateau licenses a nonexistence claim.

The final audit retains at least the existing 100,000 fresh uniform points,
structured checks, and three independent 16,384-pool / 512-start / 200-step
adversarial searches. Add multiscale neighborhoods of new maxima and
sensitivity-selected points. Cross-check signed tensors in local coordinates
and with finite differences; report J, derivative and conditioning growth.
All acceptance decisions refer to the unmodified round metric and one frozen
checkpoint. Maximum found remains a lower bound on the true supremum.

## 4. Literature and limits of the strategy

- [Rathore et al., ICML 2024](https://proceedings.mlr.press/v235/rathore24a.html)
  motivates testing preconditioning and second-order PINN optimization. Its
  PDE benchmark successes do not guarantee convergence of this nonlinear
  geometric search. The local failed L-BFGS trials are more directly relevant
  to choosing this campaign's data and trust controls.
- [Wu et al.](https://arxiv.org/abs/2207.10289) supports testing residual-adaptive
  collocation. Here the per-point gradient diagnostic additionally motivates
  sensitivity-aware coverage; that extension is our experimental proposal.
- [Parton](https://arxiv.org/abs/math/0009156) supplies explicit parallelizations
  of parallelizable sphere products. Using such a frame to evaluate a learned
  tensor is not using a known CS or replacing the learned ACS.
- [Euh–Sekigawa](https://arxiv.org/html/1301.6835) recalls Calabi's question for
  V2 × S4. Its restricted orthogonal results are not general nonexistence
  theorems for the present search. The literature checked does not provide
  a verified reference CS on S4 × T2 or a convergence guarantee to one.
- [Granja–Milivojevic](https://arxiv.org/abs/2207.12946) studies ACS topology;
  results assuming vanishing b1 do not apply directly here (b1=2). Preservation
  of homotopy class under our exponential continuation follows directly from
  scaling each generator continuously to zero.

The strongest justified claim is that the plan targets measured computational
and optimization problems while preserving the geometry and evidence
standards. No honest analysis of these checkpoints establishes that any
particular strategy will attain 0.01 on this manifold.

## 5. Local experiment artifacts

All numerical outputs and prototype scripts are ignored scratch artifacts,
not production training changes. The main original files remain untouched.

- `scratch/analyse_s607.py`: hashes, baseline, gradient/caching tests and
  small fixed-data L-BFGS experiments.
- `scratch/frame_s607.py`: equivalent six-dimensional evaluation, full
  parameter-gradient checks and repeated timing comparison.
- `scratch/leverage_s607.py`: per-point gradient concentration and finite-pool
  importance-sampling variance calculation.
- `scratch/galerkin_s607.py`: verified function-space linearization and
  regularized coarse-feature probes.
- `scratch/topology_s607.py`: converging torus-slice Chern-form diagnostic.
- `scratch/largepilot_s607.py`: larger-data L-BFGS / field-update penalty probe.
- `scratch/verify_s607.py`: parameter and live-input finite differences,
  coordinate cross-checks and frame checks at high-J-norm learned points.
- `scratch/resolved_gn_s607.py`: streamed 360-coefficient GN with 8,192
  scrambled-Sobol training points, 2,048 separate uniform development points,
  and regularization of the field variation and its covariant derivative.
- `scratch/cs_s607_analysis_0wj_or15/*.json`: numerical results.

These prototypes contain experiment-specific local paths and deliberately
small development tests. They must NOT be submitted as the next HPC campaign.
