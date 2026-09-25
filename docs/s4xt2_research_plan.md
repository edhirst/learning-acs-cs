# S4 × T2: assessment and plan for numerical ACS and CS

Date: 2026-09-25. Code reviewed: `1383f80` on `no_overlaps`.

This is a research and implementation plan, not a record of a successful CS
search. The objective is a globally defined, numerically validated almost
complex structure (ACS), followed by an integrability search on
`s4xs1xs1 = S4 × S1 × S1`. Keep the existing nine-dimensional ambient embedding
and round product metric for evaluation.

Priority update: learn the ACS first, from random initialization in the existing
ambient formulation. The detailed [ambient ACS learning and HPC plan](s4xt2_ambient_acs_plan.md)
is now the active execution plan. CS experiments are deferred until that
learned ACS passes independent validation. This update changes the plan only;
the proposed trainer changes and submission package are not yet implemented.

## 1. Recommended direction

Learn K = P M_theta P using ambient points in R9 and general 9 × 9 network
outputs, matching the other manifold ACS runs. Add three-factor polynomial
features, a genuine polynomial-only stage, ACS-adaptive sampling, and a loss
sensitive to large residuals. Validate these changes before controlled pilots
and longer, resumable HPC runs.

The analytic ACS and its frame are test fixtures and offline diagnostics only.
They will not supply training labels, initial solution coefficients, or a
replacement for learning. Known existence and a coefficient-level
representability test establish an attainable target; finding it from random
initialization is the optimization task.

Continue seed 10 in a bounded comparison while testing fresh random starts.
Completion requires an accurate raw learned field, stable polar correction,
independent validation, and a reproducible artifact. After success, return to
the CS investigation and its topological and representation choices below.

## 2. What the evidence establishes

### Existing HPC experiments

The user supplied these final summaries; the remote checkpoints were not
available for this review.

| Seed | Dense mean ACS residual | Dense maximum | Adversarial maximum found |
| --- | ---: | ---: | ---: |
| 7 | 0.049803 | 1.415174 | 1.766203 |
| 8 | 0.038817 | 1.431555 | 1.725575 |
| 9 | 0.044362 | 1.427914 | 1.696862 |
| 10 | 0.045347 | 1.450490 | 1.583377 |
| 11 | 0.046574 | 1.447873 | 1.759135 |

All completed 10,000 refinement steps without the wall-time stop. None passed
the current raw-residual gate of 1. The CS stage correctly refused those bases.
Seed 10 has the smallest adversarial maximum; seed 8 has the smallest mean.
These are different objectives and should not be conflated.

The gap between means and maxima is consistent with concentrated defects. A
residual near sqrt(2) is consistent with failure on a two-dimensional tangent
subspace, but does not diagnose the defect's rank, location, orientation, or
topology. Inspect tangent singular values and eigenvalues at saved maxima
before making that claim. At K = 0 on the entire six-dimensional tangent space,
the residual would be sqrt(6), not sqrt(2).

### Code findings

- The ambient projector and arbitrary-factor geometry are suitable for this
  manifold. A full ambient output is general, but its 81 entries include
  directions discarded by projection; a global frame needs only 36 entries.
- `extra/test_ambient_products.py` already contains an exact ACS constructed
  from a global frame. Production registration does not expose it.
- `Geometry(..., 'frame')` currently reports that this case is not
  parallelizable. This is an incorrect mathematical interpretation of a
  missing frame registration.
- ACS refinement uses uniform sampling. `--tail` and `--rar` affect CS only.
- The polynomial ladder ends at degree two. Pairwise bilinear support does
  not supply three-factor products.
- The nominal polynomial-only phase still evaluates a small random MLP
  contribution. Isolate the stage explicitly before interpreting polynomial
  representability or convergence experiments.
- The current tail builder requests groups of eight ascent starts but retains
  only one point per group. `--tail-hard 24` therefore produces three hard
  points, not 24.
- `--init-from` is designed around CS. For ACS it does not reconstruct the
  saved curriculum feature class and mark the ladder completed. Exact resume
  also restores the old learning-rate schedule; extending its step limit is
  not a clean continuation with a fresh schedule.
- Best-model selection and L-BFGS rollback are based on mean losses, which can
  prefer a model with a worse maximum residual.
- The saved step counter after L-BFGS can be set to the planned total even
  after an early training stop. Actual optimizer work needs separate counters.
- Long polynomial attempts lack inner wall-time guards and mid-attempt
  checkpoints. Routine ACS dense evaluation also computes costly Nijenhuis
  derivatives that are unnecessary for monitoring the algebraic objective.
- Adversarial maxima are lower bounds on the true supremum. The result key
  `certified_exact_acs` overstates what a finite search proves.
- A tiny chart-overlap discrepancy checks representation consistency. It is
  not the ACS residual and does not establish integrability.

### Checks performed during this review

The existing focused suite passed: 8 tests. PyTorch reported a batching fallback
for matrix exponentials, which warrants profiling rather than a correctness
conclusion.

An additional reproducible diagnostic used float64, one PyTorch thread, seed
2026, and 512 sampled points. For the existing exact reference:

- maximum raw Frobenius norm of J² + P: `1.2591548905310372e-15`;
- maximum tangency residual: `6.82368286116403e-16`;
- Nijenhuis norms on the first 64 points: mean `2.4514222857566024`,
  minimum `1.0379462454703818`, maximum `2.828426945692029`;
- agreement with a derived degree-three raw ambient representation after
  projection: maximum Frobenius difference `7.62856359580219e-16`.

These are local numerical checks of the existing reference and implementation.
They do not establish training convergence or a CS.

## 3. Literature and mathematical constraints

### ACS existence is settled

[Parton, *Explicit parallelizations on products of spheres*](https://arxiv.org/abs/math/0009156)
gives an orthonormal frame for S^m × S1 in Section 2. Adding the second circle
gives a frame for the present six-manifold. This is the relevant constructive
result, stronger for our purpose than merely knowing that an ACS exists.

In the repository's row-frame convention, put p = (x,u,v), where x is a unit
vector in R5 and u,v are unit vectors in R2. Let t_u = (-u_2,u_1) and similarly
t_v. Define

    E_a = (e_a - x_a x, x_a t_u, 0),  a = 1,...,5,
    E_6 = (0, 0, t_v).

With E the 6 × 9 matrix of these rows, direct multiplication gives
EE^T = I6 and E^T E = P. Consequently

    J_ref = E^T J0 E,     J0 = diag(J2,J2,J2),
    J2 = [[0,-1],[1,0]],

satisfies J_ref² = -P everywhere. No inverse square root is needed for this
reference. Its nonzero Nijenhuis tensor is the remaining equation to address.

### CS existence must not be assumed

[Euh and Sekigawa (2013)](https://arxiv.org/abs/1301.6835) explicitly recall
Calabi's question about complex structures on V² × S4 for a closed oriented
surface V. Their orthogonal nonexistence result concerns the specified round
even-sphere products; it is not a nonexistence theorem for arbitrary complex
structures on S4 × T2.

The search conducted for this plan did not identify a verified construction
of a CS on S4 × T2 or a theorem excluding every CS on it. Treat its CS search
as a research problem with an uncertain outcome. An ACS construction is not
an existence theorem for a CS.

The repository also treats recent S6/S2 × S4 developments too categorically.
The abstract of [Liang–Liu–Zhang (2026)](https://arxiv.org/abs/2608.29651)
explicitly makes its construction conditional on a recently proposed result.
The separate [Guo–Fang–Lu preprint](https://arxiv.org/abs/2609.05504) also
concerns S2 × S4. Neither supplies a directly usable S4 × T2 reference. Those
developments should be labelled by their stated assumptions and status, and
should not be prerequisites for this campaign.

### The target cannot be Kähler

For this manifold, H²(M;R) is generated by the torus area class a and a² = 0.
Thus every closed two-form has [omega]³ = 0. A symplectic form on a compact
six-manifold would have nonzero integral of omega³. This elementary cohomology
argument rules out symplectic, and therefore Kähler, structures.

Do not impose a closed fundamental two-form, a Kähler loss, or a parallel-J
condition. Measuring Nijenhuis with the fixed round product metric does not
impose any of those conditions. General nonorthogonal ACS must remain allowed.

### Exponential conjugation has a global restriction

For smooth matrix fields B_s,

    J_t = exp(t B_s)...exp(t B_1) J_base exp(-t B_1)...exp(-t B_s)

is a homotopy through ACS. Every finite-stage model stays in the base's
homotopy class and preserves its Chern classes. Pointwise transitivity of
GL+(6) does not prove that two globally smooth exponentials cover all fields
of a chosen orientation. Correct that claim in the existing mathematics note.

[Granja–Milivojević (2022)](https://arxiv.org/abs/2207.12946) describes the
topology of spaces of ACS using twistor bundles and Chern classes. Its results
with a vanishing-first-Betti-number hypothesis must not be applied to this
manifold, which has b1 = 2. The conjugation restriction above also follows
directly from the displayed homotopy.

Here the explicit constant-frame reference has c1 = 0. Since TM is trivial,
c1 mod 2 = 0. Write possible first Chern classes as 2k a. Moreover p1 = 0
and a² = 0 imply c2 = 0, and chi(M) = 0 implies c3 = 0. Different values of k
therefore deserve explicit exploration. Increasing the random-seed count or
the number of exponentials cannot change k.

### Optimization literature informs experiments, not guarantees

[Wu et al. (2023)](https://arxiv.org/abs/2207.10289) compare uniform and
residual-adaptive sampling on PDE benchmarks and support testing adaptive
collocation. They do not prove convergence of a Nijenhuis search.

[Rathore et al., ICML 2024](https://proceedings.mlr.press/v235/rathore24a.html)
study ill-conditioning in PINN losses, motivate Adam followed by L-BFGS, and
investigate a Newton-CG refinement. This supports retaining the current
optimizer sequence while improving initialization, validation, and objective
design first. A second-order method is a later measured escalation.

[Newlander–Nirenberg](https://www.mathnet.ru/eng/mat750) characterizes exact
integrability by vanishing Nijenhuis. A small numerical residual alone does
not give a quantitative distance to an integrable structure.

## 4. Assessment of strategies considered so far

| Strategy | Assessment | Decision |
| --- | --- | --- |
| Ambient embedding and P M P | Smooth global fields with automatic tangency; appropriate general framework | Retain for all ACS training, exports, and validation |
| More copies of the existing ACS run | Five runs have the same mean/max gap | Stop unchanged restarts |
| Polynomial curriculum through degree two | Useful bias but misses a simple three-factor reference | Add targeted cubic interactions; verify representability |
| Large Dirichlet penalty | Can prefer shrinking K and penalize valid nonparallel structures | Use only as an annealed ablation; do not make it the main repair |
| Tail/adversarial ACS continuation | Directly addresses observed large residuals, but may move or sharpen defects | Compare uniform and adaptive continuation from seed 10 alongside fresh starts |
| Polar correction | Valid under its spectral hypotheses; not a remedy for globally singular fields | Retain for learned ACS, with margin and independent numerical checks |
| Exact reference from a global frame | Supplies an independent representation and validation fixture | Tests only; cannot substitute for learning the ACS |
| Learning 6 × 6 coefficients in the frame | Removes normal-direction redundancy without excluding tangent tensors | Deferred alternative for later CS work; retain full ambient ACS output now |
| One versus two exponentials | One is a useful local model; two can increase pointwise flexibility; neither changes topology | Deferred to CS pilots |
| Strong orthogonality or Kähler constraints | Restrict the target without justification; Kähler is impossible here | Do not impose |
| Mean Nijenhuis minimization alone | Can conceal a small region of large residual | Add a controlled tail objective and stronger validation |
| Adam followed by L-BFGS | Reasonable baseline supported by PDE experiments | Keep, with fixed closure samples and improved rollback |

### Why the missing cubic features matter

For the particular J_ref above, the block from the second circle to the first
is -x_5 t_u t_v^T. It is odd under each of x -> -x, u -> -u, and v -> -v.
The projectors are even under these independent antipodal maps. A raw
polynomial of total degree at most two cannot produce a component odd in all
three factors after projection. Thus the current polynomial-only model cannot
represent this particular exact reference. This does not prove that it cannot
represent any ACS, or that the subsequent MLP lacks enough capacity.

A sparse extension is enough to include the reference: constants, linear
terms, existing pairwise terms, and all 5 × 2 × 2 = 20 products x_i u_j v_k.
To verify this, let Q be the upper-left 5 × 5 block of J0 and b = e_5. A raw
skew matrix M with blocks

    M_xx = Q,
    M_xu = (Q x) t_u^T,
    M_xv = -b t_v^T,
    M_uv = -x_5 t_u t_v^T,
    M_uu = M_vv = 0,

and negative-transpose lower blocks satisfies P M P = J_ref. This identity
was checked numerically during this review. Use it as a representation test,
not as a supervised target in the residual-only experiment.

## 5. Execution plan and acceptance gates

### Phases A–B — learn and independently validate the ambient ACS first

Follow the [dedicated ACS plan](s4xt2_ambient_acs_plan.md). Its sequence is:

1. Verify generic triple features against an exact coefficient fixture and
   isolate polynomial training from the disabled MLP.
2. Test ACS losses, gradients, adaptive sampling, narrow-defect detection,
   polar correction, and other-manifold regressions.
3. Repair and test candidate selection, restart/continuation, phase checkpoints,
   wall-time handling, and inexpensive ACS-only evaluation.
4. Run local controls and HPC smoke jobs, then up to 12 matched pilots and
   two seed-10 continuation comparisons.
5. Promote selected recipes to fresh long runs: up to 20k polynomial steps
   and 40k refinement steps, split across resumable allocations.
6. Independently audit candidates, targeting raw RMS <= 1e-4, maximum found
   <= 1e-2, and corrected algebraic error <= 1e-10, using floor-free norms.
7. Preserve the first successful learned artifact and establish repeatability
   with additional random starts. Return to CS only after this milestone.

The following phases are retained as future research options. They do not
authorize or schedule CS training during the ACS campaign. In particular,
the frame-based CS representation is a proposal to revisit, not an active
change to the ambient ACS architecture.

### Deferred Phase C — construct and identify topological starting classes

1. Extend the Chern diagnostic from S2 cycles to the torus cycle
   {x0} × S1 × S1, using periodic quadrature and resolution doubling.
2. Fix the eigenspace/sign convention: (P - iJ)/2 projects onto the +i
   eigenspace; the existing (P + iJ)/2 expression selects the -i eigenspace.
   Calibrate against a known oriented example and label the convention.
3. Confirm c1 = 0 for the constant-frame reference and invariance under
   smooth conjugation.
4. Build reference ACS in sectors with c1 on the torus equal to 0, +2, and -2;
   expand to +/-4 if pilots justify it. One constructive route is a smooth
   degree-k map T2 -> S2, the quaternionic S2 family of ACS on R4, and a
   fixed J2 on the remaining R2, transported by E. Verify degree and Chern
   number rather than inferring them from a network seed.
5. Give each CS checkpoint an immutable base and sector identifier. A
   continuation must check both, not just dimensions and architecture.

Gate: exact ACS identities for each base, convergent integer Chern integrals,
and invariance under representative deformations. Do not claim a finite list
of sectors exhausts all possibilities.

### Deferred Phase D — a better conditioned CS representation

Use ambient points p in R9 as inputs, with full dependence on both circles.
Circle coordinates already make periodicity automatic. Predict general
6 × 6 matrices B1(p), B2(p) and form

    G = exp(B2) exp(B1),
    A = G A_base G^{-1},
    J = E^T A E.

Initially reuse ambient Nijenhuis evaluation for a straightforward bridge.
After tensor-level agreement, benchmark the frame evaluator with analytic
brackets and six tangent directional derivatives. Profile both forward loss
and training gradients; speedups are not implied solely by matrix size.

Compare one and two stages with the same backbone and wall-time budget. Start
at the base or a small perturbation (for example scales 0, 0.02, 0.1), measuring
initial Nijenhuis and condition numbers. The present O(1) initialization can
spend substantial computation undoing a very distorted starting field.

Keep general GL deformations. Orthogonal conjugations are an optional
restricted ablation. Later, if redundant parameters cause trouble, compare an
18-real-component local deformation model: at an ACS A, infinitesimal
variations anticommute with A, and the component of B commuting with A has no
first-order effect. Treat such a reduced model as local and use re-centering;
do not assume it is one globally exhaustive coordinate chart.

Monitor ||A||, ||dA||, the singular values and condition number of G, and the
actual ACS residual throughout. If R_base = A_base² + I, then

    R_final = G R_base G^{-1},
    ||R_final|| <= condition_number(G) ||R_base||.

Thus algebraic error amplification is controlled by conditioning, not by a
blanket claim that the two-stage design is stable or unstable. Nijenhuis also
depends on derivatives of G and the base. Those derivatives allow the model
to change integrability; freezing a base does not freeze its Nijenhuis tensor.

### Deferred Phase E — CS training with explicit decision points

1. Prove the training implementation can reduce a deliberately perturbed,
   known integrable structure on a six-dimensional control such as S3 × S3
   or S5 × S1. Validate the control independently first. Two-dimensional
   controls alone are inadequate because their ACS are automatically
   integrable.
2. Start S4 × T2 pilots using a small smooth network, e.g. three layers of
   width 128, then compare with the existing width-256 network only when
   useful. Preserve full input dependence and cross-factor interactions.
3. Train in float64 with Adam, a scheduled learning rate, and moderate gradient
   clipping. Use a mean-plus-tail Nijenhuis objective, with a uniform batch
   component retained throughout.
4. Add soft conditioning or derivative penalties only when diagnostics show
   instability. Anneal and ablate them. A fixed strong smoothness penalty
   changes the target and could bias against the desired structure.
5. Follow with L-BFGS on a fixed collocation set within each solve. Refresh
   between solves. Reject nonfinite updates explicitly and use both mean and
   tail validation for rollback. Save the best accepted model separately
   from the most recent model and optimizer state.
6. On a plateau, inspect gradient norms, residual distribution, conditioning,
   and representation saturation. Compare Fourier features on the two circles
   and low-degree spherical features before indiscriminately adding width.
7. If the residual is smooth but optimization stalls, trial damped
   Gauss–Newton/Levenberg–Marquardt or Newton-CG with Jacobian-vector products.
   Reuse any repository optimizer components only after adapting and testing
   them against the current exact-ACS residual. Compare reduction per wall
   time and memory, including the cost of derivatives.
8. Do not infer nonexistence from a plateau. Expand topological sectors,
   parameterizations, or regularization choices in response to diagnostics.

A change of frame is harmless when coefficients remain general. Imposing
torus invariance, a product splitting, or constant coefficients throughout
would be additional restrictions; use those only as labelled ablations.

### Deferred Phase F — independent CS validation and evidence of convergence

Keep the metric fixed to the unit round product. A trainable metric must not
be allowed to make a nonzero residual look small by rescaling norms.

Report separately:

- ACS identity and tangency errors;
- Nijenhuis mean norm, RMS norm, median, 95th/99th/99.9th percentiles, and
  maximum found;
- independent adversarial searches with saved maximizers;
- norms of J and its first and second derivatives, plus conditioning;
- sector/Chern measurements, numerical precision, and actual compute cost.

Use at least 100,000 fresh points for finalist assessment, in chunks, plus
structured sphere/torus slices and targeted neighborhoods of old maxima.
If the derivative cost is high, benchmark and budget this explicitly rather
than quietly reducing validation to the training batch. Increase independent
adversarial starts and use multiple search seeds. Reusing training hard points
alone is not an independent audit.

Provisional engineering milestones under the fixed metric:

| Milestone | Criterion |
| --- | --- |
| Valid ACS | Algebraic residual <= 1e-10; tangency, periodicity, reload, and derivative checks pass |
| Useful CS pilot | At least 10x improvement in both RMS Nijenhuis and maximum found, without uncontrolled norm/derivative growth |
| Strong numerical candidate | RMS Nijenhuis <= 1e-5 and maximum found <= 1e-4, supported by independent tests |
| Convergence evidence | Residuals continue falling as representation, optimization, and validation resolution increase, with bounded conditioning and controlled regularity |

These thresholds are proposed targets, not guarantees of existence or a
theorem converting residual to distance from a CS. Calibrate attainable
precision and cost on positive controls first; record any threshold change.

For a promising finalist, fit local complex coordinate functions z on a
cover using dz o J = i dz. Normalize their derivatives to prevent constant
solutions, check full real Jacobian rank, and examine approximately holomorphic
transition maps. These would be additional tests, not an automatic consequence
of a small sampled Nijenhuis tensor. More than the original eight large
stereographic patches may be needed for well-conditioned local coordinates.

A rigorous result would require further work: a verified uniform residual
bound and an applicable analytic correction/existence argument, or an explicit
integrable construction. Finite adversarial search is not such an argument.

### Deferred Phase G — later CS campaign and operational considerations

The active ACS submission sequence and budgets are specified in the
[ACS HPC plan](s4xt2_ambient_acs_plan.md). This section concerns the later CS
campaign, after independent acceptance of the learned ACS.

Benchmark one training/evaluation block on the actual Lovelace environment
before assigning long jobs. Preserve the accepted serial resource settings
until queue information supports a change. Do not infer GPU suitability from
availability; compare float64 derivative throughput and memory first.

Possible later CS campaign, to reassess after ACS success:

| Stage | Runs | Purpose |
| --- | ---: | --- |
| Learned-base import and geometry validation | 1 | Reload and revalidate the successful ACS artifact |
| Nontrivial integrable controls | 2–3 | Establish solver and validation performance |
| CS pilots | 6 | Three initial Chern sectors × two initialization seeds |
| Stage-count comparison | 2 | One versus two exponentials in the most promising sector |
| Longer continuation | Best 2 | Fresh schedules; expanded validation; preserve originals |
| Independent final assessment | Each finalist | No parameter updates; final evidence report |

Give pilots a measured wall-time budget, with an initial 2,000–5,000-step
cap if that fits. Select continuations by progress per wall time and large-
residual behavior. Do not launch another unconditional 80,000-step grid.

Keep short site submission scripts in gitignored `scratch/`. Keep the actual
experiment settings and run documentation tracked elsewhere, using the
ambient CLI or a small validated configuration layer; legacy atlas YAML files
do not control this trainer. Include an explicit preflight for case, source,
revision, base, output directory, and environment.

Use new output directories for new experiments. Save self-contained base
recipes or frozen base artifacts with hashes. Persist adversarial replay
points, optimizer/scheduler state, RNG state, and run counters. Distinguish
exact restart from fresh-schedule continuation, and make interruption tests
part of the implementation acceptance criteria.

Afterok dependencies only check exit status. The CS script must also verify a
machine-readable ACS acceptance status and finite corrected residual. Write
an explicit incomplete/rejected status, and make the submission wrapper's
exit policy match the intended dependency behavior. Check that qsub actually
returned a job identifier before submitting a dependent job.

Stop an individual run on nonfinite fields, invalid ACS error, unstable
conditioning, or exhausted budget, saving its state and reason. Stop expanding
an unchanged strategy when repeated matched pilots cease improving. Failure
of all tested strategies is a bounded negative computational result, not a
topological nonexistence proof.

## 6. Proposed implementation order

1. Generic triple-product features and genuine polynomial isolation in the
   ambient trainer; representation, gradient, and compatibility tests.
2. ACS-specific adaptive sampling, tail objective, economical diagnostics,
   and reliable model selection.
3. Reproducible restart/continuation, phase checkpoints, truthful counters,
   and bounded wall-time handling in `extra/campaign_unit.py`.
4. Tracked ambient run configurations, an independent ACS validator, and short
   site scripts in ignored `scratch/`; test submission and interruption paths.
5. Local controls, HPC smokes, matched pilots, and selected longer ACS runs.
6. Final independent ACS audits, repeatability runs, and learned-artifact
   handoff with complete provenance.
7. Subsequently reconsider Chern diagnostics, CS representations, integrable
   controls, and CS pilots from the accepted learned base.

The first deliverable is a learned, independently validated ambient ACS and a
repeatable HPC recipe. An analytic reference export does not satisfy it.
CS learning remains a later task after that result is established.
