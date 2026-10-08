# S⁴ × T²: a scheme for a numerical complex structure

2026-09-30. Research plan only: no code, configs or runs were changed.
It replaces the incremental continuation plans (`s4xt2_cs_long/precision/refine/
s607_analysis`) as the strategy for the complex-structure (CS) search. The ACS
pipeline and audit infrastructure are kept.

---

## 0. Summary

**State of play.** The ACS problem is solved: the s307 field audits at raw RMS
9e-5, and an exact analytic ACS exists anyway (the Parton frame). The CS search
has plateaued: uniform Nijenhuis RMS is about 0.28 and the largest value found is
about 0.63. The last 43k updates improved RMS by 0.7%. The repo records show
that this plateau tells us almost nothing about whether a CS exists:

1. **The optimiser is noise-dominated at a stationary point.** At batch size 128,
   pairwise gradient cosines average −0.017. The norm of the 8-batch mean
   gradient, 1.70, matches the pure-noise prediction of about 5/√8 ≈ 1.8. The
   population gradient is therefore below the sampling noise even at 1,024
   points. Continuing small-batch Adam at LR ≤ 1e-5 cannot help.
2. **The fields have drifted far from orthogonal.** The Frobenius norm of J has
   quantiles 4.1 / 6.1 / 11.7 / 21.6. An orthogonal J has √6 ≈ 2.45. Gradient
   energy sits on a few high-|J| points: the top 1% of points carry 49% of it.
   §2.1 shows this is the signature that a mean-square objective rewards.
3. **Second-order steps overfit the collocation set.** In those trials the
   training loss fell, but held-out RMS rose from 0.277 to 0.380 and the
   held-out maximum from 0.45 to 3.7. This points to under-resolved 6D
   quadrature, not a capacity limit.
4. **No nontrivial CS has ever been recovered by this code base.** That holds in
   4D and 6D, from any non-trivial start. The chart-era S³×S³, S⁴×S², S⁶ and
   even S⁴ runs reached per-patch Nijenhuis of 1e-5 to 1e-2. But their overlap
   discrepancies were 2.3–3.3, so the patches never glued; S⁴ has no ACS at all.
   With gluing enforced (overlap weight 1000), Nijenhuis rose to 3.5–27. S²×S²
   stalled at 0.07, S³×S¹ at 0.09–0.14, and the Gauss–Newton S³×S¹ run at 0.61.
   The only CS control passed is the T⁶ pullback, which is trivial.

**The scheme in one paragraph.** (i) Redefine success. A mean-square Nijenhuis
loss can be driven to zero on any almost-complex 6-manifold by concentrating the
defect (§2.1), so it is not evidence. The deliverable is a convergent sequence
in which sup|N| → 0 while sup|J|, sup|∇J| and sup|∇²J| stay bounded; by
compactness that yields an honest CS (§2.2). (ii) Calibrate before spending
compute. Positive controls, the now-known-to-be-complex S⁶ and S²×S⁴, and a
4D negative control with a known answer (§4, Phase 0). (iii) Climb a symmetry
ladder, ODE → 2D → 4D → 5D → 6D. The analytic reference ACS is
U(2)×T²-invariant, so it lies in every rung. The T²-invariant rung is an exact
4D problem on S⁴ (§2.5, Phase 1). (iv) Use a better-posed formulation:
Cayley/Beltrami parametrisation with re-centring, a supercritical objective,
gauge awareness, function-space Gauss–Newton with Laplacian preconditioning,
batches 10²–10³ times larger, and GPUs (§3). (v) Polish with a spectral
discretisation, check the linearised operator, and identify the structure
(Phase 3). (vi) In parallel, a theory track: sector classification and
construction ideas from the new S⁶ / S²×S⁴ results.

---

## 1. Evidence base (from this repo)

| Item | Evidence | Implication |
|---|---|---|
| ACS | s307: raw RMS 9.1e-5, max found 2.1e-4, corrected identity 3e-15. Parton frame gives an exact ACS `J_ref = Eᵀ J0 E` | ACS is not the bottleneck. Its homotopy class (k = 0) is all that matters for CS |
| CS trajectory | RMS 1.54 → 0.46 → 0.44 → 0.31 → 0.289 → 0.287; max found 2.1 → 0.63 | Plateau |
| Field geometry (s607 low) | ‖J‖_F quantiles 4.1 / 6.1 / 11.7 / 21.6; ‖∇J‖ ≈ 7–8 at sampled points | Strongly non-orthogonal everywhere; distortion growth over the lineage (3.4 → 10.7 → 21.6) |
| Gradient statistics | cos ≈ −0.02 between 128-point batches; top 1% of points carry 49% of gradient energy; Spearman(‖g‖, ‖J‖) = 0.66 | Noise floor; gradient sits in distorted regions |
| L-BFGS / GN probes | train ↓, held-out RMS and max ↑ (max up to 3.7); a 1,980-coefficient Galerkin step made held-out RMS 0.27 → 1.05 | 6D quadrature is under-resolved relative to model degrees of freedom |
| Controls | Only the T⁶ pullback passes; see summary item 4 | Method not yet validated on any nontrivial CS |
| Compute | 1 CPU thread, batch 32, 10 GB / 4 h serial jobs | Orders of magnitude below what 6D derivative-heavy least squares needs |
| Topology | c₁ torus pairings ≈ 0 for base and CS fields | Search confined to sector k = 0 (§2.3) |

---

## 2. Mathematical facts that should shape the design

### 2.1 Mean-square Nijenhuis cannot be evidence (concentration)

**Claim.** In the ACS class of `J_ref` (and, by Landweber, in every class), the
infimum of ∫|N_J|ᵖ_g dvol_g over J in the class is 0 for every p < 6. This holds
whether or not a CS exists.

*Sketch.*
1. M° = M ∖ {pt} is open. The Parton frame trivialises (TM, J_ref) as a complex
   bundle. The Gromov–Phillips h-principle then gives an immersion f: M° → ℂ³
   whose differential is homotopic to that trivialisation. So J₁ = f*J_std is
   integrable on M° and homotopic to J_ref. (Landweber 1974 gives the general
   statement for open almost complex manifolds.)
2. Keep J₁ on M ∖ B for a small ball B, and extend it arbitrarily inside B. The
   only obstruction lies in π₅(ℂP³) = π₅(S⁷) = 0. The result is J₂, integrable
   off B and in the same class, since the class is detected on M ∖ B.
3. Pull J₂ back by a diffeomorphism φ_r that dilates B_r onto B and stretches a
   collar radially. Then N is supported in B_r with |N| ~ 1/r, so
   ∫|N|ᵖ ~ r^{6−p} → 0.
4. The cost is |∇J| ~ 1/r in B_r and slowly growing distortion |J|_g on the
   collar.

**Consequences.**
- **What does not count.** A falling mean-square loss, or "RMS ≤ 0.01", is not
  evidence of existence.
- **Critical exponent.** The critical exponent is p = dim M = 6. Objectives must
  be supercritical (p > 6, or a sup), or must carry explicit C¹ control.
- **What to monitor.** Distortion growth and gradient energy concentrating on a
  few high-|J| points are exactly the signatures to watch. They are what §1
  shows.
- **Negative control.** The same argument in dimension 4, with critical exponent
  4, is the basis of the negative control in Phase 0.

### 2.2 What a numerical approximation of a CS must show

Suppose J_n are ACS with sup|N_{J_n}| → 0 and ‖J_n‖_{C^{1,α}} ≤ C (equivalently,
uniform bounds on sup|J|, sup|∇J| and a Hölder or second-derivative bound).
Arzelà–Ascoli gives a C^{1,β} limit J_∞ with J_∞² = −1 and N = 0.
Low-regularity Newlander–Nirenberg (Nijenhuis–Woolf; Malgrange; Hill–Taylor
for rougher J) then makes J_∞ a complex structure on M.

So the target is a **resolution study**: as model capacity, collocation and
precision increase, sup|N| decreases at a consistent rate while
(sup|J|, sup|∇J|, sup|∇²J|) stay bounded. Report it as a frontier:

    F(Λ) = min { sup|N_J| : sup|J| ≤ Λ₀,  sup|∇J| ≤ Λ₁,  sup|∇²J| ≤ Λ₂ }

computed at increasing resolution.

| What the frontier does | What it means |
|---|---|
| F(Λ*) → 0 for some finite Λ* | Evidence for existence |
| F(Λ) → c(Λ) > 0 for every Λ tested, possibly while the L² loss → 0 | Evidence against existence *in that sector and representation* |

Neither outcome is a proof. Certification needs Phase 3.4 or a construction.

### 2.3 Topological sectors

- **Homotopy classes.** For a fixed orientation, ACS homotopy classes are
  [M, GL(6)/GL(3,ℂ)] = [M, ℂP³] = H²(M; ℤ) = ℤ. The map ℂP³ → ℂP^∞ is
  7-connected, and dim M = 6. The class k has c₁ = 2k·a, where a is the torus
  area class.
- **No characteristic-class obstruction.** In every sector, p₁ = 0 and a² = 0
  give c₂ = 0, and χ(M) = 0 gives c₃ = 0. So χ(𝒪) = c₁c₂/24 = 0 and
  Riemann–Roch does not single out a sector.
- **Conjugation models preserve k.** Every conjugation model preserves k:
  exponentials, Cayley transforms, any number of stages. The current lineage is
  k = 0.
- **The T²-invariant sector forces k = 0.** In the T²-invariant Parton frame,
  J is constant along each T²-fibre, so c₁ pairs to zero on that fibre.
- **One orientation suffices.** −J has the opposite orientation and the same
  Nijenhuis tensor.

### 2.4 Structural constraints on any CS

- **Not symplectic, not Kähler.** Every closed 2-form has [ω]³ = 0, since
  a² = 0 (already noted in the repo).
- **No nowhere-vanishing holomorphic 3-form.** Such an Ω would be closed. Since
  H³(M) = 0 it would be exact, so ∫Ω∧Ω̄ = 0, a contradiction.
- **K_M is not torsion.** Finite covers of M are again S⁴×T², so the same
  argument applies to the cover that would trivialise K^m.
- **Ruled-out ansätze.** Calabi–Yau-type ansätze are impossible, including
  SU(3) structures with dΩ = 0 and Hitchin's closed-stable-3-form flows.
- **Post-hoc check.** For any candidate, dΩ = θ∧Ω, with θ^{0,1} not ∂̄-exact,
  can serve as a check (Phase 3.3).
- **Orthogonality is not required.** Nothing forces orthogonality with respect
  to the round metric, so non-orthogonal J must stay allowed. But distortion has
  to be *bounded* along the convergent sequence of §2.2, not merely allowed to
  grow.

### 2.5 Symmetry, complexification, and a symmetry ladder

**Complexification.** If a compact group G acts holomorphically on a compact
complex manifold, the action extends to G^ℂ, because holomorphic vector fields
are complete. That makes symmetric sectors strongly constrained, and sometimes
classifiable. The Parton reference `J_ref = Eᵀ diag(J₂,J₂,J₂) E` is invariant
under U(2) × T²:
- T² rotates the two circles.
- U(2) acts on (x₁…x₄) = ℂ² and commutes with J₂ ⊕ J₂.

So the reference lies in every rung of the ladder below.

| Rung | Symmetry | Reduced problem | Notes |
|---|---|---|---|
| A | U(2)×T² (A′: SU(2)×T²) | ODE in x₅ | Cohomogeneity one; the open orbit of G^ℂ makes it classifiable. Solve exactly (ODE plus pole smoothness conditions) |
| B | T²_{S⁴}×T² = T⁴ (a maximal torus action) | 2D PDE | Probably excluded by Ishida's classification of maximal torus actions: the orbit space S⁴/T² is a bigon, which is not dual to a fan. **Verify.** If confirmed, a built-in negative control on M itself |
| C | S¹_rot × T² | 3D PDE | Rotation of S⁴ with fixed S² or with fixed poles |
| D | T² (circle factors) | **4D PDE on S⁴, exact**: network input x ∈ ℝ⁵ only. The Parton frame is T²-invariant | Forces k = 0. Holomorphic ℂ²-action; its generic orbits are complex surfaces |
| E | S¹_rot | 5D PDE | The fixed set must be complex: S²×T² (one rotation plane) or two elliptic curves (diagonal rotation). The S⁶ structure has Aut⁰ = ℂ* with fixed curves (Engel, Remark 8.3), so this type of sector is natural |
| F | none | 6D | |

**The ladder gives upper bounds.** The minimal residual is monotone along
inclusions of sectors (more symmetry gives a larger or equal minimum). Every
rung therefore bounds the next.

**Symmetry-breaking test.** At a G-invariant stationary point the residual is
G-invariant. Perturbations in a non-trivial isotypic component (e.g. torus
Fourier mode m ≠ 0) cannot lower |F|² at first order, because they are
orthogonal to F. The second-order form

    Q_m(δ) = |L_m δ|² + ⟨F, D²F(δ, δ̄)⟩

splits by mode m, and each piece is a 4D-sized problem. If Q_m has a negative
eigenvalue, symmetry breaking lowers the residual: follow that branch. If
Q_m ≥ 0 for all m below a cutoff, the symmetric point is a genuine 6D local
minimum to second order.

### 2.6 Gauge and ellipticity

**The gauge orbit.** Diff(M) preserves integrability, so integrable structures
come in infinite-dimensional orbits, and the linearised operator has the
infinite-dimensional kernel {L_X J}. The equation is first order and
overdetermined, but elliptic modulo gauge (a Dolbeault-type complex).

**Consequences.**
- Gauss–Newton normal matrices are singular in gauge directions.
- Newton needs a slice: Kuranishi gauge ∂̄*μ = 0, or damping.
- "Small residual ⇒ near a solution" holds only if the gauge-fixed operator has
  a bounded right inverse. That is measurable (Phase 3.2).

**Extra degeneracy in the current parametrisation.** On top of gauge, the
present exp(B) parametrisation has 18 pointwise null directions per stage: the
part of B that commutes with J.

### 2.7 What the 2026 constructions tell us

**The constructions.**
- S⁶ (Alpöge–Claude; Engel's independently verified account).
- S²×S⁴ (Guo–Fang–Lu; Liang–Liu–Zhang, conditional on the S⁶ result).

Both are complex 2-torus fibrations over ℙ¹ built from the (3,4,∞) family, with
multiple fibres (log transforms) and Mumford degenerations. They have algebraic
dimension 1, and the S⁶ structure has Aut⁰ = ℂ*.

**Implications.**
- The problem class is not hopeless.
- S⁶ and S²×S⁴ are now **difficulty-matched benchmarks with known answers**.
- Expected features of an S⁴×T² structure, if one exists: a holomorphic
  fibration, torus fibres, a circle action with fixed points, and strong
  non-orthogonality relative to round metrics. The derivatives may be large but
  should stay bounded.
- Engel (Remark 8.5) notes the machine is general and it is unclear which
  6-manifolds it produces. S⁴×T² is a natural test case.

**A cheap Leray filter already prunes naive variants.** Run the same machine over
an elliptic-curve base. The covector ψ dual to the ℤ-quotient period is always
monodromy-invariant, and a multiple of it extends. Then the d₂ argument forces
one of two outcomes: b₁ ≥ 3, or the pullback of the base class dies in H². For
S⁴×T² the map to the base would be homotopic to the projection, so neither
outcome is allowed. Any construction must therefore be subtler. (This is my
reasoning, not a published statement; check it.)

---

## 3. Formulation changes

### 3.1 Base points

- **k = 0.** Use the analytic Parton `J_ref`: it is exact, U(2)×T²-invariant,
  needs no learned derivatives, and has an explicit global (1,0)-frame
  Z_j = E_{2j−1} − iE_{2j}.
- **k ≠ 0.** Use J_k = Eᵀ(L_{f_k(u,v)} ⊕ J₂)E, where f_k: T² → S² has degree k
  and L_q is left multiplication by a unit imaginary quaternion on the
  (E₁…E₄) block. Verify c₁ = 2k numerically. (This is Deferred Phase C of
  `s4xt2_research_plan.md`.)
- **The learned s307 base** is in the same class as `J_ref`, so it is
  mathematically interchangeable. Keep it as a cross-check. The rule "no
  analytic input" matters for claims about *learning an ACS*; for CS existence
  only the final J matters.

### 3.2 Cayley–Beltrami parametrisation with re-centring

    J = (I + R) J_r (I + R)⁻¹ = J_r (I − R)(I + R)⁻¹,    R = ½ (M + J_r M J_r)

- **Equivalence with the Beltrami form.** R anticommutes with J_r; it is exactly
  the projection already used in `Correction`. Equivalently, R = μ + μ̄ with a
  Beltrami differential μ ∈ Hom(T^{0,1}_r, T^{1,0}_r). That is exactly 18 real
  degrees of freedom; for k = 0, μ is a 3×3 complex matrix field in the frame Z_j.
- **Bounded output.** Squash the output: R = (1 − δ) R̂ / √(1 + ‖R̂‖²_F). Then
  ‖R‖_op < 1 − δ, I + R is invertible, and the distortion of J relative to J_r
  is at most (2 − δ)/δ. This is built-in C⁰ control, with no runaway |J|.
- **Cost.** No matrix exponentials: one 6×6 solve per point.
- **Re-centring.** After each accepted block, freeze J_r ← J (cache its first
  jet, as `cs_jet_cache` already plans) and start a fresh zero-initialised head.
  Large total deformations are reached through bounded steps.
- **Residual.** In μ-coordinates, integrability is the cubic Maurer–Cartan
  equation ∂̄_r μ + (quadratic in μ) + (terms in N_r) = 0, valued in Λ^{0,2} ⊗
  T^{1,0} (9 complex components). Use it as the linearisation target in
  Gauss–Newton. Keep |N_J|_g (round metric) for reporting.

### 3.3 Objective

    L = E|N|² + β (E|N|ᵖ)^{2/p}  (+ γ · Reg)  (+ λ_g |∂̄_r* μ|²  in Newton phases)

- **Continuation in p.** Move p through 2 → 8 → 16, and hand off to an
  adversarial sup or augmented-Lagrangian witness constraints.
- **Regulariser.** Reg is a soft barrier on |∇J| > Λ₁ and |J| > Λ₀, plus an
  estimate of |∇²J| along random directions.
- **Frontier, not a single run.** Sweep Λ to trace F(Λ) (§2.2). Anneal γ only
  along the frontier, never to zero inside a single run.
- **Importance sampling.** Keep it, with unbiased 1/(Mq) weights; the repo
  measured a 1.85× variance reduction from J-norm mixtures.

### 3.4 Diagnostics (every checkpoint, fresh points)

- **Pointwise distributions.** Distributions of |N|, |J|_g, |∇J|_g and
  ρ = |N| / (|J|²_op |∇J|). ρ is a scale-free integrability fraction; a cleaner
  variant uses the J-Hermitian metric (g + JᵀgJ)/2 and its Levi-Civita
  connection (Gray–Hervella W₁⊕W₂ fraction).
- **Concentration measures.** The participation ratio (E|N|²)² / E|N|⁴, and the
  fraction of ∫|N|² carried by the top 1% of volume. A falling ratio means
  concentration.
- **Frontier and resolution.** The frontier points, and rates versus resolution.
- **Topology and symmetry.** c₁ on torus slices. The symmetry defect |L_ξJ| for
  the Killing fields ξ of each rung.
- **Spatial picture.** Structured-slice maps of where |N| and |∇J| live.

### 3.5 Optimisation

- **Exploration.** Adam with effective batches of 10⁴–10⁵, using QMC (scrambled
  Sobol through the normal inverse CDF on S⁴) plus importance sampling.
- **Precision.** Function-space Gauss–Newton or energy natural gradient
  (Müller–Zeinhofer), or Levenberg–Marquardt.
  - Use ≥ 5× more collocation points than parameters. Keep them fixed inside a
    solve and refresh between solves.
  - Accept a step only on held-out improvement, including the held-out maximum.
- **Preconditioning.** Precondition μ-updates with the inverse round Laplacian,
  which is diagonal in spherical harmonics × Fourier modes. It is spectrally
  equivalent to the gauge-fixed ∂̄-Laplacian of a smooth reference, so the
  conditioning becomes resolution-independent.
- **Function-class curriculum.** Raise the harmonic degree of μ progressively;
  this is the CS analogue of the ACS harmonic ladder.

### 3.6 Compute

- **GPUs.** Use A100/H100-class GPUs for float64, or at minimum full multi-core
  nodes. Run vmap and jvp over 10⁵ points, with float32 for exploration and
  float64 for polishing.
- **Expected gain.** Expect 10²–10³× throughput over the present single thread.
  The single-thread frame evaluator already handles 128 points in 0.15–0.28 s.

---

## 4. Phases and gates

### Phase 0 — calibration (prerequisite for interpreting anything)

**0.1 Infrastructure.**
- GPU port of the frame evaluator (`scratch/frame_s607.py` exists).
- At the s607 checkpoint, estimate the population gradient norm with batches of
  1k / 8k / 64k to confirm the stationary point.
- Measure how T²-invariant the s607 field already is: |∂_θ J| versus |∇_x J|,
  and its torus Fourier spectrum.

**0.2 New code.**
- Cayley–Beltrami head on the Parton base.
- The rung-D flag: network input x ∈ ℝ⁵ only.
- The supercritical objective and the diagnostics module.

**0.3 Positive controls**, started from generic far-from-integrable fields in the
right class:
- 4D: S³×S¹ (Hopf).
- 6D: S⁵×S¹ and S³×S³.

Two notes:
- **The frame trivialises S⁵×S¹.** In the Parton-type frame on S⁵×S¹, *every
  constant J is integrable* (the brackets are [u,v] = (x·u)v − (x·v)u, and N = 0
  identically). On S⁴×T², *no* constant J is integrable: the extra bracket
  projector forces a real eigenvector. So S⁵×S¹ only tests something
  nontrivial with non-constant planted targets (diffeomorphism pullbacks, or a
  resonant non-linear Hopf contraction) and with far starts.
- **Reduced controls.** For the reduced code, S³×S³ (Calabi–Eckmann, T² the
  maximal torus, quotient S²×S²) and S⁵×S¹ (quotient ℂP²) carry known
  T²-invariant CS. They need principal-bundle frames, so treat them as optional.

**0.4 Difficulty-matched benchmarks.** Run the identical pipeline on S²×S⁴
(c₁ = 0 class, already measured) and on S⁶ (unique class). If it plateaus there
in the same way, a plateau on S⁴×T² carries no information. If it makes real
progress, that is notable in itself: it would be the first numerical approximation
of the new structures.

**0.5 Negative control.**
- **Set-up.** T⁴ with J = L_{f(x₁,x₂)} for a degree-1 map f: T² → S² ⊂ Im ℍ.
  This class has c₁ ≠ 0. Every complex structure on T⁴ is a complex torus with
  c₁ = 0, so the class contains no CS.
- **Prediction.** The L² loss heads to 0 by concentration (p < 4), with
  |∇J| → ∞, while the supercritical frontier stays bounded below.
- **On M itself.** If Ishida's classification confirms rung B is empty, run it
  as a negative control on S⁴×T².

**Gate 0.**
- Positive controls reach sup|N| ≲ 1e-6, a ≥10⁴× reduction from the start, with
  bounded C¹ norms and a converging resolution study.
- The negative control shows the predicted dichotomy.

Until the gate passes, no S⁴×T² number is interpretable.

### Phase 1 — symmetry ladder on S⁴×T²

**1.1 Rungs A / A′ (ODE).**
- Derive the reduced system symbolically. In the invariant frame
  (∂_t, left-invariant fields, ∂θ₁, ∂θ₂), J is a 6×6 function of t alone, and
  N = 0 becomes first-order ODEs plus algebraic constraints.
- Apply the pole smoothness conditions (isotropy commutants: ±i for U(2);
  quaternionic for SU(2)).
- Find every solution, or prove there is none in the sector; compare with the
  classification of almost-homogeneous compact complex manifolds.
- Deliverable: a theorem-level statement for the sector, either way.

**1.2 Rung B (2D).**
- Check Ishida's classification first.
- Numerics: the rung-D code plus a soft T²_{S⁴}-invariance penalty, or a 2D
  spectral Newton solver.
- Expected outcome: a negative control.

**1.3 Rung D (T²-invariant, 4D): the main explorer.**
- 10⁵–10⁶ QMC points on S⁴; ≥ 32 seeds over initialisation amplitude and
  harmonic degree.
- Then spectral polishing. Represent 18 real fields as polynomials of degree
  ≤ L in x ∈ ℝ⁵ restricted to S⁴: about 157k unknowns at L = 16 and 684k at
  L = 24.
- Solve with Newton–Krylov, preconditioned by 1 + ℓ(ℓ+3).
- Run the resolution study (L, quadrature).

**1.4 Rungs C / E (S¹_rot).** 3D and 5D, both fixed-set types.

**Gate 1.**
- If some rung's frontier goes to 0: go to Phase 3, and to Phase 2 for the
  full-space check.
- If every rung shows a floor: run the mode-by-mode symmetry-breaking test
  (§2.5) at the best minimisers; negative modes seed Phase 2.

### Phase 2 — full 6D

**Route (i): symmetry-breaking continuation.** Start from the Phase 1 minimisers
along their negative second-order modes.

**Route (ii): fresh large-batch searches** in k ∈ {0, ±1, ±2} from the §3.1
bases, using:
- the Cayley–Beltrami head with re-centring;
- the p-continuation and the frontier sweep;
- Gauss–Newton polishing.

**Comparison.** Compare frontiers by sector and route. Promote only on frontier
improvement on common fresh data, never on mean-square loss alone. Keep the
receipt and hash discipline; drop the long CPU chains.

### Phase 3 — precision, verification, identification

1. **Spectral polishing.** Distil the best candidate into SH(S⁴) × Fourier(T²)
   (or the reduced basis), run Newton–Krylov to 1e-10, and check for
   exponential convergence in the degree.
2. **Linearised operator.** The smallest singular values of the gauge-fixed
   operator should be stable under refinement. That is quantitative stability,
   and it is the input a Kantorovich argument needs.
3. **Independent checks.**
   - Local holomorphic coordinates: solve ∂̄_J z = 0 on ≥ 2⁶ patches, check full
     Jacobian rank, and check that transition maps are holomorphic.
   - c₁ sector.
   - K_M not torsion: θ^{0,1} = ∂̄Ω/Ω is not ∂̄-exact (§2.4).
   - Holomorphic vector fields: the kernel of ∂̄ on T^{1,0}.
   - Signs of a fibration or of algebraic dimension.
4. **Rigorous certification.** Do the reduced sector first: an
   interval-arithmetic Newton–Kantorovich proof on the spectral approximant.
   For a genuine theorem, either this or a construction is needed.

### Track T — theory (parallel, cheap, high leverage)

- **Sector classification.** Rung A via almost-homogeneous manifolds; rung B via
  Ishida; rung E via the complex fixed-point set.
- **Construction search.** Work with algebraic geometers on Engel's machine
  (Remark 8.5) and its variants. A small integer-linear-algebra tool would help:
  monodromy and special-fibre data in, π₁ / H* / χ / w₂ / p₁ out (the recipe of
  Engel §§5–7). Use it to filter candidates against S⁴×T²: π₁ = ℤ², b = (1,2,1,0,1,2,1),
  a² = 0. Identifying the diffeomorphism type with π₁ = ℤ² needs surgery theory.
- **Feedback from numerics.** If a rung finds a structure, extract its geometry
  (the ℂ²-orbit foliation in rung D, fixed sets in rung E) to guide a proof.

---

## 5. Decision table

| Observation | Interpretation | Action |
|---|---|---|
| Gate 0 fails on positive controls | Method limitation | Fix the formulation and optimiser before any S⁴×T² compute |
| Controls pass; S⁶ and S²×S⁴ also plateau | Hard-structure regime not reachable yet | Lean on Phase 1 reductions and Track T; S⁴×T² plateaus are uninformative |
| L² falls but the frontier floor persists; participation ratio falls, \|∇J\| grows | Concentration (§2.1) | Not progress. Tighten Λ and inspect the defect location |
| Train ↓, held-out ↑ | Under-resolved quadrature | More collocation, a spectral basis, lower degree; not more width |
| A rung's frontier → 0 under refinement | Strong evidence in that sector | Phase 3; symmetry-breaking check; Track T identification |
| Symmetric floor with a negative Q_m | Symmetry is the obstruction | Phase 2 route (i) |
| Floors in all sectors and k, stable across resolution and seeds | Bounded negative result | Report as such; never claim non-existence |

---

## 6. First two weeks (concrete)

1. GPU or multi-core port of the frame evaluator. Large-batch gradient-norm
   measurement at the s607 checkpoint. T²-invariance and Fourier spectrum of the
   s607 field.
2. Cayley–Beltrami head on the Parton base in 6D, plus the rung-D (x-only)
   variant in the same code path.
3. Diagnostics module (§3.4) and the frontier driver.
4. T⁴ negative control and the S³×S¹ / S⁵×S¹ / S³×S³ positive controls under
   the new formulation.
5. Symbolic derivation of the rung A (U(2)×T²) ODE and its pole conditions.
6. Literature check of Ishida (rung B) and of almost-homogeneous classification
   (rung A).

## 7. Stop doing

- Continuing the s407 / s507 / s607 lineage with small-batch Adam at LR ≤ 1e-5.
- Accepting on mean or RMS Nijenhuis without regularity bounds.
- Requiring learned ACS bases for the CS search.
- Single-thread CPU job chains for 6D second-derivative work.

## 8. Honest assessment

- **Existence is open.** The plan is built so each phase yields a usable result
  even with no CS found: sector theorems, a calibrated method, and a bounded
  negative result. A numerical approximation of the new S⁶ or S²×S⁴ structures
  would be notable in its own right.
- **Numerics give evidence, not proof.** A genuine theorem needs either a
  Kantorovich-type certificate (realistic first in a reduced sector) or a
  construction.
- **Main risks.**
  - A CS, if one exists, may be "wild" relative to round metrics, as the
    multiple-fibre and degeneration structure of the S⁶ example suggests. That
    could push the needed resolution out of reach in 6D.
  - It may exist only without symmetry, or in sectors not yet explored.
  - 6D quadrature cost.

---

## References

- P. Engel, *Complex structures on S⁶* (2026): https://philip-engel.github.io/S6.pdf
  (independently verified account; Remarks 8.3 and 8.5 are used above).
- Claude and L. Alpöge, *The (3,4,∞) modular family of 2-tori … is a complex
  structure on S⁶* (2026): https://alpo.ge/s6.pdf
- S. Guo, E. X. Fang, J. Lu, *A complex structure on S²×S⁴*:
  https://arxiv.org/abs/2609.05504
- R. Liang, C. Liu, Y. Zhang, *Remarks on the complex structures on ℙ³ and
  S²×S⁴* (conditional): https://arxiv.org/abs/2608.29651
- P. Landweber, *Complex structures on open manifolds*, Topology 13 (1974):
  https://www.sciencedirect.com/science/article/pii/0040938374900391
- A. Phillips, *Submersions of open manifolds*, Topology (1967); M. Gromov,
  *Partial Differential Relations* (1986).
- A. Nijenhuis, W. Woolf, *Some integration problems in almost-complex and
  complex manifolds*, Ann. Math. (1963); C. D. Hill, M. Taylor, *Integrability of
  rough almost complex structures*, J. Geom. Anal. (2003).
- H. Ishida, *Complex manifolds with maximal torus actions*, J. reine angew.
  Math. 751 (2019).
- M. Parton, *Explicit parallelizations on products of spheres*:
  https://arxiv.org/abs/math/0009156
- J. Müller, M. Zeinhofer, *Achieving high accuracy with PINNs via energy natural
  gradient descent*, ICML 2023.
- Rathore et al., ICML 2024: https://proceedings.mlr.press/v235/rathore24a.html
- Euh–Sekigawa (2013): https://arxiv.org/abs/1301.6835
