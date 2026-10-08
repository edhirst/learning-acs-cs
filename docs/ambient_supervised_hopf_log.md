# Supervised reproduction of the Hopf structures in the ambient (overlap-free) framework

Running log (branch `no_overlaps`). Question: can the new chart-free framework
(`docs/ambient_pipeline_math.md`) represent the known complex structures on
S³×S¹ and S⁵×S¹ accurately enough — in the first-order (value + derivative)
sense that the Nijenhuis tensor requires — where the old chart/overlap
framework on `main` could not (supervised A-param floor ‖N‖ ≈ 0.73–0.8 on
S³×S¹ at J-fit ~1e-4; see the main-branch notes)?

All numbers below are fp64. "1e-10" in any residual is the √(·+1e-20) floor
of the code, i.e. zero to machine precision.

Code version: all runs used `extra/campaign_unit.py` @ 558d158. Upstream 7b8eb49 (Ed Hirst, s4×s1×s1 support;
upstream tip 1383f80 only touches .gitignore/scratch) changed `chart_protocol`'s sampler (`sphere_uniform_sample` instead of `SphereProduct.sample_points`) and its
overlap routine, so re-evaluating the "paper chart protocol" column today gives statistically, not bitwise,
identical numbers. Main's exact protocol (`main_band_points`, implemented in our script) is unaffected.
Re-checked after the fast-forward (2026-10-01): all 6 bases load and certify identically, `--verify` is unchanged,
and the stored S³×S¹ / S⁵×S¹ checkpoints reproduce their recorded rel. errors and ‖N‖ within sampling noise.

---

## Stage 1 — the target: Hopf structure as an ambient field  (DONE, 2026-09-29)

`src/geometry/hopf_ambient.py`. On M = S^{2n−1}×S¹ ⊂ ℝ^{2n}×ℝ², p = (q, w),
∂_θ = R(w) = (−w₁, w₀), m := I q:

    J_qq = I − m qᵀ + q mᵀ,   J_qw = c·m Rᵀ,   J_wq = −(1/c)·R mᵀ,   J_ww = 0

(J = Λ⁻¹ I Λ with Λ(v, s∂_θ) = v + c s q; Hopf manifold (ℂⁿ∖0)/(z ∼ e^{2πc} z)).
c = 1/(2π) reproduces main's `hopf_complex_structure(_fixed)` (z ∼ e·z);
I = Lᵀ J₀ L transports main's J₀ through the coordinate permutation between
`hopf_map_F` (pole coordinate first) and `ball_to_sphere` (pole coordinate last).
Every block has degree ≤ 2 in (q, w).

Verification (reproducible with `python extra/ambient_supervised_hopf.py --verify`; 3000 ambient pts +
400 chart pts per patch):

| check | S³×S¹ (n=2) | S⁵×S¹ (n=3) |
|---|---|---|
| tangency ‖P J P − J‖ (unprojected formula) | 1.8e-15 | 1.8e-15 |
| ‖J² + P‖ max | 1e-10 (floor) | 1e-10 (floor) |
| intrinsic ‖N‖ max | 1e-10 (floor) | 1e-10 (floor) |
| pull-back to the 4 ball charts vs `hopf_complex_structure_fixed`, max rel. error | **5.2e-14** | **1.6e-14** |
| ‖J‖_F (mean) | 6.44 | 6.60 |
| intrinsic Dirichlet energy E‖∇J‖² | 83.0 | 166.0 |

The target is strongly anisotropic w.r.t. the round metric (the 1/c = 2π block):
this is intrinsic to main's choice λ = e, not a chart artefact.

## Design decisions (and why)

* **Parametrization = the framework's CS stage, unchanged**: J = e^{B₂}e^{B₁} J_base e^{−B₁}e^{−B₂},
  B_s = P·mlp_s(p)·P general (not antisymmetric), GELU MLP, widths 128×3 (d=4) / 256×4 (d=6),
  last layer ×0.5 (O(1) init). J² = −P exact by construction (the ambient analogue of main's A-param).
* **Base must be homotopic to the target** (the conjugation family is a homotopy; det e^B > 0 fixes
  the orientation). S⁵×S¹: ACS classes = sections of a trivial CP³-bundle; H² = 0 and
  π_k(CP³) = π_k(S⁷) = 0 for 3 ≤ k ≤ 6 ⇒ orientation is the only invariant, and −J flips it (n = 3).
  S³×S¹: fibre S², H³(M; π₃S²) = ℤ ⇒ infinitely many classes; −J does NOT flip orientation (n = 2).
* **Primary base = the unit's own ladder-certified ACS** (target-agnostic, exactly as in the campaign),
  orientation-checked against the target (S⁵×S¹: use −J* if needed). Successful supervised convergence
  (small sup error) itself certifies homotopy. Control (only if needed on S³×S¹): a base homotopic to
  the target by construction.
* **Loss** (supervised, Sobolev, λ_D = 1 as main's Phase A):
  E_p[ ‖J − J_t‖²_F + λ_D Σ_a ‖D_{Pe_a}(J − J_t)‖²_F ] — tangential directional derivatives only
  (extension-independent; Σ over the Parseval frame Pe_a = intrinsic).
* **Metrics**: relative value/derivative errors (RMS and sup), intrinsic ‖N‖ (mean, sup, adversarial sup),
  the paper's chart protocol (per-patch hard-cutoff flat ‖N‖), and chart-coordinate J / ∂J errors on the
  same band as main — for an apples-to-apples comparison with the old Phase A.

Script: `extra/ambient_supervised_hopf.py` (`--verify` reproduces Stage 1, including the chart
derivatives: rel ∂J vs main ≤ 2.5e-13 on all patches). On the hard-cutoff band, main's chart-coordinate
‖∂J‖ has RMS 25.5 (S³×S¹) / 28 (S⁵×S¹), vs the intrinsic √E_D = 9.1 / 12.9: the charts inflate the
derivatives of the same structure by ≈ 2.2–2.8×.

## Stage 2 — bases  (DONE)

Production ACS units (campaign settings: ladder 3000 steps/class, 10k MLP steps, L-BFGS 200),
`outputs/ambient_supervised/bases/`:

| unit | accepted class | adversarial ‖K²+P‖ sup (final) | ‖N(J*)‖ mean | orientation vs J_t |
|---|---|---|---|---|
| s3xs1 seed 7 | quad_h | 1.1e-3 | 1.03 | same |
| s3xs1 seed 8 | quad_h | 9.1e-4 | 0.96 | same |
| s3xs1 seed 9 | quad_h | 1.5e-3 | 0.94 | same |
| s3xs1 seed 10 | all2 (last class; all lower classes rejected) | 2.2e-3 | 1.02 | (not used) |
| s5xs1 seed 7 | bilin (ladder sup 0.262) | 3.1e-4 | 17.06 (sup 17.22: near-homogeneous) | same |
| s5xs1 seed 8 | bilin (ladder sup 0.218) | 6.3e-4 | 17.05 (sup 17.15) | **opposite** → −J* used (auto, n odd) |

Findings:
* **The ladder is initialisation-dependent on S³×S¹.** With the full ladder, `lin` was rejected on all 4
  seeds, plain and Dirichlet retry (adversarial sup ≈ 2.00 = a degenerate wall), whereas a direct
  `--feat lin` run from a different random init certified (sup 0.167, ‖N‖ ≈ 2.74; an exploratory local
  smoke run with 1500 MLP steps, not kept).
  The doc's "S³×S¹ accepted in lin" is therefore not robust; production seeds 7–9 certify in quad_h,
  seed 10 only in the last class (all2). On S⁵×S¹ both seeds certify in bilin, as the doc states, and both
  land on essentially the same ‖N(J*)‖ ≈ 17.05 (suggesting a canonical bilinear-class ACS).
* **Orientation is a real, silent obstruction.** The `lin`-class smoke base (seed 7) induces the OPPOSITE
  orientation to J_t: from it the conjugation family (det e^B > 0) can never reach main's Hopf structure —
  and the campaign's CS stage from that base could only find a structure of the other orientation.
  Remedy implemented (`--base-sign auto`): n odd → −J*; n even → transport J* by the orientation-reversing
  isometry θ ↦ −θ, S(q, w₀, w₁) = (q, −w₀, w₁), i.e. base = S·J*(S p)·S (still certified; all residuals
  S-invariant). The S³×S¹ production bases 7–9 need no fix; the S⁵×S¹ seed-8 base had the opposite
  orientation and its run used −J* automatically (recorded as `orientation_fix: negate` in its result.json).
* Orientation test: sign det[v₁, Jv₁, …, v_n, Jv_n] in the global orthonormal Hopf frame
  F_k = (e_k − q_k q, q_k ∂_θ); J_t has sign −1 in this frame for both n = 2, 3.

## Stage 3 RESULT — S³×S¹  (DONE, 2026-09-29)

**The ambient framework reproduces main's Hopf structure on S³×S¹ ~10³× better in value, ~10³–10⁴× better
in derivative and ~30–40× better in ‖N‖ than main's best supervised run, with a 55× smaller network.**

Main's exact held-out protocol (overlap band, seed 99999, mean over the 4 patches; `eval.json`):

| S³×S¹ | main, old framework: A-param Phase A + L-BFGS (4 × 0.6M params, Fourier) | ambient, base 7 | ambient, base 8 |
|---|---|---|---|
| J-fit E‖δJ‖²_F | 1.75e-4 | **2.30e-7** | **1.04e-7** |
| ∂J-fit E‖δ∂J‖²_F | 0.147 | **7.46e-5** | **3.26e-5** |
| rel. J / rel. ∂J error (RMS) | — | 8.2e-5 / 2.8e-4 | 5.5e-5 / 1.8e-4 |
| flat-chart ‖N‖ mean | **0.73** | **0.0261** | **0.0171** |
| parameters | ≈ 2.4M | 43k | 43k |

Intrinsic (ambient, round-metric) quantities, 4000 fresh points + adversarial ascent:

| | base 7 | base 8 |
|---|---|---|
| rel. value error RMS / sup / adversarial sup | 7.2e-5 / 1.6e-4 / 1.9e-4 | 5.0e-5 / 1.3e-4 / 1.4e-4 |
| rel. tangential-derivative error RMS / sup | 2.2e-4 / 4.3e-4 | 1.6e-4 / 3.3e-4 |
| intrinsic ‖N‖ mean / sup / adversarial sup | 0.0101 / 0.031 / 0.037 | 0.0069 / 0.023 / 0.025 |
| paper chart protocol (hard cutoff) ‖N‖ | 0.026 | 0.018 |
| ‖J² + P‖ | floor (exact) | floor (exact) |
| wall time (3 CPU threads) | ≈ 2 h | ≈ 2 h |

Trajectory (base 7): ‖N‖ 0.52 (1k) → 0.20 (10k) → 0.073 (20k) → 0.016 (30k) → 0.010 (40k) — still
decreasing when the cosine schedule ends (lr → 1e-5): this is **not** a floor of the network class. The
L-BFGS finisher (4 × 100 its, 2048 pts) gained < 0.2% (two chunks reverted).

Consistency with the linearisation framework of main: ‖N‖ tracks the RMS derivative error,
‖N‖ ≈ 3.0 × √(∂J-fit) here (both bases) vs ≈ 1.9 × on main — N is dominated by 𝓜(δ∂J) in both frameworks;
what changed is that δ∂J is ~70× smaller.

### Ablation: value-only supervision (no derivative term), S³×S¹ base 7, 40k steps  (DONE)

| | Sobolev (base 7) | **value only (base 7)** | main (value+∂J Sobolev, A-param) |
|---|---|---|---|
| main protocol J-fit / ∂J-fit | 2.3e-7 / 7.5e-5 | **1.4e-6 / 4.3e-4** | 1.75e-4 / 0.147 |
| main protocol flat ‖N‖ | 0.026 | **0.064** | 0.73 |
| intrinsic rel val / rel der (RMS) | 7.2e-5 / 2.2e-4 | 1.8e-4 / 5.5e-4 | — |
| intrinsic ‖N‖ mean (adv sup) | 0.010 (0.037) | 0.025 (0.090) | — |
| rel der / rel val | 3.1 | **3.0** | ≈ 12.6–14.7 (Fourier gain, see main notes) |

**Derivative supervision helps (~2.5× in ‖N‖) but is not essential**: without it, the ambient GELU MLP
(no Fourier features) keeps its relative derivative error at ≈ 3× its relative value error — the same
ratio as with Sobolev supervision. Contrast main, where value-only fits were catastrophic (Tomas's sweep:
8% value error ↔ 51% derivative error, Nij ≈ 27) because the Fourier layer's value→derivative gain is
2π·s·√d ≈ 12.6. This is the mechanism behind the whole improvement.

### Matched comparison: main's anisotropic Hopf structure vs the round-orthogonal one  (DONE)

Same base (s3xs1 seed 8), same network, same 10k-step schedule, 100 L-BFGS its (`runs/s3xs1_sob_b8_10k_{cmain,c1}`):

| target | rel. value error | rel. derivative error | intrinsic ‖N‖ mean (adv sup) | paper chart protocol | flat ‖N‖ on main's band |
|---|---|---|---|---|---|
| c = 1/(2π) (main, λ = e; ‖J‖_F 6.44, E_D 83) | 3.35e-4 | 8.3e-4 | 0.038 (0.13) | 0.10 | 0.097 |
| c = 1 (round-orthogonal, λ = e^{2π}; ‖J‖_F 2, E_D 8) | 3.29e-4 | 8.7e-4 | **0.0046 (0.013)** | 0.015 | 0.015 |

(For c = 1 the main-protocol J-fit is against main's c = 1/(2π) structure and is not meaningful; its flat ‖N‖ is.)
**Identical relative accuracy, 8.3× smaller ‖N‖** for the orthogonal target: N is bilinear, N ~ |J|·|δ∇J| + |∇J|·|δJ|,
and the absolute errors scale with the target (‖J‖ 6.44 vs 2, ‖∇J‖ √83 vs √8 ⇒ factor ≈ 3.2 × 3.2 ≈ 10).
This quantitatively explains Phase B: at fixed network precision the round-orthogonal Hopf structure admits ~8×
smaller ‖N‖, so the Nijenhuis-only flow slides toward it. It also means main's choice λ = e made the supervised
problem ~8× harder in ‖N‖ than necessary. (10k vs 40k steps for c = 1/(2π): ‖N‖ 0.038 vs 0.007–0.010.)

## Stage 3 RESULT — S⁵×S¹  (PARTIAL: local runs stopped at the 2026-09-29 wrap-up — FINAL results: section "S⁵×S¹ — final results (2026-10-06)" at the end)

Network 256×4 GELU, 2 stages, 233k parameters; ≈ 0.55 s/step locally (8×8 matrices), hence not converged locally.

| S⁵×S¹ | short run, base 8: 8k steps, schedule fully annealed + 100 L-BFGS | 80k-step run, base 7: stopped at step 22k (lr ≈ 8e-4) |
|---|---|---|
| rel. value error RMS / sup / adversarial sup | 1.05e-3 / 2.0e-3 / 2.75e-3 | 7.4e-3 / 1.2e-2 / 1.3e-2 |
| rel. derivative error RMS | 2.6e-3 | 8.9e-3 |
| intrinsic ‖N‖ mean (adversarial sup) | 0.159 (0.43) | 0.387 (1.06) |
| main protocol J-fit / ∂J-fit / flat ‖N‖ | 4.8e-5 / 1.3e-2 / **0.39** | 1.9e-3 / 0.13 / 0.90 |
| paper chart protocol ‖N‖ | 0.39 | 0.88 |
| wall time | 75 min (4 threads) | 203 min (6 threads) |

(80k run evaluated from `state.pt` at step 22000 — `eval_partial_step22000.json`; that same checkpoint is left on
disk for the cluster to resume from. Mid-schedule, lr still high, so it is noisier than the annealed short run.)

References: old results table S⁵×S¹ CS 6.42; Tomas's chart-framework Sobolev "architecture ceiling" 1.14
(600-step fit); S³×S¹ at a comparable 10k-step schedule: intrinsic ‖N‖ 0.038. The short run is **16× below the
table and 3× below the chart ceiling after 75 min, and still improving** (‖N‖ 0.40 → 0.24 → 0.17 → 0.16 over
its last 3k steps). No sign of a floor; convergence requires the full 80k schedule — cluster.

## Stage 4 — Phase B: unsupervised (Nijenhuis-only) continuation from the supervised endpoint

S³×S¹, from `s3xs1_sob_b8` (rel val 5e-5, ‖N‖ 0.0070), loss E‖N‖² only, same network/base, Adam lr 1e-3
(main's Phase B value) cosine over 10k steps. (`runs/s3xs1_phaseB_b8_lr1e-3`)

| step | rel. distance to J_hopf (RMS / sup) | intrinsic ‖N‖ mean / sup |
|---:|---|---|
| 0 | 5.0e-5 / 1.3e-4 | 0.0070 / 0.018 |
| 250 | **0.645 / 0.655** | 0.0075 / 0.018 |
| 1000 | 0.752 / 0.757 | 0.0094 / 0.018 |
| 2000 | 0.799 / 0.804 | 0.016 / 0.035 |
| 5000 | 0.821 / 0.826 | 0.0018 / 0.0038 |
| 7250 | 0.824 / 0.83 | 0.0015 / 0.0038 |
| **final** (10k Adam + 4×100 L-BFGS) | 0.825 / 0.830 (adv 0.830) | **7.9e-4 / 2.2e-3 (adversarial sup 2.5e-3)** |

Final endpoint, main's protocol: flat ‖N‖ **2.5e-3** (main's best supervised: 0.73; main's Phase B: 0.43),
J-fit vs J_hopf 23.7 (main's Phase B: 20); paper chart protocol 2.6e-3. The most integrable field of the study.
Endpoint identification is unchanged from step 6000 (it has stopped moving): ‖J‖_F = 2.007 — since
‖J‖²_F ≥ d for any ACS with equality iff J is orthogonal, J is within 0.35% of orthogonal — E_D 8.07,
frame components constant to 6.0%, K̄ nearly an orthogonal complex structure. Close to the round-orthogonal Hopf
family but not exactly a constant member of it (plausibly a diffeomorphism image, equally integrable).

* **J leaves J_hopf immediately (within 250 steps) while staying integrable** (‖N‖ never above 0.02, and it
  ends 5× BELOW the supervised fit). The displacement is uniform (sup ≈ RMS): not a defect.
* Mechanism: the zero set of N is not isolated — N_{φ*J} = φ*N_J for every diffeomorphism φ isotopic to the
  identity, plus the moduli of Hopf structures — so J_hopf sits on an infinite-dimensional valley floor;
  Adam's per-coordinate normalisation turns the vanishing gradients along the valley into O(lr) steps.
* **Identification of the endpoint** (checkpoint step 6000; `extra/ambient_phaseB_identify.py RUN_DIR BASE_DIR`,
  which on the final endpoint reproduces the numbers quoted below for it): ‖J‖_F = 2.008
  (= √d: an ORTHOGONAL ACS), intrinsic Dirichlet energy 8.07 (target 83.0; round-orthogonal Hopf J_{c=1}: 8.00);
  closest member of the J_c family: c' = 0.96 (rel. distance 0.23). In the global Hopf frame
  F_k = (e_k − q_k q, q_k ∂_θ) its components are nearly CONSTANT (mean deviation 5.9%), with average K̄ an
  orthogonal complex structure (‖K̄²+I‖ 3.6e-3, ‖K̄K̄ᵀ − I‖ 1.0e-2) in the orientation class of I, a rotation of I.
  Every constant K in this frame descends to an integrable structure (a linear complex structure on ℂ²,
  invariant under z ↦ e^{2π} z). ⇒ **The unsupervised flow moves from main's anisotropic Hopf structure
  (λ = e) toward the round-orthogonal Hopf structures (λ = e^{2π}; an S²-family of K)** — the metrically
  simplest points of the moduli (smallest derivatives), which the network represents best.
* **Learning-rate control** (`runs/s3xs1_phaseB_b8_lr1e-4`, lr 1e-4 → 1e-6, 4000 steps, no L-BFGS):

  | step | 250 | 500 | 1000 | 2000 | 4000 (final) |
  |---|---|---|---|---|---|
  | rel. distance to J_hopf, lr 1e-4 | 0.053 | 0.080 | 0.115 | 0.147 | 0.158 |
  | rel. distance to J_hopf, lr 1e-3 | 0.645 | 0.706 | 0.752 | 0.799 | 0.818 |

  Final (lr 1e-4): intrinsic ‖N‖ 4.1e-3 (adv sup 1.7e-2), main-protocol flat ‖N‖ 1.0e-2, still decreasing.
  10× smaller lr ⇒ 12× smaller initial drift: the drift scales with lr, as predicted by Adam's per-coordinate
  normalisation acting on the (near-)flat directions of the N = 0 set.
* Contrast main: Phase B there drifted to J-fit 20 but stalled at ‖N‖ ≈ 0.43 (chart norm) — never integrable,
  because the overlap loss competed with N. Here the drift is harmless: it is motion along the space of
  complex structures.
* Consequence for the unsupervised campaign: "does the CS search recover J_hopf?" is the wrong question —
  N = 0 defines a gauge-degenerate set; the right diagnostics are integrability (‖N‖ sup) plus invariants of the
  endpoint (here: Hopf-type, λ = e^{2π}).

## Stage 3 — supervised Sobolev fit: setup

Runs `outputs/ambient_supervised/runs/`: s3xs1 from bases 7 and 8; 40k Adam steps (batch 192, lr 1e-3 → 1e-5
cosine, clip 1), then chunked L-BFGS (4×100 its, fresh 2048-pt batch per chunk, revert if the held-out loss
worsens). Net 128×3 GELU, 2 stages, 43k parameters (main's Phase A: 4 patches × [256,512,512,256] + 128
Fourier features ≈ 2.4M). Init (O(1) random conjugation of the base): rel. value error 5.8, ‖N‖ ≈ 1.5e3.
Pipeline smoke-tested (300 steps) from the reflected lin base; timing ≈ 0.17 s/step on 3 threads.
s5xs1 from base 7: 80k Adam steps, net 256×4, 2 stages; init rel. value error 5.2, ‖N‖ ≈ 8e2.

**Main's evaluation protocol replicated exactly** (`main_band_points`, from the stashed
`extra/diagnostics/precision_2026_06_16/precision_lbfgs_only.py::eval_held_out`): seeds 99999,
`SphereProduct(2n−1,1).sample_points(3000)`, keep |x| < 0.85 in both charts of each factor (overlap band;
2348 / 2379 points kept for n = 2 / 3, |x_α| ∈ [0.08, 0.85]), per patch J-fit = E‖δJ‖²_F, dJ-fit = E‖δ∂J‖²_F,
flat ‖N‖ mean; averaged over the 4 patches. Validated on the exact target: J-fit 9e-30, dJ-fit 2e-27, ‖N‖ at
floor. Reported as `main_protocol` in result.json / eval.json (`--eval-only RUN_DIR` for finished runs).

Early observation (step 1000, 2.7 min): rel. value error 1.7% / 1.4% (bases 7 / 8), **sup of the relative
value error only ≈ 3%** ⇒ no localised defect ⇒ the quad_h bases are in the target's homotopy class
(the π₃(S²)-type obstruction does not occur for them); rel. derivative error ≈ 2%; intrinsic ‖N‖ 0.52 / 0.36.

---

## Wrap-up (2026-09-29, ~20:35) — conclusions

1. **The overlap-free framework reproduces the known structure on S³×S¹; the old framework could not.**
   Under main's exact evaluation protocol: flat ‖N‖ 0.017–0.026 vs 0.73, J-fit ~10³× and ∂J-fit ~10³–10⁴×
   smaller, with a 43k-parameter network instead of 2.4M, in 2 CPU-hours. ‖N‖ was still decreasing at the end of
   the schedule, so this is not a floor of the class.
2. **Mechanism:** the ambient GELU MLP (no charts, no Fourier features) has an O(1) value→derivative error gain
   (rel der / rel val ≈ 3.0–3.1, with or without derivative supervision), vs 2π·s·√d ≈ 12.6 for main's Fourier
   charts; ‖N‖ ≈ 3·√(∂J-fit), dominated by 𝓜(δ∂J) as on main. Value-only supervision already gives ‖N‖ 0.064.
3. **S⁵×S¹ (partial):** after 75 min, fully annealed 8k steps: flat ‖N‖ 0.39 on main's protocol (table 6.42; chart
   Sobolev ceiling 1.14), still improving; the 80k schedule is needed → cluster.
   *(Final, 2026-10-06: flat ‖N‖ 0.064 ± 0.005 over 4 runs of the 80k schedule — see the last section.)*
4. **Unsupervised continuation (Phase B)** stays integrable (‖N‖ down to 7.9e-4, adversarial sup 2.5e-3) but leaves
   main's Hopf structure (λ = e) for a near-orthogonal Hopf-type structure close to λ = e^{2π}; drift ∝ lr (Adam on
   the gauge-degenerate N = 0 set). A matched comparison shows why: at equal relative accuracy the orthogonal
   structure has 8.3× smaller ‖N‖.
5. **Caveats found along the way:** the base's orientation must match the target's (silent obstruction;
   handled automatically); ladder acceptance on S³×S¹ is initialisation-dependent (lin rejected on 4/4 production
   seeds; quad_h accepted on 3/4, all2 on the 4th).

## Cluster continuation (prepared, not submitted)

`pbs/ambient_supervised.pbs` + `jobs_ambient_supervised.txt` + `jobs_ambient_supervised_phaseB.txt`:
phase 1 = s5xs1 base 7 (80k, **resumes from the local step-22000 checkpoint**), s5xs1 base 8 (80k, fresh),
s3xs1 base 8 (120k, "how low can ‖N‖ go"); phase 2 = Nijenhuis-only continuation (10k, lr 1e-3) from each.
Before `qsub`: rsync `outputs/ambient_supervised/bases/` and `outputs/ambient_supervised/runs/s5xs1_sob_b7/state.pt`
to the same paths in the cluster checkout (bases are small; missing bases are re-certified automatically).
Unit naming and syntax dry-run checked locally (GNU `xargs -a` on the cluster, as the existing campaign scripts).

## Files (committed to branch no_overlaps on 2026-10-01; outputs are gitignored and stay local)

* `src/geometry/hopf_ambient.py` — closed-form ambient Hopf structure (any n, any c).
* `extra/ambient_supervised_hopf.py` — experiment: `--verify`, supervised (`--loss sobolev|value`), Phase B
  (`--loss nij --init-from`), `--eval-only RUN_DIR` (also partial runs from state.pt), main's exact protocol.
* `extra/ambient_phaseB_identify.py` — endpoint identification (Hopf-family scan, Hopf-frame constancy).
* `pbs/ambient_supervised.pbs`, `jobs_ambient_supervised.txt`, `jobs_ambient_supervised_phaseB.txt` — cluster job.
* `docs/ambient_supervised_hopf_log.md` — this log.
* Outputs (gitignored): `outputs/ambient_supervised/bases/` (6 certified ACS units), `outputs/ambient_supervised/runs/`
  (every run: progress.jsonl, final.pt/state.pt, result.json and/or eval.json, .log).

---

## Status 2026-10-05

* Nothing has run since the 2026-09-29 wrap-up. No cluster results yet: the job above was prepared but not
  submitted.
* **Laptop timing for the paused S⁵×S¹ 80k run** (`runs/s5xs1_sob_b7`, `state.pt` at step 22000, lr 8.26e-4),
  measured on the idle laptop by resuming a scratch copy of the checkpoint for 300 steps:
  0.503 s/step with 6 threads and 0.507 s/step with 10 threads. Both figures include an evaluation every 100
  steps; the production run evaluates every 1000, so the true step time is ≈ 0.45–0.50 s. Threads saturate at
  about 6, because the batch is small.

  | remaining work | estimate |
  |---|---|
  | 58k Adam steps (22k → 80k) | 7.3–8.1 h |
  | L-BFGS finisher (4 × 100 its, 2048 points) | ≈ 0.7–1 h |
  | final evaluation (adversarial search, chart protocol, main's protocol) | ≈ 0.2–0.3 h |
  | **total** | **≈ 8.5–9.5 h** |

  A second S⁵×S¹ run in parallel (6 threads each) would cost only ≈ 10% (cf. 0.555 s/step under contention on
  2026-09-29). Run under `caffeinate -i` so the laptop does not sleep. The command resumes automatically from
  `state.pt`:
  `OMP_NUM_THREADS=6 caffeinate -i python extra/ambient_supervised_hopf.py --case s5xs1 --seed 7 --base-from outputs/ambient_supervised/bases/s5xs1_acs_s7 --out outputs/ambient_supervised/runs/s5xs1_sob_b7 --steps 80000`

  *Correction (later on 2026-10-05, from measurements; table above left as written):* the last two rows were
  guesses and the final-evaluation row was too high. The complete final report at production settings takes
  ≈ 1 min (measured off-laptop, torch 2.12); on the laptop the short run's report (reduced settings) took 1.2 min. Its L-BFGS
  (100 its at batch 768) took 3.5 min, i.e. ≈ 0.6 h for 4 × 100 its at batch 2048 by batch scaling. Corrected
  laptop total ≈ 7.3–8.1 h + 0.6 h + < 0.1 h ≈ **8–9 h**.
* Code compatibility with upstream `7b8eb49`/`1383f80` (Ed Hirst, s4×s1×s1 support) was re-verified on
  2026-10-01, before our commit `09cee7d` (see the "Code version" note at the top).

## S⁵×S¹ — final results (2026-10-06)

**Setup.**
* Schedule: 80k Adam steps (lr 1e-3 → 1e-5 cosine, batch 192), then the chunked L-BFGS finisher (4 × 100
  iterations at batch 2048, reverted if the held-out loss worsens).
* Network: 256×4 GELU, 2 stages, 233k parameters; λ_D = 1 unless stated.
* Hardware: run off-laptop, 8 CPU threads per run, torch 2.12.0+cpu. The code is `09cee7d` plus a patch that
  only makes checkpoint writes atomic (no effect on numerics). Speed ≈ 0.24 s/step, so 5.2–5.6 h per
  80k-step run.
* `s5xs1_sob_b7` is the laptop run resumed from its step-22000 checkpoint; the laptop's copy is preserved as
  `runs/s5xs1_sob_b7_laptop_step22000`. All other runs are fresh.
* Reproducibility, checked before the runs: re-evaluating the step-22000 checkpoint off-laptop reproduced the
  laptop's evaluation in 86 of 87 numbers to ≤ 5·10⁻¹⁴ relative. The exception is the paper chart protocol:
  since 7b8eb49 it samples through numpy's unseeded global RNG, so it carries ±2% Monte Carlo noise.
* Kill/resume of a checkpoint and bitwise determinism of identical runs were also verified.

**Main result — main's target (c = 1/(2π)), two bases × two seeds** (main's exact held-out protocol: overlap
band, seed 99999, mean of the 4 patches):

| S⁵×S¹ | base 7 (resumed) | base 8 | base 7, seed 17 | base 8, seed 18 |
|---|---|---|---|---|
| J-fit E‖δJ‖²_F | 1.02e-6 | 7.1e-7 | 8.0e-7 | 7.4e-7 |
| ∂J-fit E‖δ∂J‖²_F | 4.4e-4 | 2.9e-4 | 3.4e-4 | 3.0e-4 |
| **flat-chart ‖N‖ mean** | **0.072** | **0.060** | **0.063** | **0.061** |
| paper chart protocol ‖N‖ | 0.073 | 0.061 | 0.063 | 0.062 |
| rel. value error RMS / sup / adversarial sup | 1.50e-4 / 3.4e-4 / 4.4e-4 | 1.26e-4 / 2.7e-4 / 3.5e-4 | 1.32e-4 / 2.8e-4 / 4.0e-4 | 1.27e-4 / 2.8e-4 / 4.3e-4 |
| rel. derivative error RMS | 4.6e-4 | 3.8e-4 | 4.0e-4 | 3.8e-4 |
| intrinsic ‖N‖ mean / sup / adversarial sup | 0.029 / 0.067 / 0.086 | 0.0245 / 0.053 / 0.073 | 0.025 / 0.055 / 0.068 | 0.024 / 0.054 / 0.078 |
| orientation fix of the base | none | −J* | none | −J* |

**Flat ‖N‖ = 0.064 ± 0.005** (mean ± sample sd over the 4 runs). The per-patch values agree within 2%, so no
chart-dependent defect.

How this compares:

| reference | its value | ratio to 0.064 |
|---|---|---|
| this framework, short 8k run (2026-09-29) | 0.39 | 6× |
| old-framework results table, S⁵×S¹ CS | 6.42 | 100× |
| chart-framework Sobolev "architecture ceiling" | 1.14 | 18× |
| main's best supervised result on the lower-dimensional S³×S¹ | 0.73 | 11× |

For comparison, S³×S¹ in this framework (40k steps, 43k parameters) gave 0.017–0.026.

Further observations:
* **Not a floor.** ‖N‖ decreases monotonically through the annealing tail; in every run the mean over the
  last 5k steps is ≈ 9% below the 5k before. The endpoint is set by the schedule.
* **Robust.** Two independent certified bases, one of them used with the opposite sign, and two network seeds
  agree within ±8%.
* **The S³×S¹ mechanism carries over to d = 6.** rel. der / rel. val = 3.0–3.1 in every Sobolev run (3.4 at
  c = 1), the same O(1) value→derivative gain as on S³×S¹; flat ‖N‖ ≈ 3.5·√(∂J-fit) (S³×S¹: ≈ 3).

**Variants** (same schedule):

| run | flat ‖N‖ | intrinsic ‖N‖ mean (adv. sup) | rel. value error RMS (adv. sup) |
|---|---|---|---|
| round-orthogonal target c = 1, base 7 | 0.0160 | 0.0053 (0.015) | 1.83e-4 (5.4e-4) |
| round-orthogonal target c = 1, base 8 | 0.0142 | 0.0047 (0.013) | 1.64e-4 (4.6e-4) |
| λ_D = 3, base 7 | 0.076 | 0.030 (0.095) | 1.57e-4 (5.5e-4) |
| value-only supervision, base 7 | 1.27 | 0.48 (**905**) | 3.3e-3 (**2.02**) |

* **c = 1 against main's target, at comparable relative accuracy:** 4.3× smaller flat ‖N‖ and ≈ 5× smaller
  intrinsic ‖N‖ (S³×S¹, matched 10k runs: 8.3×). Main's anisotropic λ = e structure is intrinsically harder:
  its Dirichlet energy is 166, against 16 for the orthogonal one.
  (The c = 1 runs' main-protocol J-fit / ∂J-fit, 21.5 / 1.0e3, are measured against main's c = 1/(2π)
  structure, so they are not errors.)
* **λ_D = 3 brings no improvement** (0.076, against 0.063–0.072 for base 7): λ_D = 1 is adequate.
* **Value-only supervision fails in d = 6.**
  * It is fine on average (relative value error 3.3e-3).
  * The adversarial search finds a localised defect: relative value error 2.02, i.e. a strongly
    non-orthogonal conjugation, with ‖N‖ up to 905 there.
  * The L-BFGS finisher reverted 3 of its 4 chunks.
  * On S³×S¹ value-only supervision still gave 0.064; in d = 6 the derivative term is essential, not merely
    helpful.

**Phase B on S⁵×S¹** — Nijenhuis-only continuation (lr 1e-3, 10k steps + L-BFGS) from the base-7 and base-8
endpoints (`runs/s5xs1_phaseB_b{7,8}_lr1e-3`):

| step | 0 | 250 | 1000 | 4000 | 7000 | final |
|---|---|---|---|---|---|---|
| rel. distance to J_hopf (b7 / b8) | 0 / 0 | 0.799 / 0.801 | 0.816 / 0.818 | 0.817 / 0.824 | 0.819 / 0.823 | **0.817 / 0.820** |
| intrinsic ‖N‖ (b7 / b8) | 0.029 / 0.025 | 0.22 / 0.17 | 0.19 / 0.17 | 0.11 / 0.091 | 0.040 / 0.022 | **0.0114 / 0.0105** |

* Final values:
  * intrinsic ‖N‖ adversarial sup 0.027 / 0.029;
  * main-protocol flat ‖N‖ 0.034 / 0.032, with J-fit against J_hopf 22.0 / 22.2.
* All L-BFGS chunks were reverted.
* The drift is faster than on S³×S¹ (0.80 after 250 steps, against 0.645 there). ‖N‖ rises transiently to
  ≈ 0.2 during the jump, then falls 2.5× below its supervised starting value.
* **Identification** (`extra/ambient_phaseB_identify.py RUN_DIR BASE_DIR --case s5xs1`, generalised from n = 2
  on 2026-10-06; on the S³×S¹ endpoint it reproduces the numbers above exactly):

  | | b7 / b8 | reference |
  |---|---|---|
  | ‖J‖_F | 2.462 / 2.461 | √6 = 2.449, so within 0.5% of orthogonal |
  | Dirichlet energy | 16.35 / 16.30 | J_{c=1}: 16.00; target: 166 |
  | closest Hopf member | c′ = 1.06, distance 0.148 / 0.122 | to J_1: 0.147 / 0.125 |
  | Hopf-frame components | constant to 9.8% / 8.5% | |
  | average K̄ | ‖K̄²+I‖ 1.6e-3 / 2.2e-3, ‖K̄K̄ᵀ−I‖ 3.9e-3 / 4.9e-3 | an orthogonal complex structure |

  ⇒ **Same picture as on S³×S¹, only cleaner:** the unsupervised flow leaves main's anisotropic Hopf structure
  (λ = e) for the round-orthogonal Hopf family (λ = e^{2π}). It lands closer to that family than on S³×S¹
  (0.12–0.15 against 0.23), and stays integrable throughout.

**Conclusion (S⁵×S¹).**
* The overlap-free ambient framework reproduces main's Hopf structure on S⁵×S¹ to a relative error of
  ≈ 1.4·10⁻⁴ in value and ≈ 4·10⁻⁴ in derivative.
* Flat ‖N‖ = 0.064 ± 0.005 over 4 independent runs: 100× below the old framework's table, 18× below the chart
  framework's Sobolev ceiling, and 11× below main's best supervised result on the lower-dimensional S³×S¹.
* It is still improving at the end of the schedule.
* Derivative supervision is essential in d = 6, and the round-orthogonal target is ≈ 4–5× easier.
* Phase B flows to the round-orthogonal Hopf family, as on S³×S¹.
