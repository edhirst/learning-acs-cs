# Towards an explicit diffeomorphism Y ≅ S⁴ × T² (working notes, 2026-10-06 → 2026-10-07)

> **How to read this file.**
> * **Sections 1–10** are the chronological log, with plans written before they were implemented. Where a plan
>   was later changed, the old text stays and is marked *[superseded]* or *[corrected 2026-10-07]*.
> * **Section 11** described the state on the morning of 2026-10-07 (numpy code, one open issue).
> * **Section 12 is the authoritative description of the finished numerical construction** (torch code). The
>   diffeomorphism F: Y → S⁴ × T² is built and verified; J_Y is transported to the round S⁴ × T²; its Nijenhuis
>   tensor is computed with the repository's ambient code; and `J_round(q)` evaluates it anywhere. §12.12–12.13
>   analyse what a rigorous (computer-assisted) proof would need. Where §11 and §12 differ, §12 holds.
> * **Section 13 is the current programme** (agreed 2026-10-07): an explicit, sound, *symbolic* description of F
>   and of its inverse G, giving a purely symbolic computer-assisted proof. **Resume from §13**, using §12 as the
>   reference for the existing construction.

Y = the compact complex threefold of [HSSS26] (manuscript on branch `agentic_ledger`). The paper proves
Y ≅ S⁴×T² non-constructively: Y → W = Y/T² is a trivial principal T²-bundle because H²(W;ℤ²) = 0,
and W ≅ S⁴ follows from Wall + s-cobordism, or from Iwase's classification. Goal: replace both steps by
explicit maps, so that J_Y can be transported to the round S⁴ × T².

## 1. Key observation: a second holomorphic circle action

Every identification map in the construction has the form (ζ, w) ↦ (m ζ, w + u), with m, u depending only
on the base. This covers:
* g₂, g₃, g₄, the toric G₂, G₃ and the ℂ*-action on 𝒳;
* the log transform Φ (m = e^{iπ s(s−1)}, u = 2ts);
* the map Θ at ∞, since Θ(ζ, w, t) = (ζ e^{w/t}, …) is linear in ζ;
* ψ, which commutes with the ℂ*-action on 𝒳.

Hence the rotation R_θ: ζ ↦ e^{iθ} ζ, which is translation by θ/2π · e₁ (e₁ is the common vanishing
cycle), is a global holomorphic S¹-action on Y. It commutes with the T²-action, giving an effective
holomorphic T³ = S¹_a × T² acting on Y.

The induced action of S¹_a on W = Y/T²:
* **Smooth fibres.** On X_t/T² ≅ ℝ²/ℤ⟨e₁, e₄⟩ the action is translation along e₁, which is free
  because e₁ is primitive.
* **Singular fibres.** On I₁ × E it rotates I₁ and fixes exactly its node. In the toric chart
  λ·(x, y) = (λ⁻¹x, λy), the fixed point is x = y = 0.
* ⇒ the action on W is **semi-free with exactly two fixed points**, p₀ (over t = 0) and p_∞ (over t = ∞),
  and has no exceptional orbits.

**Weights, with W oriented by or(Y) = or(W) ⊕ or(T²).** The slice {w = 0} is the complex surface 𝒳, and
the T²-directions ē₂, ē₃ are positively oriented in ℂ_w iff ω(ē₂, ē₃) > 0.
* At t = 0, ω = 1 − |t|⁴ > 0, so (x, y) is positively oriented: the weights are (−1, +1).
* At ∞, ω = |s|⁴ − 1 < 0, so (x̄', y') is positively oriented: the weights are (+1, +1).

The round S⁴ ⊂ ℂ² × ℝ with λ·(z₁, z₂, x) = (λz₁, λz₂, x) has weights (1, 1) at N and (−1, 1) at S.
The data match, with p_∞ ↦ N and p₀ ↦ S.

## 2. The orbit space B = Y/T³ = W/S¹_a is the Hopf 3-sphere

Write x ∈ ℝ⁴ for the real lattice coordinates: (z, w) = Σ x_j e_j(t). T³ translates x₁, x₂, x₃, so x₄ mod 1
is T³-invariant. It is also invariant under the monodromy (e₄ ↦ e₄ + e₁ leaves x₄ fixed), so over ℂ* the
orbit space is B_mid = (annulus) × S¹_{x₄}.

The circle fibration β: B → ℙ¹ (induced by f) is a circle bundle. The fibres over 0 and ∞ are the images
of I₁, i.e. D₀/S¹ with the two ends identified at the fixed point.

**Cusps.** The zero section [ζ = 1, w = 0, t] extends over both cusps and has x₄ ≡ 0. On the slice w = 0,
x₄ = log|ζ| / P(t) (the paper's height H/P). *[clarified 2026-10-07: this is the e₄-coordinate of the slice
point in the filling, before the log transform. The middle coordinate is x₄^mid = log|ζ|/P₀(t) + ξ₄(t), and at
∞ x₄^mid = −log|ζ'|/P_∞(s); both verified to 9·10⁻¹⁶.]* Θ sends (e₁, e₂, e₃, −e₄) to the basis at ∞, so x₄^∞ = −x₄.
Before the log transform, β is therefore the trivial bundle S² × S¹.

**Log transform.** Φ translates by σ̃(s) = Σ ξ_j e_j with σ̃(s+1) − σ̃(s) = e₄. Continuing around
|t| = const therefore increases ξ₄ by exactly 1: Φ regluess the 0-cusp solid torus by a Dehn twist of
degree 1 in x₄. So β has Euler number ±1, i.e. B ≅ S³ with β the Hopf fibration.

This gives an alternative proof of W ≅ S⁴. W is simply connected with a semi-free circle action having two
fixed points and orbit space S³, so by the classification of circle actions on simply connected
4-manifolds (Fintushel 1977) it is equivariantly diffeomorphic to S⁴ with the action above. The point here
is that this route is constructive.

## 3. Construction plan (each step verified numerically)

*[superseded in part: steps 2–3 below (a Morse function X plus weight-1 Z on level sets) were replaced by an
orbit-space map Λ: B → S⁴/S¹ followed by an equivariant lift through one weight-1 function E (§9–§11).
Steps 1, 4 and 5 stand; the exact formula for step 4 is in §11.7.]*

1. Global holomorphic model of Y: charts, group actions, gluings, the S¹_a- and T²-actions. Check that S¹_a
   commutes with everything, and check the fixed points, weights, x₄-invariance and the clutching degree.
2. Explicit Morse function X on W: S¹-invariant, minimum −1 at p₀ and maximum +1 at p_∞, standard
   quadratic form at both, no other critical points.
3. Explicit weight-1 functions Z = (Z₁, Z₂) on W with |Z|² = 1 − X², nonvanishing off the fixed points,
   restricting to an equivariant diffeomorphism on each level set. Then Ψ_W = (Z, X): W → S⁴.
4. Explicit S¹_a-invariant T²-trivialisation φ: Y → T². Use (x₂, x₃) in the middle, and the slice
   coordinate near the cusps, corrected by a ψ-equivariant smooth height.
5. F = (Ψ_W ∘ π, φ): Y → S⁴ × T². Verify numerically that it is a local diffeomorphism of degree 1. Then
   push J_Y forward and compute N on the round S⁴ × T² with the ambient code. J is T³-invariant there,
   i.e. it is determined by data on S⁴/S¹ = S³.

## 4. Step 1 verified (`extra/s4t2_explicit/check_step1.py`)

| check | result |
|---|---|
| S¹_a commutes with g₂, g₃, g₄, T², Φ, Θ (relative deviation) | ≤ 4.3·10⁻¹⁶ |
| S¹_a on lattice coordinates: x₁ += θ/2π | 4·10⁻¹⁶ |
| T²: (x₂, x₃) += (c₂, c₃) | 10⁻¹⁵ |
| g₄: x₄ += 1 | 9·10⁻¹⁶ |
| branch change of log t: x₄ fixed, x₁ −= x₄ | 4·10⁻¹⁶ |
| clutching of Φ: ξ₄ gains exactly +1 around \|t\| = r (r = 0.01, 0.3, 3) | +1.000000000000 |
| Θ: x₄^∞ = −x₄ and (x₂, x₃)^∞ = (x₂, x₃), both mod ℤ | 6·10⁻¹⁵ |
| zero section | x ≡ 0 |

So (x₂, x₃) mod 1 is a T²-trivialisation over the whole regular part, and it is compatible with Θ.

## 5. Refined plan for Ψ_W: W → S⁴ (S¹_a-equivariant)

**(a) Orbit-space coordinates on the cusp solid tori.** On the slice w = 0 of the filling, define the
ψ̂-equivariant height

  μ̃ = Σ_n n w_n,   w_n ∝ |ζ|^{−2n} e^{n(n+1) P(t)}   (P = the paper's height function).

Since log|ζ| ∘ ψ̂ = log|ζ| + P(t), this gives μ̃ ∘ ψ̂ = μ̃ + 1 exactly *[corrected 2026-10-07: the sign is +1]*. It is a softmax over the vertices of
the periodic moment polygon, i.e. the moment map of the ψ-periodic theta-function Kähler potential
log Σ_n |ζ|^{−2n} e^{n(n+1)P}. Its properties:
* it is smooth across the central fibre;
* near the node P_n it equals n + a|x_n|² − b|y_n|² + O(4), with a, b > 0;
* it is strictly monotone along every orbit-space fibre.

So B₀ = D_ε × (ℝ/ℤ) with coordinates (t, μ̃), and similarly B_∞. On the overlap with B_mid,
μ̃ ≈ x₄ − ½ (degree 1).

*[superseded: (b) and (c) below are the first plan. The implemented versions are in §10–§11. There the
base map is T(t) = t·k(|t|), b_∞ ↦ (−1, 0), and the target linear models are as in §11.4.]*

**(b) Explicit Hopf identification.** Map B → S³ ⊂ ℂ² by matching the two solid tori
{|t| ≤ 1}, {|t| ≥ 1} with V₁ = {|u₁| ≤ |u₂|}, V₂ = {|u₁| ≥ |u₂|}. The angle discrepancy on |t| = 1 is
null-homotopic (the Euler numbers agree; check the sign), so it is blended by a cutoff.

**(c) Target.** S⁴ ⊂ ℂ² × ℝ with the standard action. Near the fixed points, Ψ is the linear model,
up to a positive scale:
* (z̄₁, z₂) = (x, y) at p₀ ↦ S;
* (z₁, z₂) = (y', x̄') at p_∞ ↦ N.

Away from them, Ψ is lifted from the orbit-space map. The two are interpolated equivariantly in thin shells.
This is legitimate because, in the first-order Hopf coordinates v = (2√(ab) t̄, μ̃ − n), the orbit-space map
is arranged to be conformal at b₀, b_∞. Hence the lifts are C¹-close in the shell and a straight-line
interpolation in ℂ² followed by normalisation is an equivariant diffeomorphism.

## 6. Step 2a verified (`extra/s4t2_explicit/cusp.py`, `check_step2.py`)

*[corrected 2026-10-07: `check_step2a.py` and `check_step2b.py` were merged into `check_step2.py`, which
checks both cusps.]*

The cusp height μ̃:

| check | result |
|---|---|
| μ̃ ∘ ψ̂ − μ̃ − 1 | 1.2·10⁻¹⁴ |
| torus formula vs toric-chart formula (chart k ± 1) | 7.5·10⁻¹⁵ |
| central fibre: D_k ↦ (k−1, k], D_{k+1} ↦ [k, k+1), node P_k ↦ k | exact |
| near P₀: μ̃ = |x|² − |y|² + O(4) | remainder constant ≤ 1 |

* ∂μ̃/∂log|ζ| < 0 everywhere. It tends to 0 only at the node, where the orbit-space coordinate is quadratic
  upstairs, as it should be.
* |μ̃ − (x₄ − ½)| ≤ 0.375, i.e. the two are related with degree 1.

## 7. The weight-1 function E (needed for the lift)

**Cusp 0.** In the toric chart U₀ (ζ = y, t = xy), E₀ = y + x̄ has weight 1. It vanishes exactly on
{|x| = |y|, arg t = π}: in B this is an arc ending at the node b₀, and E₀ is the linear model there.

The ψ̂-invariant extension is

  E_cusp = Σ_n (ψ̂ⁿ)*(χ(μ̃) · E₀),

where χ is a cutoff in μ̃ supported inside the slab of P₀, with Σ_n χ(μ̃ + n) > 0. On torus points:

  E_cusp(ζ, t) = Σ_n χ(μ̃ + n) · (mt)ⁿ ζ · (1 + t̄ / |(mt)ⁿ ζ|²).

Consecutive terms have aligned phases on the overlaps. For example, on D₁ the terms are χ₀ x̄₀ and χ₁/x₀,
both positive multiples of 1/x₀. So no cancellation occurs, and the zero set is the single arc.

**Cusp ∞.** The same construction with the ∞-data. *[corrected 2026-10-07: the local model at ∞ is now
E₀ = y + |ν| x̄, so that its zero set is exactly where the node height |ν|²|x|² − |y|² vanishes; see §10.]*

**Middle.** *[superseded: implemented as the theta blend of §8, with exact log-ratios.]* Use the theta-type function e^{2πi x₁} · f(arg t/2π − ½, x₄ − x₄⁰), where f is a degree-1 theta
function on the torus, with its zero placed on the same arc. Glue the three pieces through the logarithm of
the ratio of two functions with the same zero arc, which is nonvanishing and has winding 0.

## 8. Status (2026-10-06): a global explicit weight-1 function E on W

Code: `cusp.py` (generic cusp class), `glue.py` (lattice ↔ cusp-slice conversions), `middle.py` (global E_W).

Checks: `check_step2.py` (both cusps), `check_glue.py`, `check_middle.py`.
* **Cusp ∞** (ν = e⁻²). Height μ̃: equivariance 6·10⁻¹⁵; torus = chart; near the node
  μ̃ = |ν|²|x|² − |y|² + O(4). E: weight 1 to 4·10⁻¹⁵, ψ̂-invariant to 10⁻¹³, smooth across the central
  fibre. One zero per orbit-space torus at arg t = π, winding −1. Both cusps behave the same.
* **Conversions** lattice coordinates (t, x₁, x₄) → cusp slices, via Φ⁻¹ (resp. Θ) followed by a T²-move.
  Consistent with weight 1, the e₁- and e₄-periods and the branch change, to ≤ 6·10⁻¹³.
* **Degeneration near |t| = 1.** The slices degenerate there (max |F| ≈ 6·10¹³ at |t| = 0.97). So the cusp
  formulas are used only for |t| ≤ 0.6 (resp. ≥ 1/0.6). They are read at the plateau radii 0.7 and 1/0.7, and
  bridged by Θ_{z(r)} · exp(ℓ(r)).
  * The zero location z and the periodic log-ratio ℓ are interpolated, with ℓ evaluated exactly. The spectral
    approximation only selects the 2π branch (margin 0.13 ≪ π); the exact reconstruction error is 5·10⁻¹⁶.
  * **Topological consistency:** both reference windings are (a, b) = (0, 0). The zeros are at σ = ½ (arg t = π),
    with x₄ = 0.6655 (|t| = 0.7) and x₄ = 0.8275 (|t| = 1/0.7).
  * E_W is continuous at the four transition radii (≤ 7·10⁻⁹, i.e. bisection-limited) and has the twisted
    periodicities to 10⁻¹⁴. In the blend, |E/Θ_{z(r)}| ∈ [1.9, 270] > 0, so there is **exactly one zero** on each
    torus.

⇒ E is a smooth weight-1 function on W whose zero set is a single S¹-invariant 2-sphere π⁻¹(Γ_E) through
p₀ and p_∞. Γ_E is an arc in B along arg t = π. Near p₀, E = y + x̄ exactly (χ = 1 there); near p_∞,
E = y' + x̄'. *[corrected 2026-10-07: after the change of local model, E = y' + |ν|x̄' near p_∞. The two
reference numbers at |t| = 1/0.7 above were computed with the old model; the current values are in §11.6.
The caveat about smoothness of the blend at its zero is in §10.]*

## 9. Remaining construction (design)

**Lift formula.** For (v, X) = Λ(π(p)) ∈ B_S = S⁴/S¹ = {(v, X) : |v| + X² = 1}, set

  Z₁ = √((|v| + v₃)/2) · E/|E|,   Z₂ = (v₁ − i v₂)/(|v| + v₃) · Z₁,   Ψ(p) = (Z₁, Z₂, X).

Ψ is equivariant, and on W∖Fix it is smooth iff
* Λ⁻¹{z₁ = 0} = Λ⁻¹{v₁ = v₂ = 0, v₃ < 0} = Γ_E, and
* the transverse orientations agree.

It is a diffeomorphism iff the orbit-space map Λ: B → B_S is a homeomorphism (a diffeomorphism off the fixed
points).

**At the fixed points.** Smoothness there forces Λ to be conformal to first order in the Hopf coordinates
* v₀ = Hopf(x̄, y) = (2t̄, |x|² − |y|²) at p₀, and
* v_∞ = Hopf(x̄', y') = (2s̄, |x'|² − |y'|²) at p_∞.
  *[corrected 2026-10-07: implemented as v_∞ = Hopf(|ν|x̄', y') = (2|ν|s̄, |ν|²|x'|² − |y'|²), because at ∞
  the node height is |ν|²|x|² − |y|².]*

Orientation-preserving conformal maps are rotation × scale, i.e. they are induced by unitaries on ℂ². In a small
shell, Ψ is then interpolated with the exact linear model.

**Construction of Λ.** Λ = P_rad⁻¹ ∘ M ∘ R ∘ h_B, where:
* **h_B: B → S³_u** is the Hopf identification. It uses the base map T(t) = t·k(|t|), with k = 4π near 0 and
  1/(4π) near ∞, and r·k(r) increasing. *[corrected 2026-10-07: implemented with k = e²/(4π) = 1/(4π|ν|)
  near ∞.]* The fibre coordinates are μ̃_new, made exactly equal to |x|² − |y|²
  near the nodes. The degree of the clutching fixes the orientation of T.
* **R** is a fibre rotation u ↦ e^{−iγ}u that straightens h_B(Γ_E) to a standard arc.
* **M** is an explicit conformal (Möbius) map of S³ with M(b₀) = S, M(b_∞) = N, sending the standard arc to the
  z₁ = 0 meridian. *[superseded: no Möbius arc lies inside {u₁ū₂ ∈ ℝ_{<0}}, so this was implemented
  differently: b_∞ ↦ (−1, 0), a fibre rotation R_γ, stereographic projection from (−1, 0) and P_final;
  see §10–§11.]*
* **P_rad⁻¹** is the inverse radial projection onto B_S. *[superseded by P_final, §11.]* It is smooth upstairs because Λ is conformal at the
  poles: |w| is then a constant multiple of |v| = |z|².

**Then:**
* the T²-trivialisation φ (section 3, item 4);
* the assembly F = (Ψ∘π, φ);
* numerical verification: Jacobian bounded away from 0, degree 1, and N of F_*J_Y on the round S⁴ × T².

## 10. Status 2026-10-07: B ≅ S³ explicit; Λ: B → S⁴/S¹ built; one compatibility layer missing

**Hopf identification h_B** (`hopf.py`, `check_hopf.py`):
* fibre coordinates a₁ (|t| ≤ 1) and a₂' (|t| ≥ 1), using node-corrected heights μ̃_new;
* clutching δ(t) = arg t/2π − ξ₄(t), continuous across the cut to 3·10⁻¹⁶;
* base map T(t) = t·k(|t|) with |T| strictly increasing: d log|T|/d log r ∈ [0.23, 1.0]. An earlier choice
  with the transition on [−2, 2] folded; the transition now runs over log|t| ∈ [−4, 4].

| check | result |
|---|---|
| overlaps c0/mid/ci | ≤ 4·10⁻¹⁵ |
| the two Hopf charts agree on \|t\| = 1 (correct Euler class) | 2·10⁻¹⁵ |
| fibre coordinates strictly increasing along every fibre | yes |

* b₀ ↦ (0, 1) and b_∞ ↦ (−1, 0); the latter via a ½ shift of a₂' near ∞.
* Conformal at both fixed points with scale 2π (all three singular values 6.2832). The determinant is −1 at b₀
  and +1 at b_∞, as predicted by the orientation analysis.

**Zero arc Γ_E** (`curve.py`, `check_curve.py`):
* E has exactly one zero per orbit-space torus, lying exactly on arg t = π (to 4·10⁻¹⁶) at all scales
  10⁻⁶ ≤ |t| ≤ 10⁶.
* For this to work, the local model of E at the cusp ∞ is y + |ν|x̄, and b_∞ ↦ (−1, 0). Then Γ_E leaves both
  fixed points exactly along the real great circle through (0, 1) and (−1, 0).

**Λ = P_final ∘ moebius ∘ R_γ ∘ h_B** (`lam.py`, `check_lam.py`):
* R_γ rotates the Hopf fibres over |T| by the exact fibre phase of Γ_E there (Newton). This puts Γ_E onto the
  real great semicircle; R_γ is the identity near both fixed points.
* moebius is the stereographic projection from (−1, 0), translated; it is evaluated stably near (−1, 0).
* P_final has |v| = |w| near 0 and |v| = 1/|w| near ∞.
* ρ_O has det −1, which makes Λ orientation-preserving.

| check | result |
|---|---|
| Λ(Γ_E) on the z₁ = 0 meridian | 4·10⁻¹⁶ |
| X along Γ_E | strictly increasing |
| at b₀, b_∞: conformal (scales 2π, π), orientation-preserving | yes |

**Open: compatibility of E with Λ along Γ_E.**
* The lift Z₁ = √((|v| + v₃)/2) · E/|E| is smooth along π⁻¹(Γ_E) iff E's linearisation on the normal slices of
  Γ_E is complex-proportional to that of Λ*(v₁ + i v₂).
* This holds near the fixed points: there everything is linear and conformal, and E = y + x̄ (resp. y + |ν|x̄).
* It is NOT automatic elsewhere. In particular, the theta-blend E = Θ_{z(r)} e^{ℓ} of the middle region is only
  continuous at its zero, because the ratio of two sections with a common simple zero is direction-dependent
  there. Symptom: Newton for Γ_E stalls at a relative residual of ~10⁻⁷.

**Fix (explicit):**
* Along Γ_E away from the fixed points, use E_new = Λ*(v₁ + i v₂) / Ē″. Here E″ is any weight-1 function
  nonvanishing near Γ_E, e.g. the theta section with its zero at arg t = 0, or y − x̄ at the cusps. E_new is
  compatible with Λ by construction.
* Glue E_new to E by the logarithm of their ratio, in shells around the fixed points and in a tube shell around
  Γ_E. The ratio is smooth and nonvanishing there, with winding 0.
* The old middle blend then never has to be differentiated across its zero.

**Remaining after that:**
1. Verify that Ψ = (Z₁, Z₂, X) is a smooth equivariant diffeomorphism W → S⁴: 4×4 Jacobian on random and
   adversarial samples, degree 1, linear models at p₀ and p_∞.
2. The T²-trivialisation φ: (x₂, x₃) in the middle; near the cusps the slice coordinate corrected by
   −μ̃_new·R(t)⁻¹ē₄, which is ψ̂-equivariant and smooth.
3. Assemble F: Y → S⁴ × T².
4. Push J_Y forward and compute N on the round S⁴ × T². This needs F⁻¹ by Newton and second derivatives of F,
   e.g. by porting to torch/JAX.

## 11. State on the morning of 2026-10-07 (numpy code) *[superseded by §12 where they differ]*

*[2026-10-07, evening: §11.1–11.5 and §11.7 remain correct descriptions of the numpy modules model, cusp,
glue and hopf, which tcore.py ports one-to-one. Superseded: the middle blend of §11.3(b) and §11.6
(it is discontinuous, §12.2); Λ of §11.3(e) (§12.3); the fix plan of §11.8, whose step 2 (E″ as a ψ̂-sum
with zero arc at arg t = 0) is impossible (§12.2) and which was replaced by the two-chart lift (§12.4);
§11.9–11.11 (done, §12).]*

### 11.1 Goal and status

**Goal.** An explicit diffeomorphism F: Y → S⁴ × T², with Y the paper's threefold, so that J_Y can be
transported to the round S⁴ × T² ⊂ ℝ⁵ × ℝ² × ℝ² and its Nijenhuis tensor computed there.

**Route.** Y carries a holomorphic T³ = S¹_a × T² (§1). F is built T³-equivariantly:
* Ψ_W: W = Y/T² → S⁴ is an S¹-equivariant lift of an orbit-space map Λ: B = W/S¹ → S⁴/S¹;
* φ: Y → T² is a T²-trivialisation;
* F = (Ψ_W ∘ π, φ).

**Status.** Λ is built and verified. The lift Ψ_W needs one more layer (§11.8) before it is smooth along the
zero arc. φ is worked out on paper (§11.7) but not coded. F and the N computation are not started.

### 11.2 Files (all untracked; nothing committed)

Code, in `extra/s4t2_explicit/` (it has an `__init__.py`; run scripts from the repo root):

| file | content |
|---|---|
| `model.py` | lattice e₁…e₄ (branch k); `lattice_coords(z, w, t, k)`; g₂, g₃, g₄; `rot` (S¹_a); `torus` (T² on Y_mid); Φ; Θ; σ̃ |
| `cusp.py` | `bump_step`, `chi(m, inner = 0.4, outer = 0.9)`, `real_coeffs`; class `Cusp` with `phi`, `P`, `multiplier`, `psi_hat`, `mu_torus` (J = 12), `mu_chart` (J = 8), `E0_torus` (y + e0coef·x̄), `E_torus`; instances `CUSP0` (ν = 1, e0coef = 1) and `CUSPINF` (ν = e⁻², e0coef = e⁻²) |
| `glue.py` | `e_vectors`, `mid_point` (lattice point (x₁, 0, 0, x₄)), `to_cusp0` (Φ⁻¹, then T²-move to w = 0), `to_cuspinf` (Θ, then T²-move to w' = 0), `E0_W`, `Einf_W` (E in lattice coordinates (t, x₁, x₄, log t)) |
| `middle.py` | global weight-1 function `E_W(t, x1, x4, logt)`: cusp formulas for \|t\| ≤ R0 = 0.6 and ≥ 1/R0; theta blend between (§11.3b) |
| `hopf.py` | Hopf identification h_B: `h_B_c0`, `h_B_mid`, `h_B_ci`; `mu_new`; a₁, a₂, δ, `k_of`, `T_of`, `xi4` |
| `curve.py` | `gamma_point(r)` (Newton for the zero of E on the torus over \|t\| = r); `moebius(u)` (stable) |
| `lam.py` | `rho_O`, `nu_fn`, `P_final`, `_arc_table` (continuation), `_newton_damped`, `gamma_arc`, `gamma_u`, `fibre_phase`, `Lambda(u, r)`, `lift(v, X, E)` |

Checks, in the same folder (status as of 2026-10-07):

| script | what it checks | status |
|---|---|---|
| `check_step1.py` | S¹_a commutes with all gluings; lattice actions; clutching +1; Θ | PASS |
| `check_step2.py` | both cusps: μ̃ and E | PASS |
| `check_glue.py` | conversions | PASS |
| `check_middle.py` | E_W | PASS; numbers in §11.6 |
| `check_hopf.py` | h_B; item 5 updated to b_∞ ↦ (−1, 0) | PASS |
| `check_curve.py` | traces Γ_E | PASS; numbers in §11.6 |
| `check_lam.py` | Λ | items 1–2 PASS; item 3 aborts on the Γ_E Newton assertion (relative residual ~10⁻⁷ in the middle; this is the open issue of §11.8) |

### 11.3 Definitions and parameters as implemented

**Y_mid.** Coordinates (ζ, w, t), with ζ = e^{iz}.
* Lattice: e₁ = (2π, 0), e₂ = (−it, 1 − t²), e₃ = (−t, i(1 + t²)), e₄ = (−i log t, 2t).
* Lattice coordinates x ∈ ℝ⁴: (z, w) = Σ x_j e_j. Changing the branch log t ↦ log t + 2πi sends x₁ ↦ x₁ − x₄.
* T²: (x₂, x₃) += (c₂, c₃). S¹_a: x₁ += θ/2π.

**Cusps (filling slice w = 0).**
* Toric chart U_k: ζ = x^k y^{k+1}, t = xy. In particular U₀: ζ = y, x = t/ζ.
* S¹_a acts by λ·(x, y) = (λ⁻¹x, λy).
* ψ̂(ζ, t) = (ν·t·e^{−(c̃₂α₂ + c̃₃α₃)}·ζ, t), with c̃ = R(t)⁻¹ē₄, so that log|ζ| ∘ ψ̂ = log|ζ| + P(t), where
  P = log|ν| + log|t| − φ_t(ē₄).
* Data:
  * cusp 0: α₂ = t, α₃ = −it, ē₂ = 1 − t², ē₃ = i(1 + t²), ē₄ = 2t, ν = 1;
  * cusp ∞: α₂ = s, α₃ = is, ē₂ = s² − 1, ē₃ = i(1 + s²), ē₄ = −2s, ν = e⁻².

**Gluings.**
* A U₀-point q is the Y_mid point Φ(q); a U_∞-point corresponds to Θ⁻¹ of it.
* Exact fibre-coordinate relations (verified to 9·10⁻¹⁶):
  * x₄^mid = log|ζ|/P₀(t) + ξ₄(t) on the cusp-0 slice;
  * x₄^mid = −log|ζ'|/P_∞(s) on the cusp-∞ slice;

  where ξ = lattice coordinates of σ̃(t), i.e. of Φ's translation (`hopf.xi4`).

**(a) Heights.**
* μ̃ = Σ n w_n with w_n ∝ |ζ|^{−2n} e^{n(n+1)P}. Then μ̃ ∘ ψ̂ = μ̃ + 1, and μ̃ = k at the node P_k.
* Near P₀: μ̃ = a₀|x|² − |y|² + O(4), with a₀ = |ν|².
* μ̃_new (`hopf.mu_new`) equals exactly n + a₀|x₀|² − |y₀|² for rad := a₀|x₀|² + |y₀|² ≤ NODE_IN = 0.02, blended back to
  μ̃ by rad = NODE_OUT = 0.05. Here (x₀, y₀) are the U₀-coordinates of ψ̂^{−n}p, with n = round(μ̃).

**(b) The weight-1 function E.**
* Cusps: E = Σ_{n ∈ {n₀−1, n₀, n₀+1}} χ(μ̃ + n) · E₀(ψ̂ⁿp), with n₀ = −round(μ̃) and E₀ = y₀ + e0coef·x̄₀.
* Middle (R0 = 0.6, RP = 0.7, R1 = 0.8):
  * |t| ≤ R0: E = E0_W. R0 < |t| ≤ R1: E0_W read at the reparametrised radius h0(r), which reaches RP at R1.
  * R1 < |t| < 1/R1: E = e^{2πi x₁} · Θ_{z(r)}(σ, x₄) · exp(2πi(aσ + bx₄) + ℓ(r)), where:
    * σ = arg t/2π along a continuous branch;
    * Θ_{s₀,x₀}(s, x) = e^{2πi s(x₀ − ½)} f(s − s₀ + ½, x − x₀ + ½), with f(s, x) = Σ_n e^{−π(s+n)²/κ} e^{−2πi n x} and κ = 0.5;
    * z(r) = (1 − λ) z₀ + λ z_∞, with λ = bump in log r over [log R1, log 1/R1];
    * ℓ = (1 − λ) ℓ₀ + λ ℓ_∞, the exact logs of F_cusp/(Θ_{z_ref} e^{2πi(aσ + bx)}) at the reference radii RP and 1/RP. The 2π branch is selected by a 256² spectral approximation.
    * Both reference windings are (a, b) = (0, 0).
  * |t| ≥ 1/R0 and the reparametrised band: symmetric to the cusp-0 side.

**(c) Hopf identification h_B** (RA = 0.3, RB = 0.5).
* a₁ = μ̃_new⁰ for |t| ≤ RA, blended with `wrap` to x₄^slice − ½ by RB; in the middle, a₁ = x₄ − ξ₄(t) − ½.
* a₂: the same at ∞ with μ̃_new^∞; in the middle, a₂ = −x₄ − ½.
* a₂' = −a₂ + λ_∞δ + (1 − λ_∞)·½, with:
  * λ_∞ = 1 for |t| ≤ 1.2 and 0 for |t| ≥ 1.6;
  * δ(t) = arg t/2π − ξ₄(t), which is single-valued.
* T(t) = t·k(|t|): log k is blended from log(4π) to log(e²/(4π)) over log|t| ∈ [−4, 4]. Then
  d log|T|/d log r ∈ [0.23, 1.0], so |T| is strictly increasing.
* The point of S³:
  * for |t| ≤ 1: u₂ = e^{2πi a₁}/√(1 + |T|²), u₁ = T·u₂;
  * for |t| ≥ 1: u₁ = e^{2πi a₂'}/√(1 + |T|⁻²), u₂ = u₁/T.
* h_B(b₀) = (0, 1) and h_B(b_∞) = (−1, 0).

**(d) Γ_E.**
* The zero of E on the torus over |t| = r lies exactly at arg t = π.
* Representations: 'c0' for r ≤ 0.6, unknowns (arg t, log|ζ|); 'mid' for 0.6 < r < 1/0.6, unknowns
  (arg t, x₄); 'ci' for r ≥ 1/0.6, unknowns (arg s, log|ζ'|).
* `_arc_table`: 3601 radii in [10⁻⁹, 10⁹], traced by continuation. `gamma_arc` interpolates within each
  representation, then runs damped Newton (maximum step 0.05) and asserts a relative residual < 10⁻⁸.

**(e) Λ = P_final ∘ moebius ∘ R_γ ∘ h_B.**
* R_γ(u) = conj(φ(r))·u, with φ(r) = u₂/|u₂| (r ≤ 1) or −u₁/|u₁| (r ≥ 1), evaluated on h_B(Γ_E(r)).
* moebius(u) = (P₂, P₃, P₄)/(1 + P₁) − (0, 1, 0), with P = (Re u₁, Im u₁, Re u₂, Im u₂) and
  1 + P₁ = (P₂² + P₃² + P₄²)/(1 − P₁) for stability.
* P_final(w) = (e^{−ν(τ)}·ρ_O(w/|w|), sign(τ)·√(1 − e^{−ν(τ)})), with:
  * τ = log|w|;
  * ν(τ) = (1 − β)τ² + β|τ|, where β = bump((|τ| − ½)/½);
  * ρ_O(w) = (w₀, −w₂, −w₁), with det −1.

**(f) The lift.**
* Z₁ = √((|v| + v₃)/2)·E/|E| and Z₂ = (v₁ − iv₂)·(E/|E|)/√(2(|v| + v₃)); then Ψ_W = (Z₁, Z₂, X).
* Target: S⁴ ⊂ ℂ² × ℝ with λ·(z₁, z₂, X) = (λz₁, λz₂, X) and B_S = {(v, X) : |v| + X² = 1}, where
  v = (2z₁z̄₂, |z₁|² − |z₂|²).
* p₀ ↦ S (X = −1) and p_∞ ↦ N (X = +1).

### 11.4 Conventions at the fixed points

**Orientation of W.** W is oriented by or(Y) = or(W) ⊕ or(T²):
* (x, y) is positively oriented at p₀, with S¹-weights (−1, 1);
* (x̄', y') is positively oriented at p_∞, with weights (1, 1).

**Hopf coordinates.**
* v₀ = Hopf(x̄, y) = (2t̄, |x|² − |y|²). It is negatively oriented in B, as v_S is at S.
* v_∞ = Hopf(|ν|x̄', y') = (2|ν|s̄, |ν|²|x'|² − |y'|²). It is positively oriented, as v_N is at N.

**Measured differentials** (`check_lam.py`, item 2). DΛ is conformal: 2π·(rotation) at b₀ and π·(rotation) at
b_∞, both with det +1.

**Hence the linear models of Ψ_W:**
* at p₀: Ψ ≈ √(2π)·U₀(x̄, y);
* at p_∞: Ψ ≈ √π·U_∞(|ν|x̄', y');

with U₀, U_∞ ∈ U(2) whose first rows are proportional to (1, 1), because z₁ ∝ E = y + x̄, resp. y' + |ν|x̄'.

### 11.5 Verified (current code, reruns of 2026-10-07)

| check | result |
|---|---|
| S¹_a vs all gluings | ≤ 4·10⁻¹⁶ relative |
| clutching of Φ | +1.000000000000 |
| Θ: x₄^∞ = −x₄ and (x₂, x₃) preserved mod ℤ | 6·10⁻¹⁵ |
| μ̃ equivariance | ≤ 1.5·10⁻¹⁴ |
| E: weight 1 / ψ̂-invariance | ≤ 8·10⁻¹⁵ / ≤ 10⁻¹³ |
| E: one zero per torus at arg t = π, winding −1 | both cusps |
| h_B: overlaps / Hopf clutching | ≤ 4·10⁻¹⁵ / 2·10⁻¹⁵ |
| h_B: fibre derivatives | > 0 |
| h_B at b₀ / b_∞ | conformal, 2π (det −1 / det +1) |
| Λ(Γ_E) on the meridian | 4·10⁻¹⁶ |
| X along Γ_E | strictly increasing |
| Λ at b₀ / b_∞ | conformal 2π / π, det +1 |

### 11.6 Current numbers of the middle E_W and Γ_E (reruns 2026-10-07, current code)

**References** (`check_middle.py`):

| | zero | \|F(z)\| | windings (a, b) | exact-branch reconstruction | branch margin |
|---|---|---|---|---|---|
| \|t\| = 0.7 | (σ, x₄) = (0.5, 0.665524) | 3·10⁻¹⁵ | (0, 0) | 5·10⁻¹⁶ | 0.13 |
| \|t\| = 1/0.7 | (0.5, 0.724574) | 8·10⁻¹⁷ | (0, 0) | 5·10⁻¹⁶ | 0.026 |

The |t| = 1/0.7 zero was 0.8275 before the ∞ local model became y + |ν|x̄.

**E_W across the middle:**
* continuous at the transition radii 0.6 / 0.8 / 1.25 / 1/0.6: 6.5·10⁻⁹ / 0 / 0 / 5.2·10⁻⁹;
* twisted periodicities hold to ≤ 7·10⁻¹⁵;
* in the blend, |E/Θ_{z(r)}| ∈ [0.62, 31], so there is exactly one zero on each torus. All 13 radii in
  [0.55, 1.82] show one zero cluster, at σ = ½.

**Γ_E** (`check_curve.py`, 601 radii in [10⁻⁶, 10⁶]):
* zero residuals ≤ 2.7·10⁻¹⁵, and |arg t| − π ≤ 4·10⁻¹⁶ (exactly on arg t = π);
* the endpoints map to (0, 1) and (−1, 0);
* after moebius ∘ h_B, but before R_γ, the arc is exactly the ray d = (0, 1, 0) for r ≤ 10⁻³ and r ≥ 10³. In
  between it wanders: angle to d up to 150°, |w| non-monotone. That is why R_γ exists.

### 11.7 The T²-trivialisation φ: Y → T² (worked out 2026-10-07, not yet coded)

Requirements: φ(c·p) = φ(p) + c, φ invariant under S¹_a, smooth.

* **Regular part:** φ_mid = (x₂, x₃) mod 1. It is invariant under the branch change and under e₄, and it is
  compatible with Θ, since (x₂, x₃)^∞ = (x₂, x₃) mod ℤ (verified to 2·10⁻¹⁵).
* **Cusp 0, on Ũ₀ = 𝒳^ε × ℂ_w:** φ₀(p, w) = C − H·c̃(t), where:
  * C = R(t)⁻¹w ∈ ℝ², i.e. w = C₂ē₂ + C₃ē₃;
  * c̃ = R(t)⁻¹ē₄(t) = O(t);
  * H = μ̃_new at the T²-slice point of (p, w), i.e. ζ_slice = ζ·e^{−(C₂α₂ + C₃α₃)}.

  Checks:
  * G₂, G₃: φ₀ += (1, 0), (0, 1) ✓;
  * T²: φ₀ += c ✓;
  * G₄: C += c̃ and H += 1, so φ₀ is unchanged ✓;
  * φ₀ is smooth through the central fibre, because C, H and c̃ are.
* **Gluing to the middle** (U₀ ∋ q ~ Φ(q)): φ_mid ∘ Φ − φ₀ = (H − x₄^slice)·c̃(t) + (ξ₂, ξ₃)(t).
  * H − x₄^slice is G₄-invariant, real and bounded (≈ −½ + small).
  * (ξ₂, ξ₃) are single-valued: continuation only changes ξ₁ and ξ₄.
  * So the difference D is a smooth ℝ²-valued function on the overlap, and φ := φ₀ + χ(|t|)·D (mod ℤ²) is
    explicit.
* **Cusp ∞:** the same through Θ (no ξ term). Use the ∞-data and c̃^∞ = R_∞(s)⁻¹ē₄^∞.

### 11.8 The open issue and the exact fix plan

**Corrected sufficient criterion (2026-10-07; statement change).** Away from the fixed points, put
f = Λ*(v₁ + iv₂) and d = |v| − v₃ > 0 near π⁻¹(Γ_E). If f = H E with H smooth and nonvanishing
(in a local unitary frame), lift B extends by Z₁ = |H|E/√(2d), Z₂ = (H̄/|H|)√(d/2).
These are identities, checked in `sym/lemmas_1d.py` L8. Matching the transverse linearisations of E and f
alone is insufficient: E = z+z̄² with v = (2 Re(z)√(1−|z|²), 2 Im(z)√(1−|z|²), 2|z|²−1)
has df(0) = 2dE(0), but lift B is not differentiable at 0 (L8d–e; programme §2).
The earlier equivalence with first-order compatibility is withdrawn. The exact factorisation
f = (smooth nonvanishing)·E·Ē″ remains sufficient when E″ is smooth and nonvanishing near Γ_E.

**Current E.**
* It is compatible exactly at b₀ and b_∞, in the limit.
* It is not compatible elsewhere. The middle blend is only C⁰ at its zero: the ratio F_cusp/Θ_{z_ref} is
  direction-dependent at z_ref.
* Symptom: Γ_E Newton stalls (relative 10⁻⁸ to 10⁻⁷) at the blend onsets |t| ≈ 0.822 and 1.216. `check_lam`
  item 3 aborts on the assertion in `gamma_arc`.

**Fix, in order:**
1. **Make Λ exactly linear near the fixed points.** Λ_mod = DΛ + c(|v|)·(Λ − DΛ) in v₀- (resp. v_∞-)
   coordinates: c = 0 in a small ball, 1 outside a shell. DΛ is the measured conformal map (2π resp. π times a
   rotation, det +1). Γ_E is exactly the ray v₀ ∝ (−1, 0, 0) there, and both Λ and DΛ send it into the
   meridian, so Λ_mod still maps Γ_E onto the meridian. Then E = y + x̄ (resp. y + |ν|x̄) is exactly compatible
   inside the ball.
2. **Choose E″.** E″ = the same construction as E, with its zero arc at arg t = 0:
   * cusp local models y − e0coef·x̄ (zero set {μ_node = 0, arg t = 0}), with the same ψ̂-sum;
   * in the middle, the theta section with σ₀ = 0.

   E″ is nonvanishing near Γ_E except at the fixed points. Inside the linear balls, choose E″ ∝ the second row
   of U₀ (resp. U_∞) applied to (x̄, y) (resp. (|ν|x̄', y')), so that E_new is smooth at p₀, p_∞.
3. **E_new = Λ_mod*(v₁ + iv₂)/Ē″.** It has weight 1, is smooth near Γ_E, and is compatible by construction.
4. **Glue.** E_final = E_new on a tube around Γ_E, E elsewhere. In the tube shell, away from Γ_E, both are
   nonvanishing; glue by E_new·exp(c·log(E/E_new)). The ratio has winding 0 around Γ_E iff the transverse
   orientations agree; CHECK this first. Inside the linear balls, E and E_new are both compatible, so the ratio
   is smooth there too provided the exact smooth-quotient condition above holds; matching linearisations
   alone does not establish it.
5. Re-run `check_lam` (with the Γ_E assertion) and add lift checks: Ψ smooth across π⁻¹(Γ_E), by finite
   differences on both sides; the linear models at p₀, p_∞; the 4×4 Jacobian (W → S⁴ in local coordinates)
   with constant sign and bounded away from 0 on random plus adversarial samples; degree 1.

### 11.9 After the fix

1. Code φ (§11.7) and check: T²-equivariance, invariance under all gluings, smoothness at the central fibres.
2. F = (Ψ_W ∘ π, φ): Y → S⁴ × T² ⊂ ℝ⁵ × ℝ² × ℝ². Verify:
   * T³-equivariance;
   * the 6×6 Jacobian is nonsingular everywhere;
   * degree 1, e.g. by counting preimages of random points.
3. Transport J: J_round(q) = DF·J_Y·DF⁻¹ at p = F⁻¹(q). This needs F⁻¹ by Newton and second derivatives of F
   for N. Port the construction to torch (float64, autograd), or use JAX. Then compute the intrinsic ‖N‖ with
   `src/loss/ambient_loss.py`; it should be round-off. Note that J_round is T³-invariant, i.e. determined by
   data on S⁴/S¹ = S³.
4. Report and document, e.g. as a computational companion to the paper. Commit only when the user asks.

### 11.10 Pitfalls already met (do not repeat)

* **|T| = r·k(r) must be strictly increasing.** A k-transition over log r ∈ [−2, 2] folded the map; now
  [−4, 4].
* **Möbius near B' = (−1, 0).** Use the stable form 1 + Re u₁ = (|Im u₁|² + |u₂|²)/(1 − Re u₁).
* **Newton for Γ_E:**
  * interpolate initial guesses only within one representation: the l-coordinate means log|ζ| in c0/ci and x₄
    in mid;
  * use damping;
  * assert convergence.
* **E at ∞** must use the local model y + |ν|x̄. Otherwise Γ_E is not tangent to the real great circle at b_∞.
* **Relative errors.** The commutation tests need relative error, because e^{w/t} reaches e⁶⁰.
* **The cusp slices degenerate near |t| = 1** (|F| ~ 10¹³ at 0.97). Never evaluate cusp formulas beyond
  |t| = 0.6 (resp. 1/0.6) except at the plateau radii 0.7 and 1/0.7.
* **Ratios at a common zero.** The ratio of two sections with a common simple zero is direction-dependent. Never
  log-blend through a zero unless both sides are compatible (§11.8).

### 11.11 Resume checklist

1. Read §11. Run `check_step1`, `check_step2`, `check_hopf` and `check_middle` (≈ 3 min, all PASS).
   `check_lam` item 3 is expected to abort until the §11.8 fix is in.
2. Implement §11.8, steps 1 → 5, verifying each.
3. Then §11.9.

## 12. AUTHORITATIVE STATE (2026-10-07, evening): F: Y → S⁴ × T² built and verified, J on the round S⁴ × T²

### 12.0 Results in one table

| claim | evidence (script) | result |
|---|---|---|
| torch port = verified numpy code | `check_tcore.py`, 22 functions | ≤ 2·10⁻¹⁴ |
| new E, Λ, G_A well defined on W (branch of log t, x₄ + 1, weights) | `check_tlift.py` | ≤ 1.2·10⁻¹⁴ |
| chart and middle representations agree | `check_tlift.py` | ≤ 1.7·10⁻¹⁵ |
| Ψ_W: W → S⁴ lands in S⁴, S¹-equivariant, orbit map = Λ | `check_lift.py` | ≤ 1.3·10⁻¹⁴ |
| exact linear models at p₀, p_∞; seamless switch | `check_lift.py` | 1.1·10⁻¹⁶ |
| Ψ_W is C¹ across every gluing locus (jump/h independent of h) | `check_lift.py` | ratios 1.000 / 1.007 / 1.002 |
| det DΨ_W of constant sign, ≠ 0 (100k samples per chart + adversarial) | `check_lift.py` | min ≥ 1.4·10⁻² |
| shell gluing has a continuous branch: max \|Arg g\| < π | `check_lift.py` | 1.598 sampled, 1.616 adversarial |
| **deg Ψ_W = 1** (volume integral, 3 resolutions) | `check_degree.py` | 0.999949 → 0.999999 → **1.000000** (det > 0 at all 18.4M nodes) |
| F well defined, T³-equivariant, lands in S⁴ × T², charts agree | `check_full.py` | ≤ 7.5·10⁻¹³ |
| det DF has one sign in all three holomorphic charts | `check_full.py` | 60k samples per chart |
| det F = −(2π)² det Ψ_W ⇒ \|deg F\| = 1 | `check_full.py` | 5.7·10⁻¹⁴ |
| D²F by nested forward AD = finite differences of DF (all charts, incl. t = 0 and nodes) | `nijenhuis_round.py` | ≤ 7·10⁻⁸ (middle), ≤ 3·10⁻⁶ (fillings) |
| J from two charts at the same points: J and DJ agree | `nijenhuis_round.py` | 7.7·10⁻¹⁴ / 1.2·10⁻¹³ |
| **N of F_*J_Y on the round S⁴ × T²** (ambient code, 2000 samples per region) | `nijenhuis_round.py` | \|N\|/scale ≤ 3.7·10⁻¹³ (fillings, incl. t = 0 and the nodes), ≤ 3.2·10⁻¹⁰ (middle, where \|J\| reaches 1.7·10³) |
| non-integrable control J_ε through the same F | `nijenhuis_round.py` | \|N\|/scale = 3.8·10⁻³ (ε = 0.1), 3.3·10⁻⁵ (ε = 10⁻³) |
| F⁻¹ on 3000 uniform points of S⁴ × T² | `check_inverse.py` | after Newton ≤ 8.4·10⁻¹⁴ (explicit chain alone exact except in the linear-blend shells) |

### 12.1 Files (all untracked; nothing committed)

New torch code in `extra/s4t2_explicit/` (float64 / complex128, differentiable to any order; run from the repo root):

| file | content |
|---|---|
| `tcore.py` | one-to-one port of model / cusp / glue / hopf / `curve.moebius` / `lam.P_final`; plus `TCusp.E_chart` (E in U₀ coordinates, valid through the central fibre), `TCusp.reduce`, `mu_new_chart`; `lattice_coords` by Cramer's rule (§12.8) |
| `tlift.py` | new E (`E_W`, `E_chart0`, `E_chartinf`); `hB_W` (dispatch), `hB_chart0/inf`; arc model (`x4_base`, `arc_model`, `cheb_eval`); straightening `arc_phases`; `Lam_from_u`; linear models (`Q0`, `QINF`, `lin_blend`, `Lam_W`, `Lam_chart0/inf`); gauge `GA_hat_W`, `GA_hat_chart0/inf`, `psi_of_r`; lift `assemble`, `g_shell`, `linear_model0/inf`, `Psi_W`, `Psi_chart0/inf` |
| `tfull.py` | φ (`D0`, `Dinf`, `phi_mid`), `F_mid`, `F_c0`, `F_ci`, chart transitions `c0_to_mid`, `ci_to_mid`, real wrappers `F_*_real` (ℝ⁶ → ℝ⁹) |
| `tinverse.py` | F⁻¹: explicit inverse chain + Newton polish (`Finv`) |
| `jround.py` | **`J_round(q)`**: the complex structure on the round S⁴ × T² at any points q ∈ ℝ⁹, returning J (9 × 9), DJ (9 × 9 × 9) in the layout of `src/loss/ambient_loss.py`, and the inverse residual; `nijenhuis(q, J, DJ)` |
| `fit_arc.py` → `arc_fit.json` | Chebyshev model of the zero arc Γ_E (§12.3) |
| `fit_psi.py` → `psi_fit.json` | branch normaliser ψ(\|t\|) of the gauge (§12.4) |
| `check_tcore.py`, `check_tlift.py`, `check_lift.py`, `check_degree.py`, `check_full.py`, `nijenhuis_round.py`, `check_inverse.py`, `check_jround.py` | the checks of §12.0 and §12.11 |
| `symbolic_identities.py` | sympy proofs: dμ̃/dρ = −2·Var_w(n) (the cusp height is strictly monotone along fibres), and DΛ(b₀) = 2π·Q₀ exactly, det 8π³ (§12.12) |
| `check_nblind.py` | the N-check is blind to F: random other 1- and 2-jets also give N ≈ 0 (§12.12) |
| `toy_symbolic_blindness.py` | sympy toy: a two-to-one, degenerate F still gives J² + 1 ≡ 0 and N ≡ 0 identically, yet J is multivalued and unbounded (§12.13) |

The numpy modules model, cusp, glue, hopf and their checks (`check_step1`, `check_step2`, `check_glue`,
`check_hopf`) are unchanged and still pass (rerun 2026-10-07); they are the independent reference for `tcore.py`.
`middle.py`, `curve.py` and `lam.py` are superseded (marked in their docstrings) and kept as a record. Their
checks `check_middle.py`, `check_curve.py` and `check_lam.py` test the superseded code; `check_lam` item 3 aborts
by design (§11.8).

### 12.2 The weight-1 function E: the middle blend replaced

**Finding.** The theta/log blend of §8 and §11.3(b) is **discontinuous**, not merely C⁰ at its zero. On a blend
torus E = Θ_{z(r)}·exp((1 − λ)ℓ₀ + λℓ_∞), where ℓ₀ = log(F_cusp0/Θ_{z_ref0}) is evaluated at the reference zero
z_ref0 on *every* blend torus. At z_ref0 the ratio of two sections with a common simple zero is
direction-dependent, so E jumps there although z(r) ≠ z_ref0. Measured: values within 10⁻⁷ of z_ref0 spread by
up to 0.12 (|t| = 0.95). This is why Γ_E-Newton stalled at 10⁻⁷ (§11.8).

**New middle.** On R1 < |t| < 1/R1 (R1 = 0.8):

  E = (1 − λ(|t|)) · E0_W(at |t| = RP) + λ(|t|) · Einf_W(at |t| = 1/RP),  λ = bump in log|t|,

a linear combination of two sections of the same bundle, hence smooth. For R0 < |t| ≤ R1 it is E0_W read at
h₀(|t|), as before (symmetric at ∞). Exactly one zero on every blend torus (cell windings on 400² grids,
21 values of λ); the zero drifts off σ = ½ by ≤ 0.016.

**Zero arc Γ_E** (Newton, 12001 radii, relative residual ≤ 2·10⁻¹⁵): arg t = π exactly except on
(0.813, 1.232), where |arg t − π| ≤ 0.107; x₄ = ½ exactly for |t| ≤ 0.366 (there E = y + x̄ near its zero);
x₄ = the exact local-model zero at ∞ for |t| ≥ 1.65.

**No companion E″.** §11.8 step 2 (E″ = ψ̂-sum of y − e0coef·x̄, zero arc at arg t = 0) is impossible:
on the slab overlaps consecutive terms are −c·x̄₀ and 1/(m x₀) ∝ x̄₀/m, which cancel unless Re(c·m) > 0. A
ψ̂-sum of local models y + c·x̄ therefore has its zero arc at arg t ∈ (π/2, 3π/2). Measured: spurious zeros, net
winding 0 on some tori. Replaced by the gauge G_A (§12.4), which needs no ψ̂-sum.

### 12.3 The orbit-space map Λ_mod: B → B_S = S⁴/S¹

Λ = P_final ∘ moebius ∘ R ∘ h_B, with h_B, moebius, P_final as in §11.3, and:

* **Straightening R(u) = (φ̄₁u₁, φ̄₂u₂)**, with (φ₁, φ₂) = the phases of (−u₁, u₂) at h_B of the *arc model*
  at radius r = |t|. R is a diffeomorphism of S³ (r is a function of |u₁|/|u₂|, which R preserves); two
  separate phases also straighten the deviation of Γ_E from arg t = π. Windowed: R = id for |t| ≤ 0.005 and
  |t| ≥ 26 (there the arc is exactly real).
* **Arc model** (`fit_arc.py`): arg t = π + W₁(ρ)p₁(ρ), x₄ = x4_base(r) + W₂(ρ)p₂(ρ), ρ = log r.
  * x4_base = ½ blended (on [1.25, 1.6], away from \|t\| = 1 where the ∞ slice degenerates) into the exact
    ∞ local-model zero.
  * p₁, p₂: Chebyshev fits, degrees 100 / 120, max errors 1.2·10⁻¹⁰ / 5.2·10⁻⁵. Their residuals vanish
    identically outside [0.81, 1.24] and [0.366, 1.65]; the windows W₁, W₂ equal 1 on those supports.
  * The fit only affects *which* diffeomorphism R is: any smooth model gives a valid Λ.
  * *Pitfall met:* a closed-form "exact" arc (|x| = \|y\| up to R1) was wrong beyond r = 0.366, because there
    the zero's height exceeds 0.1 and neighbouring χ-terms shift it; with it Λ(Γ_E) wandered to θ_v = 112°.
* Result: **Λ(Γ_E) lies within 0.056° of the meridian M₁**; X strictly increases along Γ_E.
* **Exact linear models.** In Hopf coordinates v_q = Hopf(√a₀ x̄₀, y₀) of the reduced U₀ point (a₀ = |ν|²,
  rad = |v_q| = a₀|x₀|² + |y₀|²), Λ_mod = M·v_q exactly for rad ≤ EPS1 = 0.004, blended into Λ by
  rad = EPS2 = 0.008 (bump in rad). Derived by fitting at |v_q| ≈ 10⁻⁷; exact signed permutations, det +1:
  * M₀ = 2π·[[0, −1, 0], [0, 0, −1], [1, 0, 0]];
  * M_∞ = π·[[0, 0, 1], [0, −1, 0], [1, 0, 0]].

  Both send E's zero ray (−1, 0, 0) to the meridian direction (0, 0, −1).
* Λ_mod is a homeomorphism B → B_S (composition of homeomorphisms; diffeomorphism off b₀, b_∞): checked via
  det DΨ_W ≠ 0 (§12.5).

### 12.4 The lift Ψ_W: two charts glued on a shell

Target S⁴ = {|Z₁|² + |Z₂|² + X² = 1} with λ·(Z, X) = (λZ, X), v = Hopf(Z), and (v, X) = Λ_mod(π(p)). Write
ω_E = E/|E| and ω_A = G_A/|G_A|.

* **Lift_B (via E)**: Z₁ = √((|v| + v₃)/2)·ω_E, Z₂ = (v₁ − iv₂)·ω_E/√(2(|v| + v₃)). Smooth where E ≠ 0 and
  |v| + v₃ > 0, i.e. away from Γ_E and Λ⁻¹(M₁).
* **Lift_A (via G_A)**: Z₂ = √((|v| − v₃)/2)·ω_A, Z₁ = (v₁ + iv₂)·ω_A/√(2(|v| − v₃)). Smooth where G_A ≠ 0 and
  |v| − v₃ > 0.
* **Gluing.** θ_v is the angle between v and (0, 0, −1). c = 0 for θ_v ≤ θ_a = 15° and c = 1 for
  θ_v ≥ θ_b = 45° (C^∞ bump). Ψ_W = e^{i c Arg g}·Lift_A, where g = Lift_B/Lift_A = ω_E·ω̄_A·conj(phase(v₁ + iv₂)).
  The code evaluates Lift_A where c = 0 and Lift_B where c = 1.
* **Why this is well defined.** Lift_B is only used on {θ_v ≥ 15°}, where E ≠ 0 because Λ(Γ_E) lies within
  0.056° of M₁. Lift_A is only used on {θ_v ≤ 45°}, where G_A ≠ 0 by construction. On the shell, g has winding
  0 around the tube. On each level sphere of X, E and the target section z₁ (pulled back by either lift) each
  have exactly one zero, both with the sign of c₁ of the S¹-bundle. Arg is the
  continuous branch iff |Arg g| < π there. Measured: max 1.598 on 425k shell samples, 1.616 under adversarial
  ascent, so the margin is about 1.5 rad.
* **The gauge.** G_A = e^{iψ(|t|)}·Ĝ_A, with Ĝ_A built from three pieces:
  * y₀ − x̄₀ at cusp 0 (reduced point);
  * e^{2πix₁} in the middle, on the canonical branch arg t ∈ [0, 2π);
  * y₀′ − |ν|x̄₀′ at ∞.

  The pieces are glued by exp of weighted explicit log-ratios on |t| ∈ [0.4, 0.5] and [2, 2.5]. For example
  L₀ = 2πix₁ − log(y₀ − x̄₀), with log y₀ = i z_slice − n·log(mult(t)) and log(1 − t̄/|y₀|²) principal.
  Domain facts, measured on 400k points:
  * on the tube {θ_v ≤ 45°}, |arg t| ≥ 133.8°;
  * the cusp heights there satisfy |μ̃ mod 1| ≤ 0.32 (0) and ≤ 0.08 (∞).

  So the reduction index is constant, Log(1 − t̄/|y₀|²) is principal, and the branch cut arg t = 0 is avoided.
* **Branch normaliser ψ.** ψ = −π/2 for |t| ≤ 0.002 and π for |t| ≥ 60, exactly; these are the values making
  g ≡ 1 in the linear balls, so the glued lift equals the linear model there. In between, ψ is a degree-60
  Chebyshev fit of the circular mean of g along the tube (`fit_psi.py`). Continuation steps are ≤ 0.22 rad and
  give k = 0, i.e. ψ runs from −π/2 to +π. The parameter is log|t|, not X: the whole middle of W sits in
  X ∈ [0.86, 0.91] along the arc.
* **Linear models of Ψ_W** (smoothness at the fixed points):
  * near p₀ (rad ≤ EPS1): Ψ = √(2π)·((x̄ + y)/√2, i(x̄ − y)/√2), X = −√(1 − 2π(|x|² + |y|²));
  * near p_∞: Ψ = √π·((|ν|x̄′ + y′)/√2, (|ν|x̄′ − y′)/√2), X = +√(1 − π(|ν|²|x′|² + |y′|²)).

  The unitaries U₀, U_∞ have first row (1, 1)/√2, so Z₁ has the phase of E. Their second rows are i(1, −1)/√2
  and (1, −1)/√2, i.e. β₀ = −π/2 and β_∞ = π.

### 12.5 Verification of Ψ_W (`check_lift.py`, `check_degree.py`)

* **Basic properties** (20k points, |t| ∈ [10⁻³, 10³]):
  * |Z|² + X² = 1 to 7.8·10⁻¹⁶;
  * Hopf(Z) = v and X = X_Λ to 6.7·10⁻¹⁶;
  * S¹-equivariance to 6.6·10⁻¹⁵;
  * branch change and x₄ + 1 to 1.3·10⁻¹⁴.
* **Linear models.** The general formula equals the linear model on EPS1/2 < rad ≤ EPS1 to 1.1·10⁻¹⁶, at both
  fixed points (15k and 13k points).
* **Smoothness.** For each sample, jump(h) = |DΨ(p + h·d) − DΨ(p − h·d)|/|DΨ(p)|. Max jump/h at
  h = 10⁻⁴, 10⁻⁵, 10⁻⁶:
  * middle: 725.7 / 725.8 / 725.8, worst point at |t| = 2.13 (glue of Ĝ_A);
  * c0 chart: 3974 / 4003 / 4003 (rad = 0.0044);
  * ci chart: 5770 / 5782 / 5782 (rad = 0.0052).

  An h-independent jump/h means bounded second derivatives, i.e. C¹ with no kink. The large values sit in the
  linear-blend shells of Λ_mod and the gauge glue.
* **Jacobian.** det[Ψ, DΨ] (5 × 5), 100k samples per chart, plus 64 adversarial starts × 150 Adam steps on
  log|det|:

  | chart | coordinates | sign | min \|det\| | adversarial min |
  |---|---|---|---|---|
  | middle | (log r, arg t, x₁, x₄), \|t\| ∈ [0.2, 5] | + | 1.55·10⁻² | 1.53·10⁻² |
  | c0 | (x, y) | + | 7.99·10⁻² | unchanged |
  | ci | (x′, y′) | − | 1.43·10⁻² | unchanged |

  (x′, y′) is negatively oriented at p_∞ (§11.4), and the overall orientation is consistent (§12.6).
  On the central fibres themselves and near the nodes, with 20k points for each of the three cases at each cusp,
  det is finite and keeps the chart's sign:

  | cusp | x = 0 (central fibre) | y = 0 (central fibre) | \|q\| ~ 10⁻³ (node) |
  |---|---|---|---|
  | 0 | min 8.68 | min 8.68 | min 39.5 |
  | ∞ | min 0.170 | min 0.181 | min 0.181 |

  The ∞ chart is scaled by |ν|².
* **Degree.** deg Ψ_W = (3/8π²)·∫_W Ψ*vol_{S⁴}, a 3D quadrature over (log r, arg t, x₄); the integrand is
  S¹-invariant and decays like r^{±2}. Values:
  * grid 281 × 48 × 48: 0.999949;
  * grid 561 × 96 × 96: 0.999999;
  * grid 1121 × 128 × 128: **1.000000**, with det > 0 at all 18,366,464 nodes.

* **Conclusion.** Ψ_W is an orientation-preserving local diffeomorphism of degree 1 between compact connected
  4-manifolds, hence a diffeomorphism. W is connected and S⁴ simply connected, so a local diffeomorphism is
  a covering and therefore bijective; the degree confirms this independently. It is exactly linear near p₀
  and p_∞.

### 12.6 The T²-trivialisation φ and F (`tfull.py`, `check_full.py`)

* φ as in §11.7:
  * φ = (x₂, x₃) − (1 − χ(|t|))·D₀ − (1 − χ(1/|t|))·D_∞, with χ = bump on [0.3, 0.5],
    D₀ = (H − x₄^slice)·c̃(t) + (ξ₂, ξ₃)(t) and D_∞ = (H′ − x₄′^slice)·c̃^∞(s);
  * in the fillings, φ₀ = C − H·c̃.
* F = (Ψ_W ∘ π, φ): Y → S⁴ × S¹ × S¹ ⊂ ℝ⁵ × ℝ² × ℝ² (unit circles). Three holomorphic charts:
  * **middle (z, w, t)**: lattice coordinates by Cramer's rule;
  * **filling (x, y, w) at 0**: T²-slice C = R(t)⁻¹w, (x_s, y_s) = (x e^{C·α}, y e^{−C·α});
  * **filling (x′, y′, w′) at ∞**: Y_mid point Θ⁻¹.
* Results:
  * charts agree on the overlaps (10⁻³ ≤ |t|, |s| ≤ 0.25): 4.9·10⁻¹⁴ (c0) and 6.7·10⁻¹⁴ (ci);
  * invariance under the lattice translations e₁…e₄: 7.5·10⁻¹³; under the branch of log t: 1.5·10⁻¹⁴; under
    G₂, G₃: 5.0·10⁻¹⁴ and 6.4·10⁻¹⁴;
  * T²-equivariance: 6.5·10⁻¹³; S¹_a-equivariance: 2.2·10⁻¹³; F lands in S⁴ × T² to 3.3·10⁻¹⁶;
  * det[n_{S⁴}, n_{S¹}, n_{S¹}, DF] (9 × 9) in holomorphic coordinates has the **same sign (−) in all three
    charts** (60k samples each), so F is a local diffeomorphism compatible with the complex orientation;
  * in lattice coordinates det F = −(2π)²·det Ψ_W (to 5.7·10⁻¹⁴): F is a T²-equivariant map covering Ψ_W, so
    |deg F| = deg Ψ_W = 1.
* Hence **F is an explicit diffeomorphism Y → S⁴ × T²**, T³-equivariant, real-analytic up to the C^∞ bump
  functions.

### 12.7 J on the round S⁴ × T² and its Nijenhuis tensor (`nijenhuis_round.py`, `check_inverse.py`)

**Method.** For a sample p in a holomorphic chart (J_Y = J₀ standard there):
* A = DF(p) (9 × 6) and D²F by nested forward-mode AD.
* J = A·J₀·A⁺ (tangential, J² = −P), with A⁺ = pinv(A) (SVD).
* DJ[a, r, c] = Σ_k (∂_k J)[a, r]·A⁺[k, c], where:
  * ∂_k J = (∂_k A)·J₀·A⁺ + A·J₀·∂_k(A⁺);
  * ∂_k(A⁺) = −A⁺(∂_k A)A⁺ + (AᵀA)⁻¹(∂_k A)ᵀ(I − P(q));
  * I − P(q) is the exact normal projector of S⁴ × T².
* P and DP come from `AmbientManifold((4, 1, 1))`, and N from `src/loss/ambient_loss.ambient_nijenhuis`. The
  per-sample norm is the intrinsic round-metric norm. "scale" is the sum of the norms of the four bracket
  terms.

**Results** (final run, 2000 samples per region):
* **D²F precondition.** Nested AD agrees with finite differences of DF, all (k, l): ≤ 4.6·10⁻⁸ in the middle,
  ≤ 2.0·10⁻⁶ and 7.9·10⁻⁷ in the fillings (FD-limited near the nodes).
* **Middle region, |t| ∈ [0.25, 4]:** J² + P ≤ 7.8·10⁻¹⁰; |J| ∈ [3.8, 1659], |DJ| ≤ 1.3·10⁶; |N| median
  1.1·10⁻⁸, max 4.0·10⁻⁶; |N|/scale ≤ 3.2·10⁻¹⁰.
* **Filling at 0** (including t = 0 exactly and the node): J² + P ≤ 1.5·10⁻¹³; |J| ∈ [2.45, 23.9]; |N| median
  8.8·10⁻¹⁵; |N|/scale ≤ 1.4·10⁻¹⁴.
* **Filling at ∞:** J² + P ≤ 1.5·10⁻¹²; |J| ∈ [2.45, 7.7]; |N| median 1.5·10⁻¹⁴; |N|/scale ≤ 3.7·10⁻¹³.
* **Same points through two charts** (c0 and middle): q to 4.8·10⁻¹⁵, J to 7.7·10⁻¹⁴, DJ to 1.2·10⁻¹³.
* **Control:** J_ε = M J₀ M⁻¹, with M = I + ε sin(p₀)E₁₃, gives |N|/scale = 3.8·10⁻³ at ε = 0.1 and 3.3·10⁻⁵ at
  ε = 10⁻³. N is detected and is linear in ε.

**Caveat (important for any write-up).** N(F_*J₀) = 0 holds for *every* smooth F. At a point it depends only
on the 2-jet of F, and vanishes for every symmetric second derivative. So the N computation validates the
transport and the ambient code (and the N-control shows its sensitivity); it does **not** test F. The
substantive facts are that F is a global diffeomorphism (§12.5–12.6) and that J is the same tensor from
different charts (control above). This is why the run before the `linalg.solve` fix (§12.8), with a *wrong*
D²F in the middle chart, still gave |N|/scale ≈ 3·10⁻¹¹; the chart-agreement control is what exposed that
error (DJ differed by O(1)).

**F⁻¹.** Points of the round S⁴ × T² are inverted by an explicit chain followed by a Newton polish
(`tinverse.py`):
1. invert P_final (radially) and moebius (stereographic projection);
2. invert |T(r)| (1D);
3. invert R and h_B;
4. recover the fibre coordinate (explicit in the middle, 1D Newton at the cusps);
5. recover the S¹ phase from Z;
6. recover the T² part from φ;
7. polish with p ← p + DF⁺(q − F(p)) in the output chart.

On 3000 uniform points the chain alone is exact except in the thin linear-blend shells (max 1.9·10⁻³ there), and
the polished residual is ≤ 8.4·10⁻¹⁴ (median 2.5·10⁻¹⁶). At these points J² + P ≤ 1.1·10⁻¹⁰ and
|N|/scale ≤ 3.9·10⁻¹². So J_round can be evaluated anywhere.

**Geometry of F (quality, not correctness).**
* F is far from volume-uniform: of 3000 uniform points of S⁴ × T², 2659 come from the filling at 0
  (|t| < 0.2), 268 from 0.2 ≤ |t| ≤ 5, and 73 from the filling at ∞.
* J_round is strongly non-orthogonal in the middle region: |J| (Frobenius; √6 for an orthogonal J) has median
  43–100 for every |t| ∈ [0.25, 4], max about 1.7·10³; in the fillings it is 2.45–24.
* The cause is structural. F maps the lattice coordinates x₁ (S¹-orbit, circumference 2π|Z|) and x₂, x₃ (unit
  circles of T²) to round circles, while J_Y pairs lattice directions of very different lengths: near t = −1,
  i·e₁ = 2π·e₂, so |J| ~ (2π)²/|Z|.
* This matters only if J_round is to be used as a supervised target. Possible improvements: rescaled torus
  coordinates, or a dilation of the moebius chart (it only rescales the linear models).

### 12.8 Pitfalls met in this round (do not repeat)

* **`torch.linalg.solve` has wrong nested forward-mode derivatives (torch 2.10).** First derivatives are fine;
  ∂²/∂t² of the lattice coordinates was off by O(1). `lattice_coords` now uses Cramer's rule on the 3 × 3 real
  system (e₁ is real). After the fix, D²F from nested AD matches finite differences of DF in all charts.
  Isolated tests on torch 2.10: `linalg.solve` (13% error) and `linalg.det` (100%) are wrong in fwd-over-fwd
  and rev-over-fwd; `inv`, `lstsq`, `pinv`, `solve_triangular`, `cholesky` and `eigh` are correct in all modes.
  Never take a second derivative through `solve` or `det` when the inner derivative is forward-mode. A
  repository search found no current path doing so.
* **`z.abs() ** 2` has a wrong second AD derivative at z = 0** (0 instead of 2): it differentiates through the
  non-smooth |z|. First derivatives are right. This hit exactly the central fibre (x = 0 or y = 0) and the nodes,
  where D²F in the filling charts was off by O(1). All squared moduli now use `tcore.abs2(z) = Re² + Im²`.
  `nijenhuis_round.check_D2F` (nested AD vs finite differences of DF, all (k, l), all charts, including t = 0
  and the node) is the guard: it now gives ≤ 7·10⁻⁸ (middle) and ≤ 3·10⁻⁶ (fillings, FD-limited).
  Like the solve bug, this was invisible to N.
* **Pseudo-inverse.** Normal equations square cond(A), which reaches ~3·10⁴ in (z, w, t) coordinates. Use the
  SVD pinv and the exact normal projector I − P(q).
* **A log-ratio of two sections with a common zero** must never be evaluated where that zero is not a zero of
  the result (the old blend).
* **The closed-form arc |x| = |y|** is exact only for r ≤ 0.366.
* **A ψ̂-sum of local models y + c·x̄** forces its zero arc to arg t ∈ (π/2, 3π/2).
* **X is a bad parameter along the arc** (the middle of W is compressed in X); use log|t|.
* **Torch masking idioms.**
  * Evaluate branches only on their index sets (`_scatter`), and keep unused inputs finite (`where_safe`).
  * Round to the nearest node with `.detach()`.
  * `bump_step` sets exp(−1/u) to 0 below u = 10⁻³, which is identical to the numpy version.
* **Run `check_tcore` after any change to `tcore.py`**: the numpy modules are the reference.

### 12.9 Status, open points, follow-ups

* **Done:** the explicit F: Y → S⁴ × T², its verification, J on the round S⁴ × T², N with the ambient code,
  and F⁻¹.
* **Nothing is committed.** `docs/s4t2_explicit_diffeo.md` and `extra/s4t2_explicit/` (including
  `arc_fit.json` and `psi_fit.json`) are untracked. Commit only when the user asks; no HPC material is involved.
* **Next (chosen 2026-10-07): the symbolic programme of §13.** Other possible steps:
  * reduce the metric distortion of F (§12.7);
  * use J_round as a supervised target or a certification benchmark on S⁴ × T² (it is T³-invariant, i.e.
    determined by data on S⁴/S¹ = S³);
  * a write-up as a computational companion to the paper. The ζ ↦ λζ action, the third proof of W ≅ S⁴ via
    Fintushel's classification, and the explicit F are manuscript flags for the authors, not yet acted on.

### 12.10 Final reruns (2026-10-07, after both AD fixes)

All PASS, in this order:
* `check_tcore`, `check_tlift`, and the unchanged numpy checks `check_step1`, `check_step2`, `check_glue`,
  `check_hopf`;
* `check_full`: as in §12.6;
* `check_lift`: as in §12.5;
* `check_degree` at three grids: 0.999949, 0.999999, 1.000000;
* `nijenhuis_round 2000`: as in §12.7;
* `check_inverse 3000`: as in §12.7.

Run times on the laptop: check_lift 2 min, check_degree on the finest grid 2.5 min, nijenhuis_round 2000 1.7 min,
check_inverse 3000 under 1 min. The fits are deterministic: `fit_arc.py` takes 6 s, `fit_psi.py` 7 s (seeded).
Rerunning them reproduces `arc_fit.json` and `psi_fit.json`.

**Resume checklist.** Read §12.0–12.1. Run `check_tcore`, `check_tlift`, `check_full` (seconds), then
`check_lift` and `nijenhuis_round` if needed. Any change to E, Λ or the gauge requires rerunning `fit_arc.py`,
then `fit_psi.py`, then all checks. Never edit `tcore.py` without `check_tcore.py`.

### 12.11 J_round as a black-box field on the round S⁴ × T² (`jround.py`, `check_jround.py`, 2026-10-07)

`J_round(q)` is intended to evaluate the complex structure at points q of S⁴ × S¹ × S¹ ⊂ ℝ⁹ (unit S⁴, unit circles):
* p = F⁻¹(q), in a holomorphic chart of Y chosen automatically;
* J = DF·J₀·DF⁺, and DJ from D²F.

**Implementation review (2026-10-07).** The original polar inverse returned no chart entries at either exact
pole, so `J_round` silently returned J = 0 there. `tinverse.linear_ball_inverse` now implements the closed
inverse of §12.4 on both linear caps, including the poles (programme C8); `Finv` checks chart coverage.
The added 96-point inverse and 16-point J/DJ checks include poles, axes and cap boundaries. However, inverse
evaluation on the central fibres outside the caps is not yet reliable: the polar chain can reach a singular
Newton input (reproduction in ledger 2026-10-07-EH-4; campaign T14). Thus the sampled results below do not
establish an evaluator at every point. The forward formulas and their sampled checks are unaffected by this issue.

Cost on the laptop: about 20 ms per point with DJ. `check_jround.py` tests it the way an independent user would
(300 uniform points, plus 100 points on each chart switch |t| = 0.2 and |t| = 5; all PASS, 14 min):

| test | uniform | switch \|t\| = 0.2 | switch \|t\| = 5 |
|---|---|---|---|
| inverse residual \|F(F⁻¹q) − q\| | 3.5·10⁻¹⁴ | 3.8·10⁻¹⁵ | 1.7·10⁻¹⁴ |
| J² + P | 5.7·10⁻¹² | 5.4·10⁻¹³ | 9.9·10⁻¹² |
| black-box DJ (central differences of q ↦ J_round(q) along a tangent frame, through F⁻¹) vs AD DJ; h = 10⁻⁵ / 10⁻⁶ | 5.0·10⁻⁴ / 5.0·10⁻⁶ (h², truncation) | 1.7·10⁻⁷ / 3.9·10⁻⁹ | 6.2·10⁻⁷ / 1.2·10⁻⁸ |
| \|N\|/scale from the finite-difference DJ only (no D²F), h = 10⁻⁶ | 2.4·10⁻⁶ | 1.1·10⁻⁸ | 1.1·10⁻⁸ |
| \|N\|/scale from the exact DJ | 3.4·10⁻¹² | 7.9·10⁻¹⁴ | 1.2·10⁻¹¹ |

* **No jump at the chart switches.** At the switch points the finite-difference stencils straddle the switch, and
  the finite-difference derivative still matches AD to 10⁻⁸.
* **Control.** J_ε = R·J_round·R⁺, with R = P(1 + ε sin(3q₀)E)P, through the same black box with finite
  differences only:
  * ε = 10⁻¹: |N|/scale median 0.11 (max 0.60);
  * ε = 10⁻³: median 1.1·10⁻³ (max 2.3·10⁻²).

  Linear in ε, and far above the level of J_round.
* **T³-invariance.** J_round(g·q) = dg·J_round(q)·dgᵀ for g = (e^{iθ} on (Z₁, Z₂), rotations of both circles):
  2.9·10⁻¹³. So J_round is a T³-invariant complex structure for the standard isometric T³-action on the round
  S⁴ × T², and is determined by its values on a 3-dimensional slice.
* **∞ chart vs middle chart** at the hand-over |t| ∈ [4, 5]: q 7.6·10⁻¹⁵, J 9.8·10⁻¹³, DJ 5.8·10⁻¹⁰. For
  |t| ∈ [5, 20] this degrades to J 1.4·10⁻¹¹, DJ 1.3·10⁻⁸. That band is outside the middle chart's operating
  range: there w ~ t²w′ and cond DF reaches 6·10⁴ (≤ 2·10² in the ∞ chart). It is round-off amplification;
  J_round never uses the middle chart there.

**What "exactly, explicitly, N vanishes numerically" means (the precise statement):**
1. **Exactness.** J_round := F_*J_Y is an *exact* complex structure on the standard round S⁴ × T², not an
   approximation, provided that:
   * J_Y is integrable (the paper: Y is a complex manifold built from holomorphic charts and holomorphic
     gluings), and
   * F is a diffeomorphism (§12.5–12.6).

   It is the complex structure of [HSSS26], transported. The novelty is having it on the standard model, not a
   new structure.
2. **Explicitness.**
   * F is a composition of closed-form maps: the paper's lattice, Φ, Θ, ψ; theta-type heights; finite sums of
     local models; Hopf, Möbius and radial maps; explicit gauges; and φ.
   * The only non-closed-form ingredients are two fixed smooth functions of one variable given by stored
     Chebyshev coefficients (arc model, ψ). They only select *which* diffeomorphism F is, and the diffeomorphism
     checks validate the stored choice.
   * The height μ̃ is a convergent series, evaluated to below double precision.
   * F⁻¹ is not closed-form, so J_round(q) at a *given* q needs a Newton solve. Equivalently, J_round is given
     explicitly in the parametrisation q = F(p): J_round(F(p)) = DF(p)·J₀·DF(p)⁻¹.
3. **The numerical N-check.** N vanishes to round-off, from exact derivatives and from black-box finite
   differences. Because J_round is a pushforward of J₀ in holomorphic charts, N = 0 is guaranteed once the charts
   are holomorphic and F is a local diffeomorphism. The check therefore confirms the implementation (transport,
   derivatives, chart consistency, the ambient code); it is not an independent proof of the paper's
   integrability, which our formulas take from the paper. What is genuinely established beyond the paper is the
   explicit identification Y ≅ S⁴ × T², i.e. that J_round lives on the standard S⁴ × T².
4. **Rigour.**
   * The diffeomorphism property rests on dense sampling, adversarial searches, the degree integral (1.000000
     with det > 0 at all 18.4M nodes), and the argument "local diffeomorphism + compact + simply connected
     target ⇒ diffeomorphism".
   * This is strong numerical evidence, not a computer-assisted proof. A formal certificate would need interval
     arithmetic or analytic lower bounds for det DF.
   * Smoothness (C^∞) holds by construction; C¹ and the D²F consistency are confirmed numerically.

### 12.12 Can this be done symbolically, as a rigorous proof? (analysis, 2026-10-07)

**Question (user).** Can the construction be done symbolically, so that, together with an explicit computation
showing N ≡ 0, it gives a rigorous existence proof of the complex structure?

**1. A symbolic N computation would prove nothing about F.** For J = F_*J_Y, N_J(F_*X, F_*Y) = F_*N_{J_Y}(X, Y)
(naturality), and N_{J_Y} = 0 because Y is glued from holomorphic charts by holomorphic maps. So N ≡ 0 holds for
*every* local diffeomorphism F, and at a point it only sees the 2-jet of F. `check_nblind.py` demonstrates this:
* replacing our F's 1- and 2-jets by those of a random other local diffeomorphism (random A·M, random symmetric
  second derivatives satisfying the sphere constraints) gives |N|/scale ≤ 8.7·10⁻¹³ with the same code;
* only non-symmetric "second derivatives", which are not the 2-jet of any map, give N ≠ 0 (0.18).

A symbolic simplification of our N to 0 would just re-derive the chain rule, with enormous expression swell (theta
series, exp(−1/u) bumps, Chebyshev polynomials of degree 120, logarithms).

**2. Integrability of J_Y is identities, and is already symbolic.** It is not an N computation:
* the gluing maps are holomorphic;
* the cocycle and compatibility conditions hold;
* the fillings are complex surfaces.

The paper's `manuscript_draft/code/verify_construction.py` (agentic_ledger branch) checks these with sympy, 15
checks. Its only inequality, det_R(e₁, …, e₄) > 0, holds identically.

**3. Where rigour is needed: F is a diffeomorphism.** This is a global statement, made of inequalities over a
continuum plus topology. Symbolic algebra alone cannot prove such inequalities. A rigorous version is a
**computer-assisted proof** with three layers:
* **(a) Symbolic identities (sympy).**
  * Equivariance; well-definedness (lattice, branch, G₂–G₄, ψ̂); chart consistency.
  * Hopf(Z) = v, |Z|² + X² = 1; exact linear models at p₀, p_∞; the T²-part φ.
  * Already proven symbolically in `symbolic_identities.py`:
    * dμ̃/dρ = −2·Var_w(n), so the cusp height is strictly monotone along the fibres;
    * DΛ(b₀) = 2π·Q₀ exactly (det 8π³).
* **(b) Certified inequalities (interval arithmetic)**, a short list, with the numerical margins found:

  | | inequality | numerical margin |
  |---|---|---|
  | I1 | the Hopf fibre coordinate is strictly monotone in the cusp blend regions (the node correction and the λ-blend to x₄^slice; μ̃ itself by (a)) | positive throughout (`check_hopf` item 4) |
  | I2 | det ≠ 0 for Λ_mod in the two linear-blend shells (elementary closed-form maps there) | — |
  | I3 | E ≠ 0 outside the tube {θ_v < 15°} | the arc lies within 0.056° of M₁ |
  | I4 | gauge domain on {θ_v ≤ 45°}: \|arg t\| > 90° and \|μ̃ mod 1\| < ½ | 134° and 0.32 |
  | I5 | \|Arg g\| < π on the shell | 1.62 |

  Rigorous tail bounds are also needed for the truncated μ̃ series. The domain is compact: B, 3-dimensional, is
  covered by three compact charts.
* **(c) Topology by hand.**
  * Λ = P_final ∘ moebius ∘ R ∘ h_B is a composition of diffeomorphisms: |T(r)| is increasing by hand
    (|d log k/d log r| ≤ 0.77 < 1), R has an explicit inverse, moebius is classical, P_final needs ν increasing,
    and h_B needs I1. Λ_mod adds I2.
  * With I3–I5 the lift is a smooth equivariant bijection covering Λ_mod, and a local diffeomorphism because DΛ is
    invertible.
  * F is a T²-equivariant map covering Ψ_W, hence a diffeomorphism; and J = F_*J_Y is integrable by naturality.
  * The global Jacobian sampling and the degree integral become unnecessary in a proof.
* **Cost estimate.**
  * In pure-Python `mpmath.iv`, one certified enclosure of the cusp height on one box takes 0.9 ms. A full F is
    about 100× that, so 10⁶ boxes cost about 25 CPU-hours (parallelisable on Viper).
  * Arb (python-flint, not installed) or a vectorised outward-rounded numpy implementation would be 10²–10³×
    faster.
  * Steep regions (the fixed-point shells, the gauge glue) need fine boxes.
  * Estimate: 1–3 weeks of focused work for a complete certificate. The fitted Chebyshev coefficients are no
    obstacle: they are part of the definition of the certified F.

**4. Existence is already rigorous.**
* The paper's Main Theorem (unconditional) proves it: Y is a complex manifold, and Y ≅ S⁴ × T² via H²(W) = 0 plus
  Wall/s-cobordism or Iwase.
* The S¹-action observation (§1–2) gives a short third rigorous proof of W ≅ S⁴ with **no inequalities**: B = W/S¹
  is two solid tori glued with a clutching of degree 1, hence S³ (topological). The action is semi-free with 2
  fixed points (χ = 2), and W is simply connected, so W ≅ S⁴ by Fintushel's classification.
* A computer-assisted certificate of F would therefore add constructiveness, not existence: an explicit, certified
  diffeomorphism, and hence a certified explicit complex structure on the standard round S⁴ × T².

**5. Closed form on S⁴ × T².** J is closed-form in the parametrisation q = F(p): J(F(p)) = DF(p)·J₀·DF(p)⁻¹. As a
function of q it needs F⁻¹, which has no closed form (it inverts μ̃, |T| and the blends). For a proof that is
irrelevant; pointwise, F⁻¹ can be certified by interval Newton (Krawczyk).

### 12.13 A purely symbolic computer-assisted proof (no interval arithmetic)? (analysis, 2026-10-07)

**The proposal (user).** Use F to build a smooth (1,1)-tensor J on all of S⁴ × T². Then verify symbolically that
J² = −1 and N_J ≡ 0, bypassing any proof of properties of F.

**Opinion: not feasible as stated.** The obstruction is logical, not computational.

1. **"A smooth J on the whole S⁴ × T² built from F" is equivalent to the diffeomorphism properties.**
   * Pushforward, J(q) = DF(p)·J₀·DF(p)⁻¹: every q needs exactly one preimage p, with DF(p) invertible. That is,
     F must be a diffeomorphism. A covering would do only with deck-invariance.
   * Pullback through an explicit G: S⁴ × T² → Y, J = G*J_Y: this is the minimal version. It needs only
     well-definedness of G (identities) and **det DG ≠ 0 at every point** (no injectivity).

   Either way an "at every point" non-degeneracy remains. For our transcendental formulas that is an inequality
   over a continuum, not an identity.
2. **Given such an assignment, J² = −1 and N ≡ 0 are automatic** (J₀² = −1 and naturality of N). A symbolic check
   cannot distinguish a good map from a bad one: `toy_symbolic_blindness.py` (sympy) takes
   F(x, y, u, v) = (x², y + xu, u, v) on ℝ⁴.
   * F is two-to-one, and degenerate on {x = 0} (det DF = 2x).
   * Nevertheless J² + 1 ≡ 0 and all 64 components of N_J ≡ 0 identically.
   * Yet J is not a tensor field on the target: the two preimages of one point give different J, and
     J[1, 0] = (u² + 1)/(2x) is unbounded.

   The numerical analogue is `check_nblind.py` (§12.12). In addition, brute-force symbolic zero-testing of
   expressions with exp, log, abs and sin is undecidable in general (Richardson 1968). So even the automatic part
   should rest on the naturality lemma, not on CAS simplification.
3. **The only ways to a purely symbolic certificate** make the non-degeneracy an identity or manifest. Each is a
   research redesign of uncertain outcome:
   * **(i) Explicit inverse.** An explicit G with F∘G = id verified symbolically gives DF·DG = I, hence a
     diffeomorphism, by the chain rule. Obstacle: F⁻¹ needs transcendental inversions:
     * the theta-type height μ̃, which resolves the central fibres;
     * |T(r)| and ν;
     * the λ- and linear blends.

     The closed-form-invertible fibre coordinate x₄^slice = log|ζ|/P collapses the central fibre onto the node
     (it is only continuous there). A ψ̂-equivariant fibre coordinate smooth through the central fibre is an
     Abel-type function for ψ̂, built by C^∞ blends. Blends destroy both closed-form invertibility and manifest
     monotonicity.
   * **(ii) Manifest non-degeneracy.** Every blend designed so that its Jacobian is manifestly nonzero (e.g. sums
     of manifestly positive terms). Our F has five non-manifest conditions (I1–I5 of §12.12). The global ones (E ≠ 0
     off the tube, gluing phase < π) are the least amenable.
   * **(iii) A J not obtained by transport**, given by globally manifestly smooth formulas, e.g. polynomial in the
     ambient coordinates. Then J² = −P and N = 0 are polynomial identities modulo the ideal of S⁴ × T², decidable by
     Gröbner bases. This is genuinely symbolic, but finding such a J is open. The paper's J is intrinsically
     transcendental (log t, elliptic periods), and a polynomial representative is implausible.
4. **Recommendation.**
   * For rigour of existence: the paper, or the short S¹ + Fintushel proof (no inequalities).
   * For a rigorous *explicit* J on the round S⁴ × T²: validated numerics for I1–I5 (§12.12).
   * Symbolics remain the right tool for all the identity-type steps, which are numerous.

## 13. CURRENT PROGRAMME (agreed 2026-10-07): a symbolic description of F and of its inverse G

**Goal (user, 2026-10-07):** "a completely symbolic description, explicit and sound, of F and its inverse G". It
should give a purely symbolic computer-assisted proof (no interval arithmetic) that the round S⁴ × T² carries
the complex structure F_*J_Y = G*J_Y.

> **Progress record (2026-10-07, after compaction): the programme is being carried out in
> `docs/s4t2_symbolic_programme.md`** (block-by-block symbolic specification of F and G, sympy scripts in
> `extra/s4t2_explicit/sym/`, exact facts found on the way, status of I1–I5). Read that file first; this section
> remains the record of the scheme and of the audit as agreed before the work started. Corrections established
> there: the smooth-step lemma s′ ≤ 2 is proven (L1); ν′ > 0 is manifest (L3); |T| increasing is proven (L4);
> the λ-blend needs no `wrap` (so the wrap-bound lemma of 13.3.3 is unnecessary); φ(t) = 2|t|²/(1 − |t|²) exactly
> at both cusps (N1).

### 13.1 The two clarifications, and the agreed proof scheme

**Q1 (user):** does building J(q) = DF(p)·J₀·DF(p)⁻¹ require a missing piece, the inverse of F?

*Answer:* the formula needs p = F⁻¹(q).
* We can *evaluate* F⁻¹ numerically (`tinverse.py`, about 10⁻¹³).
* The local identities (J² = −1, N = 0) can be written in terms of p without any inverse, as in the toy.
* What is missing is a *proof* that F⁻¹ exists and is smooth, i.e. that F is a diffeomorphism: every q has exactly
  one preimage, and DF(p) is invertible everywhere.

**Q2 (user):** suppose we have F and its inverse G, build J blindly, and check that (1) the assignment is smooth,
well defined and respects the global structure, (2) J² = −1, (3) N_J = 0. Is that a proof?

*Answer:* yes, by Newlander–Nirenberg, with two refinements:
* Checks (2) and (3) can never fail. J is the pullback of J_Y by a local diffeomorphism, so they hold by naturality;
  they would only catch coding errors.
* Check (1) is the entire content, and with an explicit G it becomes symbolic, as the scheme below shows.

**Proof scheme (target).** Use finitely many pieces covering S⁴ × T², with an explicit formula on each.
* **(S1) G: S⁴ × T² → Y is manifestly smooth.**
  * Every denominator is manifestly ≠ 0, every log/sqrt argument manifestly > 0, and no branch cut is crossed. The
    domains manifestly cover S⁴ × T².
  * On overlaps the pieces give the same point of Y. These are identities through Y's chart transitions (Φ, Θ,
    G₂, G₃, G₄, the toric chart changes).
* **(S2) F is explicit and manifestly smooth** near G(S⁴ × T²), and **F∘G = id_{S⁴×T²} holds as a symbolic
  identity.**
* **Consequences** (one-line arguments, no computation):
  * DF(G(q))·DG(q) = I, so DG is invertible everywhere.
  * J := G*J_Y, which in each chart is DG⁻¹·J₀·DG = (DF∘G)·J₀·DG, is a smooth (1,1)-tensor field with J² = −1 and
    N_J = G*N_{J_Y} = 0. By Newlander–Nirenberg S⁴ × T² carries a complex structure.
  * G is injective (it has a left inverse) and a local diffeomorphism from a compact manifold to a connected one,
    hence a diffeomorphism, and G∘F = id. This is not needed for existence.
* **Relaxation, "certified implicit" pieces.** A piece may be defined implicitly by f(x, y) = 0 for y, provided:
  * ∂_y f is manifestly invertible, and
  * a solution manifestly exists (e.g. a strictly monotone function with full range).

  Smoothness then follows from the implicit function theorem. Example: the inverse of the cusp height μ̃ along the
  fibres, since dμ̃/dρ = −2·Var_w(n) < 0 manifestly (all weights > 0, at least two values of n) and μ̃ → ∓∞ as
  ρ → ±∞.
* **Integrability of Y** enters only through the paper's construction (holomorphic charts and holomorphic gluings,
  identities in `verify_construction.py`).
* **Not needed:** global Jacobian sampling, degree integrals, symbolic N computations, interval arithmetic.

### 13.2 The current construction against the scheme (F-direction, block by block; audit of 2026-10-07)

| block (code) | non-degeneracy certificate | inverse | status |
|---|---|---|---|
| lattice coordinates, Φ, Θ, G₂–G₄, toric chart changes (`tcore`, `tfull`) | holomorphic, explicit | explicit | ✓ (identities; the paper) |
| T²-trivialisation φ (`tfull.phi_mid`, φ₀ = C − H c̃) | ∂φ/∂(x₂, x₃) = I | explicit given the W-point: (x₂, x₃) = φ + corrections | ✓ |
| S¹ phase of the lift | identity | explicit given the orbit point | ✓ |
| cusp height μ̃ (theta-type series; `mu_torus`, `mu_chart`) | dμ̃/dρ = −2 Var < 0, **manifest** (`symbolic_identities.py`) | certified implicit | ✓ |
| node correction μ̃_new = μ̃ + c(rad)(m₀ − μ̃) (`node_correct`) | **not manifest** (the c′ term), I1 | — | **redesign** |
| λ-blend a₁ = μ̃_new + λ(\|t\|)·wrap(x₄^slice − ½ − μ̃_new) (`a1_c0`, `a2_ci`) | equals the convex combination (1 − λ)μ̃_new + λ(x₄^slice − ½), monotone since the weight depends on \|t\| only, **provided** \|x₄^slice − ½ − μ̃_new\| < ½ on 0.3 ≤ \|t\| ≤ 0.5. Otherwise `wrap` jumps. Numerically ≤ 0.375 (`check_step2`); a theta-function lemma is needed | certified implicit | lemma, once μ̃_new is fixed |
| base profile \|T(r)\| = r k(r) (`k_of`, `T_of`) | d log\|T\|/d log r ≥ 1 − 0.77 > 0, needs the lemma max s′ = 2 for the smooth step s | certified implicit | lemma |
| Hopf identification h_B (`hopf_from_a1/a2`) | from the two rows above, plus the clutching identity | explicit given r | ✓ |
| straightening R(u) = (φ̄₁(r)u₁, φ̄₂(r)u₂) (`arc_phases`; Chebyshev polynomials with fixed coefficients) | det = 1 manifestly: in (\|u₁\|, arg u₁, arg u₂) it is triangular with unit diagonal (a twist) | explicit given r | ✓ |
| moebius (`tcore.moebius`) | conformal | explicit (stereographic) | ✓ |
| P_final (`nu_fn`) | ν strictly increasing in \|τ\|: one-variable lemma for the blend on [½, 1] | certified implicit (ν⁻¹) | lemma or redesign |
| linear blends Λ_mod near b₀, b_∞ (`lin_blend`) | **not manifest**, I2 | Newton only | **redesign** |
| E (cusps: χ-weighted ψ̂-sums; middle: linear blend at radii 0.7 and 1/0.7) | E ≠ 0 off the tube, I3: the cusp part has a near-manifest "aligned phases" argument (§7); the middle blend is **not manifest** | — | **redesign or prove** |
| gauge Ĝ_A with ψ(\|t\|) (`GA_hat_W`, `psi_of_r`, Chebyshev fit) | domain conditions, I4 | — | **redesign** |
| shell gluing e^{ic Arg g} (`assemble`) | \|Arg g\| < π, I5: **not manifest** | explicit given the orbit point | **redesign** |

### 13.3 Hard points and design ideas (to be developed)

1. **Fixed points (I1, I2).** The current implementation ensures smoothness at p₀ by making Λ
   exactly linear in the Hopf coordinates v_q near b₀. Exact linearity is sufficient, not necessary:
   programme §5.2 and `sym/lemmas_1d.py` L9 give a smooth nonlinear radial example.
   The plain theta height is not even C² there: along the axis
   μ̃ = v₃ − v₃|v₃| + O(v₃³) (verified 2026-10-07: μ̃ − v₃ = −X² on D₁, +Y² on D₀). Hence the node correction and
   the linear blend. Options:
   * (a) define Ψ near p₀ directly upstairs, through the linear model √(2π)·U₀·(x̄, y), and design the transition
     to the orbit-space construction with manifest non-degeneracy;
   * (b) choose h_B, R, moebius and P_final so that their composition is exactly conformal-linear near the fixed
     points, which removes the blend;
   * (c) a node correction whose monotonicity is manifest.
2. **The lift (I3–I5), the least amenable.** Options:
   * (a) make the transition between the two lifts an explicit e^{iθ} with θ an explicit real function, so that no
     branch condition arises;
   * (b) adapt Λ to E so that one lift suffices (exact compatibility along the zero arc);
   * (c) use symmetry. The zero arc of E lies exactly on arg t = π (to 4·10⁻¹⁶) wherever E is a pure cusp formula,
     i.e. |t| ≤ 0.8 and |t| ≥ 1.25; it deviates only in the linear blend. This suggests a reflection symmetry, to be
     identified. A symmetric middle blend might keep the arc exactly on arg t = π, and make R explicit without a
     fit;
   * (d) prove E ≠ 0 off the arc in the cusps symbolically by the aligned-phases argument of §7: on slab overlaps
     consecutive terms are positive multiples of the same phase.
3. **One-variable lemmas.**
   * The smooth step s(u) = 1/(1 + e^{1/u − 1/(1−u)}) has s′(½) = 2 (sympy). The global maximum is 2 (numerically,
     at u = ½); this is still to be proven.
   * ν′ > 0 on [½, 1] for P_final.
   * The λ-blend bound |x₄^slice − ½ − μ̃_new| < ½ on 0.3 ≤ |t| ≤ 0.5 (numerically ≤ 0.375). It is a two-variable
     estimate on a theta-type mean, e.g. via Poisson summation; or redesign the blend so that no `wrap` is needed.
   * Alternatively, redesign the profiles with manifest positivity: e.g. densities bounded by exp(1 − 1/(1 − u²)) ≤ 1,
     which is manifest.
   * Useful: s has a closed-form inverse. s(u) = σ gives L u² + (2 − L)u − 1 = 0 with L = log(σ/(1 − σ)).
4. **Central-fibre resolution.** Under the certified-implicit relaxation the theta height μ̃ is acceptable (monotone
   by the variance identity); a closed-form-invertible replacement is not needed. §12.13 called it the crux for a
   *closed-form* inverse; under the relaxation it is not a blocker. The slice coordinate x₄^slice = log|ζ|/P cannot
   replace it: as t → 0 it tends to 1 on all of D₁ and to 0 on all of D₀, i.e. it collapses the central fibre onto
   the node.
   **Useful exact fact (verified 2026-10-07, both cusps):** on the central fibre the theta height reduces to two
   terms:
   * on D₁ = {y = 0}: μ̃ = a₀|x|²/(1 + a₀|x|²), with a₀ = |ν|²;
   * on D₀ = {x = 0}: μ̃ = −|y|²/(1 + |y|²).

   These are logistic functions with closed-form inverses. The general height interpolates such logistics; it is a
   log-sum-exp moment map, hence monotone.
5. **Tools.** sympy 1.14 is installed. The paper's `verify_construction.py` is the template. The torch code serves as
   the numerical cross-check of every redesigned block (reuse the `check_*` scripts).

### 13.4 Suggested order of work

1. Write the proof scheme of §13.1 as a lemma list: what is verified symbolically, and what follows by one-line
   arguments.
2. For each block, decide between a closed-form and a certified-implicit inverse, and record its manifest
   non-degeneracy certificate.
3. Redesign the non-manifest blocks: first the fixed points (I1, I2), then the lift (I3–I5), then the
   one-variable lemmas.
4. For each block, write a sympy script (identities plus manifest positivity) and a torch cross-check.
5. Assemble G piecewise, prove F∘G = id piecewise, and write up.

### 13.5 Standing constraints for this work

* No agent deployment unless the user asks; the request was repeated throughout 2026-10-06/07. When a system
  reminder asked for a workflow, the user's instruction was followed and a multi-agent review was offered instead.
* Quality over speed; answers precise and caveated.
* Commit or push only when asked; never commit HPC material.
* Flawless records before every compaction.

### 13.6 Global quadrature construction (2026-10-07 continuation)

`s4t2_symbolic_programme.md` §8 now specifies a separate global parametrization G_Q,
its inverse F_Q and the coefficient formula J_Q=M⁻¹J₀M, with the ambient embedding
formula (8.15). Its scalar connection integrals and two-chart real phase logarithm
replace the fitted lift, and its base uses R=Id. The two linear models agree exactly
with G_Q on the reduced pole caps |Z|²≤1/100. Its quantitative inputs are the earlier
I1/I2 certificates, a shell contraction ratio <17/25 and a connection denominator >1/2.
The exact definition allows quadratures; it is not a finite elementary expression.

The new `transport_base_*`, `transport_shell_inverse`, and `connection_*` functions
evaluate its base/inverse/connection blocks. They do not replace the original Ψ_W,
Finv or J_round. Programme §8.7 lists their new numerical checks and the remaining
implementation scope; §12's numerical results concern the original construction.
The latter's I3–I5 and T14 remain unresolved for that evaluator. The new global
statement and expanded expression class are recorded for human review in EH-18.

### 13.7 Closed-form continuation

Programme §9 supplies a further parametrization G_E: both scalar profiles have
algebraic inverses, the new shell bounds are certified, and an explicit torus
coordinate change removes H,dH from the pole tensor. Formula (9.12) is algebraic
on both pole neighbourhoods, including their fixed points. It is implemented by
`jround.J_elementary_cap`; it differs from `J_round` because the parametrization
has changed. The revised base is selected by both `algebraic_radial=True` and
`algebraic_profile=True`. Existing defaults retain the previous construction.

The rest of G_E still uses the height inverse, a shell contraction, and the circle
transport/phase quadratures. Programme §9.4 proves a restricted obstruction to
algebraic equivariant heights, not impossibility of a global elementary J. It also
explains why finite theta truncation, a universal endpoint potential, and unrestricted
principal Arg do not close the remaining gaps. No full global closed-form result
or complete numerical evaluator of G_E/J_E is claimed.

### 13.8 Literature-led scalar and bundle reductions (2026-10-08)

Programme §10 records the inspected primary literature and three independent designs.
The finite elementary height H_f has an exact finite forward formula, a cubic node
remainder constant 35 and corrected slope >14/25. The `finite_height=True` option
is implemented in the independent base, inverse, connection and torus-correction
blocks; its inverse still uses scalar solves, and it does not change the legacy defaults.

A different height H_e uses Jacobi amplitude and has an explicit inverse in standard
incomplete elliptic F, including a redesigned full cusp collar. Its exact axes are
elementary; programme §10.5 proves the central extension by an independent arctangent
series. `cusp.py` supplies scalar reference evaluators. The 120 high-precision
value/inverse checks, 72 collar checks and 36 seam-jet checks pass. This design is
not elementary, and has not been joined to a globally certified base and circle lift.
Applying the old node correction would reintroduce a scalar implicit inverse.

The specified algebraic target f_alg:S⁴→S² and its two critical values are also
implemented; a fibrewise identification with W is not supplied by those formulas.
Finally, finite phase-logarithm and rank-one projection operations are proved and
implemented with explicit derivative/nonvanishing hypotheses. Their tests include
phases spanning more than three turns, both stereographic charts and exact reversal
of projection chains. A concrete bundle embedding and certified derivative bounds
for this construction remain necessary before they can replace its transport.

The current result is therefore a set of verified reductions, not a full finite
global expression for the complex structure or its diffeomorphism/inverse pair.
Programme §10.6 gives the remaining matching and bundle-lift tasks. No general
closed-form impossibility theorem has been proved.

### 13.9 Finite elementary global tensor (2026-10-08 continuation)

Programme §11 now gives a different inverse-first parametrization G_cl and the
finite coefficient formula (11.19) for J_cl on the whole standard S⁴×T².
The elementary cusp inverse has a closed collar, matching deck jets and smooth
open central-fibre charts. A closed inverse pole patch has derivative defect
<1/1000. A concrete ten-component circle-bundle embedding supplies a rank-one
projector with proved derivative bound L=10¹¹⁴ on the specified slab.
Fixed finite projection and phase formulas use N=4L+1 and
K=8·32(N+1)L+1; exact phase matching joins them to both closed pole formulas.
The previously proved finite forward height restores T². This gives a smooth
bijection with invertible derivative; pulling back the holomorphic atlas proves
J²=−Id and N_J=0, with no implicit solve, quadrature or infinite series in J.
The derivative recipe includes the central fibres and uses local argument charts.
This new global mathematical statement requires human review.

The finite expression is enormous: its counts have 115 and 232 digits. It has
not been expanded or numerically executed at those counts. Scalar implementations
are in `cusp.py`, with an explicit work guard; `elementary_closed_parameters()`
returns the exact constants. `check_inverse.py --closed` and
`--closed-nijenhuis` check the new blocks, torus and branch identities, cap/slab
matching, chart-independent local tensors and finite-difference integrability.
Their free-region tests use smaller counts and do not certify those counts
globally. Programme §11.8 records the measured margins alongside the proof bounds.
The complete legacy safety net, including the 18.4-million-node degree check and
all 300-point J_round checks, was rerun successfully; all four symbolic suites
also pass. The manuscript and frozen modules remain untouched.

This closes the **tensor expression** in the finite elementary class. It does not
supply a finite elementary formula for F_cl=G_cl⁻¹: inverses of the transition
height and inverse pole patch may remain implicit, and are not evaluated in J_cl.
The stronger inverse-pair goal, a compact practical formula, and human review
therefore remain. The earlier designs and their open legacy I3–I5/T14 flags are
retained. No closed-form impossibility theorem is claimed.
