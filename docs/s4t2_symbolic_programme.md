# The symbolic programme: F: Y → S⁴ × T² and its inverse G in symbolic form (started 2026-10-07)

> **Status line (keep current).** Sections 0–6 written on 2026-10-07; all four scripts in `extra/s4t2_explicit/sym/`
> PASS. Every claim marked *(sympy)* is verified by one of them; numerical facts carry their source script.
> Companion of `docs/s4t2_explicit_diffeo.md` (§12 = the numerical construction, §13 = the agreed scheme).
> **Resume point: §6.3 (next steps), after reading §6.1–6.2.**

## 0. What is being built, in one paragraph

F = (Ψ_W ∘ π, φ): Y → S⁴ × T² is T³-equivariant. Its inverse G is therefore determined by (i) the inverse of the
orbit-space map Λ: B → B_S = S⁴/S¹ (B = Y/T³ ≅ S³), (ii) the S¹-phase, read off from Z against the lift at a
reference point of the orbit, (iii) the T²-phase, read off from φ. The programme writes every block of Λ, of the
lift and of φ as a closed-form map with a closed-form inverse, or as a *certified-implicit* block (one-variable,
strictly monotone with manifest derivative sign and manifest range), and verifies all identities with sympy. The
residual non-degeneracy conditions are listed in §5 with their exact logical status.

## 1. Reading guide

* §2 one-variable lemmas and elementary inverses (done, `sym/lemmas_1d.py`).
* §3 block-by-block symbolic specification of F and G (the core; written block by block).
* §4 the fixed points: exact facts and the Taylor-remainder certificates needed for I1, I2.
* §5 the residual global conditions (I3–I5) and the architecture options that reduce them.
* §6 status and resume checklist.

## 2. One-variable lemmas and elementary inverses (`sym/lemmas_1d.py`, all PASS)

* **L1 (smooth step).** s(u) = 1/(1 + e^{1/u − 1/(1−u)}), p = u(1−u), q = (1−2u)/p. Identities (sympy):
  s′ = s(1−s)(1/u² + 1/(1−u)²), s(1−s) = 1/(4cosh²(q/2)), 1/u² + 1/(1−u)² = (1−2p)/p², p²q² = 1 − 4p.
  Hence s′ = (1−2p)/(4p²cosh²(q/2)). Since cosh x ≥ 1 + x²/2 (all Taylor coefficients ≥ 0), cosh² ≥ 1 + x², so
  4p²cosh²(q/2) ≥ 4p² + p²q² = (1−2p)², and **s′ ≤ 1/(1−2p) ≤ 2**, with equality only at u = ½ (p = ¼).
  This is the lemma the base profile |T| needed (§13.3 of the notes); it is now proven.
* **L2 (inverse of s).** s(u) = σ ⇔ Lu² + (2−L)u − 1 = 0 with L = log(σ/(1−σ)); **u = (L − 2 + √(L²+4))/(2L)**,
  u → ½ as L → 0, u → 1 as L → ∞, u → 0 as L → −∞ (sympy; the root identity and the chain q = −L are symbolic,
  plus a 50-digit numerical confirmation). So every smooth-step blend with a weight that depends on one variable
  only has a closed-form inverse for the weight.
* **L3 (P_final profile).** ν(τ) = (1−b)τ² + bτ, b = s(2τ−1), on [½, 1]: ν′ = 2τ(1−b) + b + 2s′(2τ−1)τ(1−τ), three
  manifestly nonnegative terms, so ν′ ≥ min(2τ, 1) = 1. On [0, ½] ν′ = 2τ, on [1, ∞) ν′ = 1. **ν is strictly
  increasing on (0, ∞) with range (0, ∞), manifestly.** Certified-implicit inverse on [¼, 1] (closed form outside:
  √ and identity).
* **L4 (base profile).** |T|(r) = r·k(r), log k = (1−λ)log K₀ + λ log K_∞, λ = s((log r + 4)/8), K₀ = 4π,
  K_∞ = e²/(4π): d log|T|/d log r = 1 + s′·(1 − log 4π)/4 ≥ 1 − 2(log 4π − 1)/4 = (3 − log 4π)/2 ≈ 0.2345 > 0,
  using L1 and 4π < 16 < e³. **|T| is strictly increasing, range (0, ∞).** Certified-implicit inverse.
* **L5 (P_final).** On |τ| ≤ ½: |v| = e^{−τ²}, X = sign(τ)√(1 − e^{−τ²}); inverse τ = sign(X)√(−log|v|), smooth in X
  (−log(1−X²) = X²(1 + X²/2 + …)). On |τ| ≥ 1: |v| = e^{−|τ|}, inverse τ = −sign(X) log|v|. Between: τ =
  sign(X)·ν⁻¹(−log|v|) (L3). The direction: v/|v| = ρ_O(w/|w|), ρ_O a signed permutation, det −1.
* **L6 (moebius).** The stereographic projection from (−1, 0) ∈ S³ ⊂ ℝ⁴ shifted by (0, 1, 0), and the inverse
  u = ((1−|y|²)/(1+|y|²), 2y/(1+|y|²)) with y = w + (0,1,0), are mutually inverse (sympy, both compositions).
* **L7 (lifts).** For both lift formulas, Hopf(Z) = v and |Z|² = |v| identically (sympy in the parametrisation
  v = n(sinθcosφ, sinθsinφ, cosθ)). Hence (Z, X) ∈ S⁴ whenever (v, X) ∈ B_S.

## 3. Symbolic specification of F and G, block by block

Notation: a *block* is a map with its inverse. "closed" = both directions closed-form; "implicit(1D)" = the inverse
is a certified-implicit one-variable inversion (strictly monotone, manifest derivative sign, manifest range; the
implicit function theorem gives smoothness). Every identity marked (sympy) is in `sym/chart_identities.py` or
`sym/lemmas_1d.py`.

### 3.1 Regular part Y_mid (closed)

* **Coordinates.** (ζ, w, t) ∈ ℂ* × ℂ × ℂ*, z = −i log ζ; lattice e₁ = (2π, 0), e₂ = (−it, 1−t²), e₃ = (−t, i(1+t²)),
  e₄ = (−i log t, 2t) for a branch of log t. Lattice coordinates x ∈ ℝ⁴: (z, w) = Σ x_j e_j.
* **F-direction.** x = Cramer's rule on the real 3 × 3 system (Im z, Re w, Im w) for (x₂, x₃, x₄) (e₁ is real),
  then x₁ = (Re z − Σ_{j>1} x_j Re e_j^z)/2π. (sympy C1: reproduces x identically; the determinant of the system
  is the explicit polynomial printed by the script, nonzero on ℂ* for the branch in use — it is the paper's
  det_ℝ(e₁…e₄) > 0 up to the factor 2π.)
* **G-direction.** (z, w) = Σ x_j e_j(t), ζ = e^{iz}. Closed.
* **Group actions in these coordinates** (sympy C2, and the paper): S¹_a: x₁ += θ/2π; T²: (x₂, x₃) += c; g₄:
  x₄ += 1; branch change log t → log t + 2πi: x₁ → x₁ − x₄.
* **B_mid coordinates:** (t, x₄ mod 1), with x₁ mod 1 the S¹-fibre and (x₂, x₃) mod 1 the T²-fibre. W_mid =
  ℂ* × T²_{(x₁, x₄)}, B_mid = ℂ* × S¹_{x₄}. Both trivial: the reference point of an orbit is x₁ = x₂ = x₃ = 0.

### 3.2 Cusp fillings (closed, plus one implicit(1D) block: the height)

* **Toric chart U₀ at cusp 0:** ζ = y, t = xy; S¹_a: (x, y) ↦ (λ⁻¹x, λy). The cusp generator (after the
  T²-slice) is ψ̂(x, y) = (1/(my), mxy²) with m(t) = ν·e^{−(c̃₂α₂ + c̃₃α₃)}, c̃ = R(t)⁻¹ē₄ (real coefficients of ē₄
  in ē₂, ē₃); ψ̂ preserves t, multiplies ζ by m·t, and ψ̂⁻¹(x, y) = (mx²y, 1/(mx)) (sympy C3). Hence
  log|ζ| ∘ ψ̂ = log|ζ| + P(t), P = log|ν| + log|t| − φ(t), φ = Re(c̃₂α₂ + c̃₃α₃).
* **Chart ↔ regular part (closed, holomorphic):** cusp 0: Φ(ζ, w, t) = (z + π s(s−1), w + 2ts, t), s = log t/2πi;
  cusp ∞: Θ⁻¹: t = 1/s, w = w′t², z = z′ + iw/t. Both are the paper's maps; inverses explicit.
* **T²-slice of a filling point (x, y, w) (closed):** w = C₂ē₂ + C₃ē₃ (real_coeffs, a 2 × 2 real solve),
  (x_s, y_s) = (x e^{C·α}, y e^{−C·α}). Inverse: given (x_s, y_s, C), (x, y) = (x_s e^{−C·α}, y_s e^{C·α}), w = C·ē.
* **Height μ̃ (implicit(1D) along fibres; closed on the central fibre).**
  Torus form: μ̃ = Σ_n n w_n, w_n ∝ e^{−2nρ + n(n+1)P}, ρ = log|ζ|. Chart form in U₀: μ̃ = Σ_j jT_j / Σ_j T_j with
  T_j = e^{(j²+j)g} X^{j(j+1)/2} Y^{j(j−1)/2}, X = |x|², Y = |y|², g = log|ν| − φ(t); an entire series in (X, Y),
  super-exponentially convergent. Facts:
  * μ̃ ∘ ψ̂ = μ̃ + 1 (index shift n ↦ n+1 under ρ ↦ ρ + P: exact);
  * dμ̃/dρ = −2·Var_w(n) < 0 (sympy C5): strictly decreasing along every torus fibre, range ℝ (the extreme terms
    dominate as ρ → ±∞). **Inverse: certified-implicit.**
  * On the central fibre (sympy C4): on D₁ = {y = 0}: μ̃ = e^{2g}X/(1 + e^{2g}X), inverse X = μ̃e^{−2g}/(1−μ̃);
    on D₀ = {x = 0}: μ̃ = −Y/(1+Y), inverse Y = −μ̃/(1+μ̃); derivatives manifestly ≠ 0. At t = 0, e^{2g} = a₀ = |ν|².
  * Near the node the series starts T₁ = e^{2g}X, T₋₁ = Y, then T₂ = e^{6g}X³Y, T₋₂ = e^{2g}XY³ (sympy C4e): the
    remainder after the two-term truncation is O(rad⁴) in rad = a₀X + Y; see §4.
* **T²-trivialisation at the cusp (closed):** φ₀ = C − H·c̃(t), H = μ̃ at the slice point. G₄-invariant (C += c̃,
  H += 1; sympy C6), T²-equivariant. Inverse: C = φ₀ + H·c̃. Gluing to the middle: φ = (x₂, x₃) − (1 − χ(|t|))·D₀
  with D₀ = (H − x₄^slice)·c̃ + (ξ₂, ξ₃)(t), a smooth ℝ²-valued function of the W-point; inverse (x₂, x₃) = φ +
  (1 − χ)D₀, closed given the W-point.
* **The S¹-phase (closed):** given the W-point's reference p_ref (x₁ = 0 in the middle; y real positive in the
  chart) and Z ∈ ℂ², e^{iθ} = ⟨Z, Ψ_W(p_ref)⟩/|⟨…⟩|, and the point is e^{iθ}·p_ref (x₁ = θ/2π). The pairing is
  nonzero because Z and Ψ_W(p_ref) lie on the same S¹-orbit of S⁴ (Hopf(Z) = v = Hopf(Ψ_W(p_ref)), L7).
  chart) and Z ∈ ℂ², e^{iθ} = ⟨Z, Ψ_W(p_ref)⟩/|⟨…⟩|, and the point is e^{iθ}·p_ref (x₁ = θ/2π). The pairing is
  nonzero because Z and Ψ_W(p_ref) lie on the same S¹-orbit of S⁴ (Hopf(Z) = v = Hopf(Ψ_W(p_ref)), L7).

*Correction to 3.1:* the determinant of the lattice system is det = (|t|² + 1)(log|t|·(|t|² − 1) + 2|t|²)
(sympy C1d), branch-independent and **manifestly positive on ℂ***: log|t| and |t|² − 1 have the same sign.

### 3.3 The orbit-space map Λ = P_final ∘ moebius ∘ R ∘ h_B and its inverse

* **Fibre coordinate on B₀ = {|t| ≤ ½} (implicit(1D) along fibres).**
  a₁ = (1 − λ(|t|))·μ̃_new + λ(|t|)·(x₄^slice − ½), λ = s((|t| − 0.3)/0.2), x₄^slice = ρ/P(t).
  **Design change (symbolic version): no `wrap`.** Both μ̃_new and x₄^slice − ½ are globally defined functions of
  ρ on the fibre, each shifting by exactly +1 under ψ̂, so the convex combination is equivariant without any
  integer choice. The code's `wrap(x₄^slice − ½ − μ̃_new)` equals the plain difference wherever |…| < ½, which holds
  numerically on 0.3 ≤ |t| ≤ 0.5 (≤ 0.375, `check_step2`), so the two definitions agree on the current code's range;
  the wrap-bound lemma of notes §13.3 is **no longer needed**. Then
  da₁/dρ = (1 − λ)·dμ̃_new/dρ + λ/P(t) < 0 manifestly where μ̃_new = μ̃ (P < 0 on |t| ≤ ½: P = log|ν| + log|t| − φ(t) with
  φ = 2|t|²/(1 − |t|²) ≥ 0 exactly (§4.1), so every term is ≤ 0 and log|t| < 0), since λ depends on |t| only. The node correction μ̃_new ≠ μ̃ only for
  rad = a₀|x₀|² + |y₀|² < 0.05, i.e. |t| ≤ rad/2 < 0.025, where λ = 0; its certificate is in §4 (I1).
  Inverse: given (t, a₁), ρ by the certified-implicit inversion; then ζ = e^{ρ}·(phase from the S¹-step).
  On the central fibre (t = 0): the closed-form logistic inverses of §3.2 (with the node correction, exactly m₀ =
  a₀X − Y for rad ≤ 0.02: X = m₀/a₀ on D₁, Y = −m₀ on D₀).
* **Middle fibre coordinates (closed):** a₁ = x₄ − ξ₄(t) − ½ (|t| ≤ 1), a₂ = −x₄ − ½ (|t| ≥ 1), a₂′ = −a₂ +
  λ_∞(|t|)·δ(t) + (1 − λ_∞)·½, δ = arg t/2π − ξ₄(t) (single-valued), λ_∞ = 1 − s((|t| − 1.2)/0.4). Inverses:
  x₄ = a₁ + ξ₄ + ½; x₄ = −a₂ − ½ with a₂ = −a₂′ + λ_∞δ + (1−λ_∞)/2. ξ₄, (ξ₂, ξ₃) = lattice coordinates of σ̃(t)
  (Cramer, closed).
* **Cusp ∞ (|t| ≥ 2):** the same as B₀ with the ∞-data (s = 1/t, ν = e⁻², λ_∞ = 0).
* **h_B: (t, a) ↦ u ∈ S³ (closed).** T = t·k(|t|) (L4). |t| ≤ 1: u₂ = e^{2πia₁}/√(1+|T|²), u₁ = T·u₂;
  |t| ≥ 1: u₁ = e^{2πia₂′}/√(1+|T|⁻²), u₂ = u₁/T. Inverse: |T| = |u₁|/|u₂| ⇒ r = |T|⁻¹ (L4, implicit(1D));
  t = (u₁/u₂)/k(r); a₁ = arg u₂/2π, a₂′ = arg u₁/2π (mod 1: the fibre coordinates are circle-valued). The two
  descriptions agree on |t| = 1 (the clutching identity, `check_hopf` item; symbolic: a₂′ − a₁ = δ there by
  the definitions, since λ_∞ = 1).
* **Straightening R (closed given r).** R(u) = (φ̄₁(r)u₁, φ̄₂(r)u₂), unit phases φ₁, φ₂ depending on r only;
  R preserves |u₁|/|u₂|, so r is read off from R(u) as from u. Inverse: multiply by (φ₁, φ₂). The phases are
  fixed data (currently from the Chebyshev arc model; any smooth choice gives a valid diffeomorphism). det = 1.
* **moebius (closed, L6)** and **P_final (closed except the implicit(1D) ν⁻¹ on ¼ < −log|v| < 1, L3/L5).**
* **Λ_mod near the fixed points.** v ↦ (1 − c(rad))·M₀v_q + c(rad)·v, c = s((rad − 0.004)/0.004), with v_q =
  Hopf(√a₀x̄₀, y₀) and M₀ = 2πQ₀ (M_∞ = πQ_∞), exact signed permutations (sympy, `symbolic_identities.py`).
  Inverse: for rad ≤ 0.004, v_q = M₀⁻¹v (closed, and the lift is linear there, §3.4); for rad ≥ 0.008, Λ⁻¹ as
  above; in the shell, a **3-variable implicit inverse** whose certificate is the Jacobian bound of §4 (I2).

### 3.4 The lift Ψ_W and its inverse

* **Linear balls (closed).** rad ≤ 0.004: Ψ = √(2π)·U₀(x̄₀, y₀), X = −√(1 − 2π|q|²), U₀ = [[1, 1], [i, −i]]/√2;
  at ∞: Ψ = √π·U_∞(|ν|x̄′, y′), U_∞ = [[1, 1], [1, −1]]/√2. Inverse: q = U₀⁻¹Z/√(2π) (U₀ unitary).
* **Lift A / lift B (closed given b and the functions E, G_A at the point):** formulas of L7. Inverse of the
  lift = the S¹-phase step of §3.2.
* **Shell gluing Ψ = e^{ic(θ_v)·Arg g}·Lift_A**, g = ω_E·ω̄_A·conj(phase(v₁ + iv₂)) (closed given E, G_A).
* **E and G_A** (closed, transcendental): E = χ-weighted ψ̂-sums of y₀ + e0coef·x̄₀ at the cusps, the linear blend
  of E0_W(RP) and Einf_W(1/RP) in the middle; Ĝ_A = y₀ − x̄₀ | e^{2πix₁} | y₀′ − |ν|x̄₀′ glued by weighted logs;
  G_A = e^{iψ(|t|)}Ĝ_A. Their non-degeneracy conditions are the global ones, §5.

**Summary of G.** Given q = (Z, X, e^{2πiφ₂}, e^{2πiφ₃}): v = Hopf(Z); w = P_final⁻¹(v, X) [L5, L3];
u′ = moebius⁻¹(w) [L6]; r = |T|⁻¹(|u′₁|/|u′₂|) [L4]; u = R⁻¹(u′; r); (t, a) = h_B⁻¹(u); the W-reference point
p_ref from (t, a) [middle: closed; cusps: μ̃-inversion, implicit(1D); central fibre: closed logistic]; the
S¹-phase from ⟨Z, Ψ_W(p_ref)⟩ (closed); the T²-phase (x₂, x₃) = φ + (1−χ)D (closed). Near the fixed points,
inside rad ≤ 0.004: q = U⁻¹Z/√(scale) directly (closed). **Every step is closed-form or implicit(1D) with a
manifest certificate, except the two Λ_mod shells (3D implicit, §4) and the smoothness of Ψ_W(p_ref) in b, which
is exactly the global conditions (§5).**

## 5. The global conditions (I3–I5): analysis and the architecture that reduces them

### 5.1 Why they are not identities

The lift needs a nowhere-degenerate equivariant choice of a point on each target circle. The current design uses
two local lifts: lift B, smooth where E ≠ 0 and off Γ′ := Λ⁻¹(antipodal meridian); lift A, smooth where G_A ≠ 0 and
off Γ_E = Λ⁻¹(M₁); glued on the shell 15° ≤ θ_v ≤ 45° by e^{ic·Arg g}. The conditions are:
* **I3:** E ≠ 0 on {θ_v ≥ 15°};  **I4:** G_A ≠ 0 on {θ_v ≤ 45°};  **I5:** |Arg g| < π on the shell.

Each is an inequality for a transcendental function on a compact 3-dimensional set. Symbolic algebra cannot decide
such inequalities (Richardson), so they are either (a) reduced by design to statements that hold by structure, or
(b) certified by explicit analytic bounds, or (c) certified by validated numerics. The numerical margins are
large (|Arg g| ≤ 1.62; the arc within 0.056° of M₁; |arg t| ≥ 134° on the tube).

### 5.2 What is structural, and what is not (analysis of 2026-10-07)

* **Zero sets of ψ̂-sums at the cusps (toward I3).** With w_n = y₀⁽ⁿ⁾ + c·x̄₀⁽ⁿ⁾ the local model at ψ̂ⁿp, the exact
  ratio of consecutive terms is w_{n+1}/w_n = (c + |m|²|y|²t)/(m̄(|y|² + c t̄)) (from ψ̂(x, y) = (1/(my), mxy²)).
  On the overlap of two χ-supports (μ̃ + n ∈ (−0.9, −0.1), the D₀ side, where |y|² ≳ 0.1 for small t), both factors
  have small arguments when |t| ≪ |y|², so the two active terms differ in phase by less than π and cannot cancel:
  E = 0 forces w_n = 0, i.e. the point lies on the arc {|y₀|² = c|t|, arg t = π}. For small |t| this is a clean
  argument with explicit bounds (|arg(|y|² + ct̄)| ≤ arcsin(c|t|/|y|²), |arg(c + |m|²|y|²t)| ≤ arcsin(|m|²|y|²|t|/c),
  |arg m| ≤ |Im(c̃·α)|). For |t| ∈ [0.1, 0.6] the overlap can have |y|² < |t| and the argument needs the full
  series: not structural. **Status: I3 is structural only for small |t|; the band and the middle blend are not.**
* **The sign of the local model is forced.** For E₀ = y − cx̄ the ratio is ≈ −c/(m̄|y|²): anti-aligned, hence
  spurious zeros (notes §12.2). So a ψ̂-sum gauge with its zero arc at arg t = 0 does not exist; two ψ̂-sums with
  zero arcs at arg t = π ± β₀ (local models y + e^{±iβ₀}x̄) are both admissible for |β₀| < π/2.
* **The middle region admits a manifest zero set.** With E = e^{2πix₁}·f(σ, x₄), f(s, x) = Σ_n e^{−π(s+n)²/κ}
  e^{−2πinx} = e^{−πs²/κ}·ϑ(−x + is/κ; i/κ) (a Jacobi theta function in z = −x + is/κ, times a positive factor),
  f has **exactly one simple zero per cell**, at (s, x) = (½, ½), by the classical zero count of degree-1 theta
  functions. The twisted periodicity f(s+1, x) = e^{2πix}f(s, x) is what the branch change x₁ ↦ x₁ − x₄ requires.
  (This was the §8 middle; the problem was only the *gluing* to the cusp sums, which was done through a log-ratio
  at a common zero.)
* **The fixed points force exact linearity of Λ in v_q** (notes §13.3.1). This is re-derived in §4: with
  Λ = M₀v_q + δ, δ = O(|q|⁴), lift B gives Z₁ = ω_E·√(4π|e|² + (v₃/|v|)δ₃ + …), not smooth across the arc, and
  lift A gives Z₂ = √(2π)(U₀q)₂·√(1 + δ-stuff/|(U₀q)₂|²), only C^{1,1} at p₀. Hence the linear blends.

### 5.3 A simpler lift that removes the gauge normaliser: the minimal-rotation lift (design, not implemented)

Let 𝔈 = (E, G_A): W∖Fix → ℂ²∖0 (weight 1; it is nonvanishing iff E and G_A have no common zero) and let
n̂_𝔈 = Hopf(𝔈)/|𝔈|² ∈ S², a smooth map of B. Let n̂_v = v/|v| be Λ's direction. Define
  **Ψ := (√|v| · U(b) · 𝔈/|𝔈|, X)**, U(b) ∈ SU(2) the lift of the minimal rotation taking n̂_𝔈 to n̂_v
  (axis n̂_𝔈 × n̂_v, angle arccos(n̂_𝔈·n̂_v); U = cos(α/2)·1 − i sin(α/2)·(k̂·σ)).
* Equivariant; Hopf(Z) = R_U n̂_𝔈·|v| = v, |Z|² = |v| (L7-type identities). Smooth wherever 𝔈 ≠ 0 and α < π.
* **At the fixed points it is exactly the linear model:** there 𝔈 = U₀q and v = M₀v_q = 2πQ₀v_q, so n̂_𝔈 =
  Q₀v_q/|v_q| = n̂_v, U = 1, Z = √(2π)·U₀q. No ψ, no shell, no Chebyshev fit of ψ.
* The conditions become **(C1)** E and G_A have no common zero, and **(C2)** n̂_𝔈·n̂_v > −1 everywhere. (C2) is
  implied by |Arg g| < π wherever both lifts are defined (the current I5), and is strictly weaker: antipodality
  needs the phase condition *and* |E|²−|G_A|² = −v₃|𝔈|²/|v|. Near Γ_E: n̂_𝔈 = (0,0,−1) = n̂_v (aligned); near the
  zeros of G_A: n̂_𝔈 = (0,0,1), antipodal only on Γ_E, excluded by (C1).
* So I4 and I5 collapse into one scalar inequality (C2) on the compact B, and I3 into (C1). Both remain global
  inequalities (type (b)/(c) of §5.1), but the construction loses one fitted function and two gluing devices.

### 5.4 The "frame" alternative, and why the orbit-map-first architecture stays

Normalising a frame directly, Ψ = (√(1−X²)·𝔈/|𝔈|, X) with a height X on B, is the §3 plan of the notes: it
trades the lift conditions for "on each level sphere of X, [E : G_A]: Σ_c → ℂℙ¹ is a diffeomorphism", which is
harder (the level spheres of a suspension structure on B are transverse to the torus fibration by |t|). The
orbit-map-first architecture, where Λ's bijectivity is a composition of one-variable monotone maps and classical
diffeomorphisms, is the right one for a symbolic treatment.

## 4. The fixed points: exact facts, and the certificates for I1 and I2 (`sym/node_remainder.py`, `sym/i2_margins.py`)

### 4.1 Exact facts (sympy, N1–N2)

* **φ(t) = 2|t|²/(1 − |t|²) at both cusps**, exactly (c̃ = R(t)⁻¹ē₄, φ = Re(c̃₂α₂ + c̃₃α₃); det R = ±(1 − |t|⁴)).
  Hence P(t) = log|ν| + log|t| − 2|t|²/(1 − |t|²) and m(t) = ν·e^{−(c̃·α)} are closed-form; P < 0 on |t| ≤ ½ for
  both cusps (log ½ < 0; at ∞ also log|ν| = −2), manifestly.
* **Node variables.** X̂ = a₀|x₀|², Y = |y₀|², rad = X̂ + Y, |t|² = X̂Y/a₀ ≤ rad²/(4a₀). The chart series is
  μ̃ = N/D, N = Σ_j jT_j, D = Σ_j T_j, **T_j = e^{−j(j+1)φ}·X̂^{j(j+1)/2}·Y^{j(j−1)/2}**, of total degree j² in
  (X̂, Y), with e^{−j(j+1)φ} ≤ 1 (j(j+1) ≥ 0, φ ≥ 0). The node height is m₀ = X̂ − Y (= a₀|x₀|² − |y₀|²).
* **Identity (N2d):** N − m₀D = X̂(e^{−2φ} − 1) − (X̂ − Y)(e^{−2φ}X̂ + Y) + R₁ − (X̂ − Y)R₂, with
  R₁ = Σ_{|j|≥2} jT_j, R₂ = Σ_{|j|≥2} T_j. Expansion: μ̃ − m₀ = Y² − X̂² + O(3) in (X̂, Y) (N2b′); the two-term
  truncation (e^{−2φ}X̂ − Y)/(1 + e^{−2φ}X̂ + Y) differs from μ̃ by O(4) (N2c).

### 4.2 I1: the node correction is monotone, with explicit constants

a₁ = μ̃ + c(rad)(m₀ − μ̃), c = 1 − s((rad − 0.02)/0.03). Along a torus fibre (ρ = log|y₀|, t fixed), with
′ = d/dρ (N3): d(rad)/dρ = −2m₀, m₀′ = −2rad, μ̃′ = −2Var, and
  **a₁′ = (1 − c)μ̃′ + c·m₀′ − 2c′(rad)·m₀·(m₀ − μ̃).**
The first two terms are manifestly ≤ −2rad·min(Var/rad, 1); the last is the only non-manifest one.
**Bound chain (rad ≤ R₀ = 0.05; explicit constants).**
1. |c′| ≤ 2/0.03 (L1), |m₀| ≤ rad. So |2c′m₀(m₀ − μ̃)| ≤ (4/0.03)·rad·|m₀ − μ̃|.
2. From N2d with D ≥ 1: |μ̃ − m₀| ≤ 2φX̂ + rad² + |R₁| + rad·R₂. Tails: Σ_{|j|≥2}T_j ≤ 2Σ_{k≥2} rad^{k²} ≤
   2rad⁴/(1 − rad⁵), Σ_{|j|≥2}|j|T_j ≤ 2Σ_{k≥2} k·rad^{k²} ≤ 4rad⁴/(1 − rad⁵)·(1 + …) ≤ 4.01rad⁴.
   φ ≤ 2|t|²/(1 − |t|²) with |t|² = X̂Y/a₀, and X̂²Y ≤ (4/27)rad³: 2φX̂ ≤ (16/27)(1+ε)rad³/a₀.
   ⇒ |μ̃ − m₀| ≤ C·rad² with C = 1 + (16/27)(1+ε)R₀/a₀ + 4.01R₀² + 2.01R₀³: **C ≤ 1.04 (cusp 0), C ≤ 2.63 (cusp ∞)**
   (sampled maxima 0.999 and 1.871, N5).
3. Var ≥ (e^{−2φ}X̂ + Y)/D − μ̃² ≥ rad(1 − 2φ)/(1 + rad + 2rad⁴) − rad²(1 + C·rad)² ≥ rad(1 − C₂rad) with
   **C₂ ≤ 2.2 (cusp 0), 2.8 (cusp ∞)** (sampled 2.0, 2.5).
4. Hence −a₁′ ≥ 2rad(1 − C₂rad) − (4/0.03)·C·rad³ = 2rad[(1 − C₂R₀) − (2/0.03)·C·R₀²] ≥ 2rad·[0.89 − 0.17] (cusp 0)
   and 2rad·[0.86 − 0.44] (cusp ∞): **a₁′ < 0 on the whole node ball**, by explicit inequalities only.
   (Sampled minima of −a₁′/(2rad): 0.852 and 0.826, N5.)

So **I1 is certified symbolically** (series tail bounds with positive terms; no interval arithmetic), given the
constants above. What is machine-checked: the identities N2d, N3, the term structure N2e, and the sampled margins;
the inequality chain 1–4 is elementary and written out here.

### 4.3 I2: the linear-blend shell of Λ_mod

Λ_mod = (1 − c(rad))M₀v_q + c(rad)Λ(v_q), c = s((rad − 0.004)/0.004). Identity (N4):
  **DΛ_mod = (1 − c)M₀ + c·DΛ + (Λ − M₀v_q)⊗∇c**, so
  σ_min(DΛ_mod) ≥ scale − c·|DΛ − M₀| − |∇c|·|Λ − M₀v_q|, |∇c| ≤ 2/0.004 = 500 (L1), scale = 2π (resp. π).
Λ is real-analytic near b₀ with Λ(b₀) = 0 and DΛ(b₀) = M₀ exactly (sympy, `symbolic_identities.py`), so
|Λ − M₀v_q| ≤ C_Λ|v_q|² and |DΛ − M₀| ≤ C′_Λ|v_q| on the shell, with constants that are explicit second-derivative
bounds of the closed-form composition h_B → moebius → P_final in the variables (t, a₁) — derivable, not yet
written out. Required: 500·C_Λ·(0.008)² + C′_Λ·0.008 < scale. Measured constants and margins: see the table
appended below (`sym/i2_margins.py`).

**I2 margins (`sym/i2_margins.py`, 4000 points per shell, rad ∈ [0.004, 0.008]):**

| cusp | max \|Λ − Mv_q\|/\|v_q\|² | max \|DΛ − M\|₂/\|v_q\| | scale | bound from N4: scale − max\|DΛ − M\| − 500·max\|Λ − Mv_q\| | measured min σ_min(DΛ_mod) |
|---|---|---|---|---|---|
| 0 | 39.5 | 79.0 | 2π | **4.42** | 5.64 |
| ∞ | 19.7 | 39.5 | π | **2.20** | 2.82 |

The constants scale exactly as a quadratic remainder (C′_Λ = 2C_Λ). So the certificate for I2 needs only
C_Λ ≤ 60 and C′_Λ ≤ 120 (cusp 0; half that at ∞), i.e. a bound on the second derivatives of the closed-form chain
on the ball |v_q| ≤ 0.008, three times the measured values. **I2 is reduced to an explicit second-derivative bound
of a closed-form map on a small ball**, the same type of statement as I1, with the identity N4 machine-checked.

## 6. Status (2026-10-07) and resume checklist

### 6.1 The lemma list of the proof scheme (notes §13.1), with its status

| # | statement | status |
|---|---|---|
| A | Y's charts, gluings, lattice coordinates, Φ, Θ, ψ̂, T³-actions: identities | symbolic (the paper's `verify_construction.py`; here C1–C3, C6) |
| B | lattice system nonsingular: det = (\|t\|²+1)(log\|t\|(\|t\|²−1) + 2\|t\|²) > 0 | **manifest** (C1d) |
| C | μ̃ strictly monotone along fibres (dμ̃/dρ = −2Var), range ℝ, equivariant, closed on the central fibre | **manifest** (C4, C5) |
| D | λ-blend a₁ = (1−λ)μ̃ + λ(x₄^slice − ½) monotone, no wrap needed | **manifest** (§3.3; P < 0 by §4.1) |
| E | node correction monotone (I1) | **explicit-constant certificate** (§4.2; identities N2d, N3 machine-checked) |
| F | \|T\| strictly increasing (needs s′ ≤ 2) | **proven** (L1, L4) |
| G | ν strictly increasing (P_final) | **manifest** (L3) |
| H | R a diffeomorphism (fibre twist, det 1), moebius, P_final zones: closed inverses | symbolic (L5, L6) |
| I | exact linear models M₀ = 2πQ₀, M_∞ = πQ_∞ | symbolic (`symbolic_identities.py`) |
| J | linear-blend shells nonsingular (I2) | **reduced to an explicit 2nd-derivative bound** (§4.3; N4 machine-checked; margins 3×) |
| K | lifts land in S⁴, cover Λ, equivariant; linear in the balls | symbolic (L7, §3.4) |
| L | E ≠ 0 off the tube (I3); G_A ≠ 0 on the tube (I4); \|Arg g\| < π on the shell (I5) | **global inequalities, not symbolic** (§5); structural for small \|t\| at the cusps and for the middle theta section; the minimal-rotation lift (§5.3) merges I4 + I5 into one scalar condition and removes ψ |
| M | φ: Y → T² equivariant, invariant, closed inverse | symbolic (C6, §3.2) |
| N | F∘G = id on each piece; pieces of G agree on overlaps | follows block by block from A–K once L is certified: every block's inverse is closed or certified-implicit |
| O | DG invertible ⇒ J = G*J_Y smooth, J² = −1, N = 0 ⇒ Newlander–Nirenberg | one-line arguments (notes §13.1) |

### 6.2 What this session established

* The symbolic specification of F and G is complete at the level of blocks (§3): every block is closed-form or a
  certified-implicit one-variable inversion, except the two Λ_mod shells (3D implicit, certificate J) and the
  lift's dependence on E, G_A (L).
* Three lemmas of notes §13.3 are now proven or manifest (F, G, D); one is unnecessary (the wrap bound).
* New exact facts: φ(t) = 2|t|²/(1−|t|²) at both cusps; the lattice determinant; the node-variable form of the
  theta series with degree j² terms; the two-term truncation accurate to O(rad⁴).
* I1 and I2 are reduced to explicit Taylor-remainder bounds with large margins (I1 written out; I2's constants
  measured, bound to be written).
* The honest status of the global conditions: **I3–I5 remain inequalities over a continuum**; a purely symbolic
  proof of the whole theorem is not available with this architecture, and §5 explains why. They are certifiable by
  validated numerics (notes §12.12) or by analytic bounds in the structural regions.

### 6.3 Not done / next steps

1. Write out the I2 second-derivative bound with explicit constants (same style as §4.2).
2. Implement the minimal-rotation lift (§5.3) in torch, rerun `check_lift`/`check_degree`/`check_full`, and measure
   (C2) = min n̂_𝔈·n̂_v; if its margin is large, retire ψ, `psi_fit.json` and the shell.
3. Decide the treatment of L: analytic bounds for small |t| (the aligned-phase argument with the explicit ratio of
   §5.2) plus validated numerics on the compact remainder; or accept the numerical certificate.
4. Assemble a sympy "F∘G = id" script per piece (N), using the block inverses of §3 (mechanical once L is fixed).

### 6.4 Resume checklist

Read this file top to bottom (20 min). Run from the repo root: `python3 extra/s4t2_explicit/sym/lemmas_1d.py`,
`.../chart_identities.py`, `.../node_remainder.py` (each ≤ 1 min, all PASS) and `.../i2_margins.py` (2 s).
Then continue with §6.3. Standing constraints: notes §13.5 (no agents unless asked; commit only on request).

## 7. The certificate goal (user, 2026-10-07) and what "100% symbolic" can mean

**Goal (user):** a symbolic implementation of the (1,1)-tensor J on S⁴ × T² with checks of (1) well-definedness,
(2) J² = −1, (3) N_J = 0, so that correctness reduces to "the object is correctly defined and the checks are
correctly implemented" — a reviewable certificate, without Lean and without reading the paper in full.

**Assessment.**
* The format is right, with one precision: for J = G*J_Y, checks (2) and (3) hold for *every* smooth G with
  DG invertible (naturality; notes §12.12–12.13, toy). They are consistency checks of the implementation and are
  best verified once, generically (J = A J₀ A⁻¹ with A a Jacobian of undetermined functions), not per piece. The
  whole theorem sits in (1): G smooth on each piece, DG invertible, pieces agreeing on overlaps, domains covering.
  Newlander–Nirenberg is cited, not re-proved.
* "Well-defined" for formula-defined tensors means "every denominator nonzero and every implicit inversion
  certified" — exactly the lemma list of §6.1. The theorem has quantitative content (an explicit map is a
  diffeomorphism), and quantitative content is inequalities. So **identity-only is impossible in principle**; the
  achievable optimum is **identities + a short list of elementary inequalities with explicit constants** (series
  with positive terms, convex combinations, one-variable bounds), each verifiable by hand or by exact rational
  arithmetic. §4.2 (I1) is the model.
* With the current architecture (two lifts glued by Arg) this optimum is not reachable: I4 and I5 are intrinsically
  global (a winding-1 transition between two trivialisations, matched by a partition of unity).
* **Architecture (C), single compatible lift — the candidate for the optimum:** one weight-1 function E with a
  manifest zero set and linearisation compatible with Λ along the whole arc; then lift B alone is smooth off the
  poles (it is smooth across Γ′ since its only denominator vanishes on Γ_E), the gauge, ψ, the shell, `fit_psi`
  disappear, and `fit_arc` disappears because the arc is explicit. Ingredients: (i) in the middle, E =
  e^{2πix₁}·ϑ with an r-dependent modulus τ(r) ∈ ℍ and explicit zero (one simple zero per cell, classical);
  (ii) at the cusps, the single-term local model y₀ + x̄₀ near the arc (exact for |t| ≤ 0.3) with zero at
  {|x₀| = |y₀|, arg t = π} and the rigid induced structure z_E = log|y₀| + (i/2)arg t (independent of the model's
  coefficients α(t), β(t)); (iii) gluing cusp ↔ middle by a log-ratio on an annulus where both have the same zero
  and the same structure (τ(r) matched to z_E there), with the ratio's windings computed from the explicit
  local model; (iv) compatibility of Λ with E along the arc by an explicit anti-conformal twist
  (v₁+iv₂) ↦ A(λ)(v₁+iv₂) + B(λ)·conj(v₁+iv₂) in the w-chart (closed-form A, B from derivatives of closed-form
  maps; a diffeomorphism iff |B| < |A|, a one-variable sign condition along the arc). Residual inequalities:
  I1, I2 (explicit constants, done/structured), I3-cusp = E ≠ 0 off the arc on |t| ≤ 0.3 (the aligned-phase
  argument; reducible to explicit-constant series bounds by narrowing the χ-supports so that the overlaps have
  |y₀|² manifestly > |t|), and the 1D sign condition of (iv). Everything else identities.
* Fallback for any residual inequality that resists: validated numerics (notes §12.12), stated as such.

**What the symbolic J would be:** per piece, J(q) = DF(G(q))·J₀·DG(q), with the implicit 1D inverses as symbols
carrying their defining identities and IFT derivative rules; checks: F∘G = id per piece, DF(G)·DG = 1, overlap
agreement through the chart transitions, domain covering (manifest unions), the generic naturality lemma, and the
lemma list of §6.1 with explicit constants.
