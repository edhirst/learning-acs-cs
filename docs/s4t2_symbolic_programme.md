# The symbolic programme: F: Y → S⁴ × T² and its inverse G in symbolic form (started 2026-10-07)

> **Status line (keep current).** Sections 0–6 written on 2026-10-07; all four scripts in `extra/s4t2_explicit/sym/`
> PASS. Review additions: I1 and I2 rational certificates, C7 (unwrapped blend), C8 (closed pole-cap inverse),
> C9 (complete tensor on the caps), C10 (six-dimensional naturality), and L8 (smooth division;
> first-order compatibility is insufficient). Numerical facts carry their source script.
> **Global expression:** §8 now supplies G_Q, its inverse and J_Q by exact one-dimensional
> quadratures, with C11–C15 and the shell contraction I2o–p. This is a new parametrization
> of Y, independent of the fitted lift's I3–I5; it is not a finite elementary closed form.
> The numerical additions evaluate its base/connection blocks, not the complete global tensor.
> **Closed-form continuation:** §9 eliminates both scalar profile inverses, certifies the revised
> shells, and gives a genuinely algebraic tensor on two pole neighbourhoods after a torus gauge
> change. The global height inverse, shell inverse and transport/phase quadratures remain.
> There is a proof excluding algebraic equivariant heights, not a proof excluding elementary J.
> **Literature-led continuation:** §10 supplies a finite elementary forward height, an alternative
> height/collar with a standard elliptic inverse, an algebraic S⁴→S² target, and finite bundle
> operations with explicit hypotheses. The elliptic block is not yet assembled into a global
> base/lift/tensor. The finite bundle operations still need a concrete embedding and derivative bounds.
> **Finite global tensor (2026-10-08-EH-2):** §11 supplied a different, elementary inverse-first
> height, smooth open central-fibre charts, a certified closed pole replacement, a concrete C¹⁰
> bundle embedding and fixed-count projection/phase formulas. Equations (11.18)–(11.19) give
> G_cl and J_cl with no implicit solve, quadrature or infinite series. The counts are enormous;
> this is a finite elementary expression, not a practical expanded formula. The global statement
> requires human review. Section 11.3 gives a finite elementary inverse for the pole patch.
> **Finite inverse pair (2026-10-11-TSG-6):** §11.2 gives a two-way elementary height,
> §11.4 its revised derivative certificate, and §11.6 the global reversal F_cl.
> C29a–h, I2v, independent audits of the type-(c) constants and forward formula,
> 276 sampled height compositions and 38 reduced-count end-to-end checks pass.
> The finite elementary inverse-pair certificate is complete, subject to human
> review. Its enormous fixed products have not been numerically executed.
> Legacy I3–I5/T14 remain flags for the old evaluator.
> **Direct second-jet audit (2026-10-12-TSG-1):** C30 checks exact-count
> projection/phase derivative steps and cusp restoration; N8 differentiates the
> actual map on both inner pole caps, and the pure-middle check computes its
> actual local tensor. All pass. The exact free-slab and seam branch circuit
> needed for a direct check of the specific global \(J_{\mathrm{cl}}\) is open.
> **Further regional audit (2026-10-12-TSG-2):** C31 and C32 follow an actual
> pure-middle projector cell and a two-active-pair cusp projector cell with
> their second jets; N9 identifies the source cap projector and checks both
> cap/slab overlaps using the stated bounds;
> the active south collar check differentiates its finite torus correction;
> and a south raw-pole annulus check constructs a six-real-coordinate
> second-jet circuit for the complete local map and its specific \(J\).
> These pass on their stated cells, but T16's complete global branch circuit
> is still open.
> **Radial-stage audit (2026-10-12-TSG-3):** N10 checks the exact elementary
> radial compression and its inverse, their source-linked implementation
> formulas, all 16 first and 64 second derivatives in each direction, and
> the positive Jacobian factor on the nonidentity pole shell. The full
> shell map and its six-dimensional tensor remain to be composed and checked.
> **North annulus audit (2026-10-12-TSG-4):** the north raw pole inverse,
> round cap and finite torus restoration give a second source-linked local
> six-real \(J\) circuit, with exact value, first and second derivative
> nodes. Its compositional square and Nijenhuis checks pass on the stated
> annulus; the radial shell and global free slab remain open.
> **Dominant-profile audit (2026-10-12-TSG-5):** N11 checks the ordered
> raw-reference source path through both positive-root ratio branches,
> exact round-meridian witnesses, their finite first and second jets, and
> the removable dominant-coordinate selector seam for the effective \(t\).
> It is a scalar generator check; the remaining projector and free-slab
> composition remain open.
> **Cusp projector and log-seam audit (2026-10-12-TSG-6):** C33 checks
> the negative-height shifted cusp frame and its second jets at both
> cusp signs. C34 checks a negative-real-log seam in the active south
> and north collars through the source Cramer coordinates, angular
> phase, cusp reconstruction and their local chart-coordinate jets.
> These are regional projector checks;
> the complete free-slab map remains open.
> **South transition audit (2026-10-12-TSG-7):** the two nonidentity
> radial maps, south raw inverse, round cap and torus restoration now
> compose to a source-linked six-real \(J\) circuit on one strict open
> radial-transition chamber. Its compositional second-jet and
> \(J^2/N_J\) checks pass; the rest of the shell and free slab remain open.
> **All-mesh free-slab audit (2026-10-12-TSG-18):** a strict near-south
> chamber now has the source's actual fixed-count projection expressed as
> a finite indexed product, with a six-real value, Jacobian, six Hessians,
> and named \(J=M^{-1}J_0M\). The local \(J^2\) and Nijenhuis reductions
> pass. Other free-slab cells and global overlap coverage remain open.
> **Cusp selector audit (2026-10-12-TSG-19):** C38 checks an open
> argument-selector overlap at both cusps. The two source references
> and projectors agree with their second jets; the restored chart
> change is an explicit holomorphic deck map. This is one selector
> family, conditional on nonzero projection overlaps elsewhere.
> **Evaluator branch audit (2026-10-12-TSG-20):** an opt-in source
> mode now follows the six exact finite-formula branches where the
> ordinary numerical evaluator takes working-precision identity
> shortcuts. Its source guards and context restoration pass; global
> T16 coverage remains open.
> **Dominant-selector overlap (2026-10-12-TSG-21):** C39 checks
> an open band across the raw reference's dominant-coordinate
> selector. Its two positive-root formulas give the same effective
> \(t\), source reference, projector and restored output through
> second order. Exact round witnesses lie on both selector sides;
> other source seams and global T16 coverage remain open.
> **North flat-selector join (2026-10-12-TSG-22):** C40 checks
> one north angular/height chamber at the source
> \(|\mathrm{small}|=1/2\) selector. The infinity-side radial
> and embedding weights are flat there. The final
> phase/restoration formulas have matching second jets
> modulo deck conditional on matching finite-fold endpoint
> jets and nonzero phase. The intermediate projector joins
> and nonzero projection denominators remain premises.
> **Cusp/middle join (2026-10-12-EH-1):** C41 checks the
> source selector \(|t_c|=3/10\) at both cusp signs along the
> whole seam, not in one chamber. The core reference is a unit
> multiple of the \(\lambda=0\) middle continuation, so the
> ten-slot projectors agree identically. The right-hand formulas
> differ from that continuation only by flat terms. For every
> \(\phi\), the c0/ci restoration is the paper's holomorphic
> chart change \(\Phi\)/\(\Theta\) of the middle restoration.
> Propagation through the fold still assumes nonzero projection
> denominators and checks of the other seams.
> **Pure-middle seams (2026-10-12-EH-2):** C42 proves one lemma
> on every open subset of \(1/2<|t|<2\). Modulo the moves
> \(x_1\mapsto x_1+1\), \(x_4\mapsto x_4+1\) and the log-branch change,
> the source \(W\)-data equal data built from continuous branches
> of \(\arg u_1\), \(\arg u_2\), \(\arg t\) and \(\log t\). The
> projector is therefore unchanged. For a final fold vector that is
> smooth there (a premise resting on the other seams), the restored
> output changes by holomorphic lattice translations. This covers the
> \(|t|=1\) selector and all argument cuts, including their corners.
> **Argument selector beyond the pure cusp (2026-10-12-EH-3):** C43
> extends C38 to every \(|t_c|<1/2\) at both signs. On the lower
> overlap, where the two charts differ by \(a\mapsto a+1\), they give
> the same ten-slot projector. For a smooth final fold vector, the
> restored outputs differ by the cusp deck \(T\) in the core and by
> \(\pm e_4\) in the collar. The key new identity is an exact deck shift
> of the four-seed height on \(1/4<Y<2\) for every \(q<1/8\), replacing
> C38's computation for \(q<1/100\).
> **Collar/middle joins (2026-10-12-EH-4):** C44 checks both
> collar/pure-middle selectors, \(|t|=1/2\) (south) and
> \(|\mathrm{small}|=1/2\) (north), along their whole seams. At
> \(\lambda=1\) the collar height is exactly \(\rho=P(a+1/2)\). The
> collar and pure-middle \(W\)-data then have the same \(x_4\) modulo an
> integer, and \(x_1\) differing by a circle gauge. So the projectors
> agree, and for a smooth final fold vector the outputs differ by a
> lattice translation. N13 and C40 are special chambers of C44.
> **Height moment and cells (2026-10-12-EH-5):** C45 certifies the
> inverse-first height moment and the cusp height cells as one smooth
> function. Its stable rewrites (in \(R\), \(R^{-1}\) and the two branches
> of the terminal map) are exact identities. The joins at \(|b|=2/5\)
> and \(1/5\) are flat. The core and collar height-cell seams are
> equalities of analytic cell formulas on open intervals.
> **Pole-cap internals (2026-10-12-EH-6):** C46 certifies the whole
> cap branch at both signs as one smooth map. Its 32 conditionals are
> classified. The raw inverse is a single branch with identical
> rationalized alternatives, and its guards are inactive with explicit
> margins. The production radial maps are exactly
> \(\psi^{-1}T\varphi\) and \(\varphi^{-1}T^{-1}\psi\), with flat-factor
> forms at their outer seams.
> **Phase band (2026-10-12-EH-7):** C47 certifies the slab phase band
> at both signs. The band cutoff is a \(C^\infty\) join. The mismatch is
> independent of the hemisphere Hopf reference. The two stereographic
> charts of the finite logarithm agree exactly, and its start constant is
> one fixed point per band.
> **Direct certificate (2026-10-13-EH-1):** §11.9 assembles C38–C50,
> N9 and C30i/C10 into a global statement: under the listed written
> hypotheses, \(J_{\mathrm{cl}}\) is a well-defined smooth tensor on the
> round \(S^4\times T^2\), with \(J^2=-\mathrm{Id}\) and \(N_J=0\) exactly.
> No guard on the path fires (C48), \(DG_{\mathrm{cl}}\) is invertible
> everywhere (C49), and the pinned source is classified test by test,
> cut by cut and operation by operation (C50). It replaces T16's
> per-branch circuit by a generic identity, a change of route that needs
> the owner's agreement. This is a mathematical statement for human
> review.
> Companion of `docs/s4t2_explicit_diffeo.md` (§12 = the numerical construction, §13 = the agreed scheme).
> **Standalone implementation (2026-10-14-TSG-8):** §12.12 supplies
> `direct.complex_structure(point)`, its exact local expression graph, and
> `sym/direct_certificate.py`. The shared formula bodies, differentiated
> scalar assignments, global tangent frame and square/Nijenhuis identities
> are checked separately from the written global contracts.
> **Resume point: review the ten block contracts in §12.12 against §12.11 and its cited proofs; the new callable and focused verifier are the practical entry points.**

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
* **L8 (smooth extension across the arc; corrected criterion, type (a)).** Put f = v₁ + iv₂ and
  d = |v| − v₃ > 0 near the negative meridian, away from the fixed points. A sufficient condition is the
  **exact factorisation f = H E with H smooth and nonvanishing**, in a local unitary frame. Lift B then equals
  Z₁ = |H|E/√(2d), Z₂ = (H̄/|H|)√(d/2), including at E = 0. These formulas are smooth, have Hopf image v,
  and agree with lift B off the arc; the identities are checked in `sym/lemmas_1d.py` L8a–c.
  Matching transverse linearisations alone does **not** suffice. On |z| < ½ take E = z + z̄² and
  v = (2 Re(z)√(1−|z|²), 2 Im(z)√(1−|z|²), 2|z|²−1). Here E has only the simple zero 0
  (|E| ≥ |z|(1−|z|) > 0 off 0, type (b)), and df(0) = 2dE(0). But the second component of lift B is
  √(1−|z|²)·z̄E/(|z||E|); its directional derivatives at 0 along 1, i, 1+i are 0, i, −i, respectively
  (L8d–e), so it is not differentiable. Even Z₁ is not C²: its second directional derivatives along
  1, i, 1+i, 1−i are 0, −2, 2−2i, 2+2i, violating the parallelogram identity for a quadratic form (L8f).
  Thus T4 must establish smooth division on a neighbourhood of the
  arc. One way to enforce it is E_new = f/conj(E″) with E″ smooth and nonvanishing there; extending or
  gluing this choice globally remains open. This corrects the first-derivative criterion in notes §11.8.

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
* **B_mid coordinates:** (t, x₄ mod 1), with x₁ mod 1 the S¹-fibre and (x₂, x₃) mod 1 the T²-fibre.
  B_mid = ℂ* × S¹_{x₄}. The coordinates (x₁,x₄) on W_mid are local: changing log t by 2πi
  changes x₁ by −x₄. Thus x₁ = x₂ = x₃ = 0 is a reference on a chosen log branch, not a
  global product trivialization of W_mid (C2, C11d).

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
  T_j = e^{(j²+j)g} X^{j(j+1)/2} Y^{j(j−1)/2}, X = |x|², Y = |y|², g = log|ν| − φ(t).
  This series converges normally on compact subsets of the cusp domain |xy|<1, including
  its axes, by the explicit tail bound in §8.2; it is not asserted to be entire in (X,Y). Facts:
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
* **The S¹-phase (closed on a regular phase patch):** given the W-point's reference p_ref (x₁ = 0 on a
  middle log branch; y real positive in a cusp patch with y≠0), put W=Ψ_W(p_ref). For Z≠0 with
  Hopf(Z)=Hopf(W), set e^{iθ}=(Σ_j Z_j W̄_j)/|Σ_j Z_j W̄_j| and take e^{iθ}·p_ref.
  This uses the explicit pairing convention and gives x₁=θ/(2π) in the middle. C15b–c show that
  Z=λW with |λ|=1, hence the numerator is λ|W|²≠0. At the fixed points use §3.4's closed cap inverse.

*Correction to 3.1:* the determinant of the lattice system is det = (|t|² + 1)(log|t|·(|t|² − 1) + 2|t|²)
(sympy C1d), branch-independent and **manifestly positive on ℂ***: log|t| and |t|² − 1 have the same sign.

### 3.3 The orbit-space map Λ = P_final ∘ moebius ∘ R ∘ h_B and its inverse

* **Fibre coordinate on B₀ = {|t| ≤ ½} (implicit(1D) along fibres).**
  a₁ = (1 − λ(|t|))·μ̃_new + λ(|t|)·(x₄^slice − ½), λ = s((|t| − 0.3)/0.2), x₄^slice = ρ/P(t).
  **T7: no `wrap` (identities (a) and monotonicity (b), C7).** Both μ̃_new and x₄^slice − ½ are globally defined functions of
  ρ on the fibre, each shifting by exactly +1 under ψ̂, so the convex combination is equivariant without any
  integer choice. In fact the old wrapped formula agrees exactly with this formula for the infinite series.
  Set σ = ρ/P − ½ and a = −P > 0. Completing the square gives weights proportional to e^{−a(n−σ)²}.
  Their mean M(σ) satisfies M′ = 2a Var > 0. At every half-integer b, reflection n ↦ 2b−n gives M(b) = b.
  Thus σ and M(σ) belong to the same interval between consecutive half-integers, giving
  **|M(σ)−σ| < ½**. The Gaussian and its derivatives converge uniformly on compact subsets of a > 0,
  justifying differentiation and pairing; C7 checks the square completion, reflection, and derivative identities.
  The same bound holds for the code's symmetric finite truncation about k = round(σ): put δ = σ−k.
  For 0 ≤ δ ≤ ½, pairing j with −j shows the centred mean is ≥ 0, since their squared-distance
  difference is 4jδ ≥ 0. Pairing j with 1−j for j = 1,…,J shows the mean is < ½: the squared-distance
  difference is (2j−1)(1−2δ) ≥ 0, and the unpaired index −J contributes strictly negatively to the
  moment centred at ½. Reflect for δ < 0. These identities are C7h; hence wrap removal also agrees
  for finite centred sums in exact arithmetic. The smooth definition still uses the infinite series.
  The node correction is absent wherever λ > 0: its support has
  |t| ≤ rad/(2√a₀) ≤ e²/40 < 1/5 < 0.3 at either cusp (e² < 8 follows from the rational bound e⁴ < 55
  in §4.2). Hence `wrap(σ−μ̃_new) = σ−μ̃_new` throughout the blend annulus. The finite-series implementations
  are compared numerically in `check_hopf.py`; this comparison is a floating-point safety check. Then
  da₁/dρ = (1 − λ)·dμ̃_new/dρ + λ/P(t) < 0 manifestly where μ̃_new = μ̃ (P < 0 on |t| ≤ ½: P = log|ν| + log|t| − φ(t) with
  φ = 2|t|²/(1 − |t|²) ≥ 0 exactly (§4.1), so every term is ≤ 0 and log|t| < 0), since λ depends on |t| only. The node correction μ̃_new ≠ μ̃ only for
  rad = a₀|x₀|² + |y₀|² < 0.05, hence |t| < 0.025/√a₀ < 0.2 at either cusp, where λ = 0; its certificate is in §4 (I1).
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
  descriptions agree on |t| = 1 (the clutching identity, `check_hopf` item; symbolic: a₂′ − a₁ = arg t/(2π) mod 1 there by
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
  **Implemented closed inverse (type (a), C8):** on the target caps X < 0, |Z|² ≤ 2π·0.004 and
  X > 0, |Z|² ≤ π·0.004 respectively, set (q₁,q₂) = U*Z/√scale,
  x_s = q̄₁/√a₀, y_s = q₂, t = x_s y_s. Then H = μ̃(x_s,y_s), C = φ+Hc̃(t),
  x = x_s e^{−C·α}, y = y_s e^{C·α}, w = C₂ē₂+C₃ē₃.
  The linear lift, slice, and torus-trivialisation compositions are identities (C8); the slice determinant
  is nonzero since |t| ≤ e²·0.004/2 < 0.016 < 1 (type (c), using §4.2).
  `tinverse.linear_ball_inverse` evaluates this formula without Newton, including Z = 0.
  The torus arguments select local representatives; changing a representative by an integer is the T² lattice action.
  `check_inverse.py` tests both cap boundaries and coverage; `check_jround.py` tests J, DJ and finite
  differences at the exact poles and nearby axes. These numerical tests supplement the identities.
  **Pole values of J (type (a), C8f–g).** At Z = 0 the differential splits into the complex Z-plane
  and the two torus angles (δX = 0). At the south pole,
  J(δZ₁,δZ₂,δφ₂,δφ₃) = (−δZ₂,δZ₁,−δφ₃,δφ₂);
  at the north pole it is (−iδZ₂,−iδZ₁,δφ₃,−δφ₂).
  Indeed the Z-block is U diag(−i,i)U*, while w = φ₂+iφ₃ at 0 and w = −φ₂+iφ₃ at ∞.
  The slice corrections have zero first derivative at the node: t = x_s y_s is quadratic,
  α(0) = c̃(0) = 0, H(0) = 0, and ē₂,ē₃ have no linear term in t.
  **Complete tensor on these caps (identities (a), C9).** The following gives every coefficient of J,
  without differentiating F or solving for its inverse. It uses the explicit convergent series for H;
  it is a symbolic series expression, not a finite algebraic expression.
  Put σ=+1 at 0 and σ=−1 at ∞, and retain x_s,y_s,t and scale from C8. Set
  f=2|t|²/(1−|t|²), W=φ₂−iσφ₃. The inverse simplifies to
  A=tW+Hf, x=x_s e^{−A}, y=y_s e^A,
  **w=σ(W̄−Wt²+2Ht)**, since c̃₂=2 Re(t)/(1−|t|²) and c̃₃=2σ Im(t)/(1−|t|²).
  Thus the three complex one-forms
  η₁=dx_s−x_s dA, η₂=dy_s+y_s dA,
  η₃=σ[dW̄−t²dW+(2H−2Wt)dt+2t dH]
  equal e^A dx, e^{−A}dy, dw pulled back by G. All their coefficients are explicit:
  dt=y_s dx_s+x_s dy_s, df=4 Re(t̄ dt)/(1−|t|²)², and
  dA=Wdt+t dW+Hdf+f dH. For dH, put X̂=a₀|x_s|², Y=|y_s|²,
  T_j=e^{−j(j+1)f} X̂^{j(j+1)/2}Y^{j(j−1)/2}, D_h=Σ_j T_j, H=Σ_j jT_j/D_h. Then
  **dH=Σ_j(j−H)dT_j/D_h**, where, with a=j(j+1)/2 and b=j(j−1)/2,
  dT_j=e^{−j(j+1)f}[aX̂^{a−1}Y^b dX̂+bX̂^aY^{b−1}dY−j(j+1)X̂^aY^b df].
  Terms with a=0 or b=0 are omitted before evaluating a zero coordinate. Here
  dX̂=2a₀ Re(x̄_s dx_s), dY=2 Re(ȳ_s dy_s). The series and these derivatives converge uniformly
  on the caps: T_j≤rad^{j²}, rad≤0.004, and differentiation introduces only polynomial factors in j
  and lowers the monomial degree by at most one; f and df have bounded denominators there.
  In particular D_h≥T₀=1, and the formulas include both axes and the poles.

  Let r=(Re Z₁,Im Z₁,Re Z₂,Im Z₂,φ₂,φ₃), and let M be the 6×6 real coefficient matrix of
  (Re η₁,Im η₁,Re η₂,Im η₂,Re η₃,Im η₃) relative to dr. With J_std=diag(J₂,J₂,J₂),
  J₂=[[0,−1],[1,0]], the complete answer in these coordinates is
  **J_cap(r)=adj(M) J_std M / det(M)**, with
  **det(M)=−σ(1−|t|⁴)/(scale²a₀) ≠ 0**.
  To verify the determinant, first change φ to C=φ+Hc̃; this is a real shear of determinant 1.
  At fixed C, the complex normal block has determinant 1. Subtract
  (C₂ē₂′+C₃ē₃′)(y_sη₁+x_sη₂) from the w-row, leaving ē₂dC₂+ē₃dC₃, whose real determinant
  is σ(1−|t|⁴). The real determinant of Z↦(x_s,y_s) is −1/(scale²a₀).
  Multiplying dx,dy by e^A,e^{−A} changes the real determinant by |e^A e^{−A}|²=1. These factors
  are C9a,b,e,l; the bound |t|<0.016 from C8 proves nonvanishing. Hence J_cap²=−Id, and the
  coframe identifies its (1,0)-forms with differentials of local holomorphic coordinates: it is integrable.

  The matrix in the r-coordinates is independent of the torus angles. Indeed W↦W+Θ corresponds to
  (x,y,w)↦(e^{−Θt}x,e^{Θt}y,w+σ(Θ̄−Θt²)), a holomorphic change of chart for constant Θ (C9d).
  Thus one may set W=0 at evaluation while retaining dW=dφ₂−iσdφ₃. There is no angular branch choice.
  For the ambient ℝ⁹ tensor, put E=∂(Z,X,cos2πφ₂,sin2πφ₂,cos2πφ₃,sin2πφ₃)/∂r and R=dr∘P,
  with P the orthogonal tangent projector. Explicitly E has Z-block I₄, X-row −Z_ℝᵀ/X, and circle columns
  2π(−sin2πφ_j,cos2πφ_j); R has the first four rows of the S⁴ projector and circle rows
  (−sin2πφ_j,cos2πφ_j)/(2π). Then **J_ambient=E J_cap R**, RE=Id, and J_ambient²=−P.
  This expression is implemented separately as `jround.J_linear_ball`, with the existing finite height
  sum j=−8,…,8; no truncation error certificate is asserted here. Its value uses neither automatic
  differentiation nor Newton. `check_jround.check_symbolic_linear_balls` independently compares the
  formula with the forward-map evaluator, determinant, axes, boundary points and finite differences.
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
* **The linear balls are a sufficient fixed-point design** (notes §13.3.1). For the uncorrected chain, with
  Λ = M₀v_q + δ, δ = O(|q|⁴), lift B gives Z₁ = ω_E·√(4π|e|² + (v₃/|v|)δ₃ + …), not smooth across the arc, and
  lift A gives Z₂ = √(2π)(U₀q)₂·√(1 + δ-stuff/|(U₀q)₂|²), only C^{1,1} at p₀. Hence the linear blends.
  These expressions alone do not exclude smooth nonlinear choices: the divisibility of the remainder matters
  (L8). For example, with w = Uq and r = |w|², set Z = √scale·√(1+r²)w and
  X = ±√(1−scale·r(1+r²)). This is smooth near q = 0 and its orbit map is
  v = (1+|v_q|²)Mv_q, which is nonlinear. The Hopf identities are L9a–b (type (a)); the height radicand
  is > 9/10 for r ≤ 1/100 and scale ≤ 7 (type (c), L9c), covering both existing scales π and 2π.
  Thus exact linearity is a convenient sufficient choice, not a necessary smoothness condition.

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
**Corrected bound chain, type (c) (0 < r = rad ≤ R = 1/20, both cusps).** All constants below are rational;
`sym/node_remainder.py` I1a–i verifies their algebra and endpoint comparisons. The earlier constants 2.63 and
2.8 were not established by the displayed estimates; the complete bounds below replace them.
1. Put A = 55. Then 1/a₀ ≤ e⁴ < A: the positive Taylor series through k = 16, plus
   (4¹⁷/17!)/(1−4/18), is a rational upper bound for e⁴ and is < 55.
   Consequently |t|² = X̂Y/a₀ ≤ Ar²/4 ≤ 55/1600 < 1.
2. For k = ℓ+2 ≥ 2, k² ≥ 4+5ℓ, because the difference is ℓ(ℓ−1) ≥ 0.
   Since T_j ≤ r^{j²}, setting q = r⁵ gives
   R₂ ≤ 2r⁴/(1−q) ≤ (201/100)r⁴ and
   |R₁| ≤ 2r⁴(2−q)/(1−q)² ≤ (401/100)r⁴.
   The two rational functions increase for 0 ≤ q ≤ R⁵ (their positive power series show this), so their
   endpoint bounds suffice. Write d₀ = 201/100 and d₁ = 401/100.
3. From 4r³/27−X̂²Y = (X̂−2Y)²(4X̂+Y)/27 ≥ 0 and §4.1,
   2φX̂ ≤ Br³, where **B = (16A/27)/(1−AR²/4) = 281600/8343**.
   N2d, D ≥ 1, |m₀| ≤ r, and 0 ≤ 1−e^{−2φ} ≤ 2φ now give
   |μ̃−m₀| ≤ r²+Br³+d₁r⁴+d₀r⁵ ≤ Cr², with
   1+BR+d₁R²+d₀R³ < **C = 27/10**.
4. D ≤ 1+r+d₀r⁴ and e^{−2φ}X̂+Y ≥ r−Br³ > 0. Also |μ̃| ≤ r(1+Cr).
   Using 1/(1+u) ≥ 1−u for u ≥ 0, and discarding the positive cross term, gives
   Var/r ≥ (1−Br²)/(1+r+d₀r⁴) − r(1+Cr)²
   ≥ 1−r−d₀r⁴−Br²−r(1+Cr)² ≥ 1−4r,
   since **1+d₀R³+BR+(1+CR)² < 4**.
5. Finally |c′| ≤ 2/(3/100) by L1. Therefore
   **−a₁′/(2r) ≥ 1−4R−2CR²/(3/100) = 7/20 > 0**.
   This proves strict monotonicity along the scaling flow (x₀,y₀) ↦ (e^{−s}x₀,e^s y₀) on the punctured node
   ball, including its axes: N3c gives dT_j/ds = −2jT_j directly, without log coordinates.
   Where y₀ ≠ 0, s is the increment of ρ. At the node the flow is stationary; the local polynomial model applies.

Thus **I1 is certified with explicit constants**. For comparison, N5 (200000 random points per cusp,
seed 0, numerical evidence) gives max |μ̃−m₀|/r² = 0.9993 / 1.8713 against the common bound 2.7,
max (1−Var/r)/r = 1.998 / 2.479 against 4, and min −a₁′/(2r) = 0.8521 / 0.8261 against 0.35.

### 4.3 I2: the linear-blend shell of Λ_mod

Λ_mod = (1 − c(rad))M₀v_q + c(rad)Λ(v_q), c = s((rad − 0.004)/0.004). Identity (N4):
  **DΛ_mod = (1 − c)M₀ + c·DΛ + (Λ − M₀v_q)⊗∇c**, so
  σ_min(DΛ_mod) ≥ scale − c·|DΛ − M₀| − |∇c|·|Λ − M₀v_q|, |∇c| ≤ 2/0.004 = 500 (L1), scale = 2π (resp. π).
The following certificate proves |Λ − Mv_q| ≤ C_Λ|v_q|² and |DΛ − M| ≤ C′_Λ|v_q| on the whole ball
|v_q| ≤ 0.008. Its constants are **(C_Λ,C′_Λ) = (89,178) at 0 and (89/2,89) at ∞**. These are larger than
the earlier proposed, unproved targets (60,120) and (30,60), but give strict margins **σ_min > 2 and > 1**.
The identities and rational comparisons are `sym/i2_margins.py` I2a–n; the argument below is type (c).

1. **Reduce the composition (identities (a)).** Put z = 2πv_q = (a,b,c),
   d = √(1+a²+b²), p = a cos c+b sin c, D = d+p, and
   B(z) = (a sin c−b cos c, −sin c, D−cos c)/D.
   Write Qz = (−b,−c,a), I(a,b,c) = (a,−b,c), O(a,b,c) = (−b,−a,c).
   On rad = |v_q| ≤ 1/125, the node height is exactly v_{q,3}, and the fitted twist R is the identity:
   |t| ≤ 0.004 < 0.005 at 0, and |t| ≥ 250e⁻² > 26 at ∞. The local formulas are
   **Λ₀(v_q) = B(2πv_q)** and **Λ_∞(v_q) = ½ O B(2π I K(v_q))**, where
   K(v_q) = (βv_{q,1},βv_{q,2},v_{q,3}), ρ = √(v_{q,1}²+v_{q,2}²), and
   β = K_∞/k(2e⁻²/ρ), extended by 1 on ρ ≤ 2e⁻⁶.
   For the first formula the S³ coordinates are (p, a sin c−b cos c, cos c, sin c)/d.
   At ∞ they are (−cos c, sin c, −a cos c+b sin c, a sin c+b cos c)/d after inserting β.
   Shifted stereographic inversion has the identity
   w/|w|² = (P₁,P₂−1−P₀,P₃)/(2(1−P₂)) for P ∈ S³.
   This gives the second formula without a singular denominator at the node. I2a–f check these identities
   and DB(0) = Q. The radial zones used by P_final are justified in step 2.
2. **One derivative bound for B.** Put s = |z| ≤ S = 51/1000 and
   N = (a sin c−b cos c, −sin c, D−cos c), so B = N/D. The elementary bounds
   |sin c| ≤ |c|, 0 ≤ 1−cos c ≤ c²/2, |c−sin c| ≤ |c|³/6 follow by successive integration;
   d−1 = (a²+b²)/(d+1) ≤ (a²+b²)/2. Consequently
   |N−Qz| ≤ s²+s³/6, ‖∂N−Q‖ ≤ 2s+s²/2,
   |D−1| ≤ s+s²/2, D ≥ 1−s, and |∇D| ≤ 1+s+s²/2.
   To check the derivative estimate, split N−Qz into the rotated (a,b)-pair minus its value at c=0,
   the component c−sin c, and the component d−cos c. The first derivative has norm at most
   √(c²+a²+b²)=s (Cauchy–Schwarz); the latter two have bounds c²/2 and
   √((a²+b²)/d²+sin²c) ≤ s. I2f verifies the rotation and gradient identities used here.
   The quotient rule now gives
   ‖∂B−Q‖/s ≤ (2+s/2)/(1−s) + (1+s/2)/(1−s)
   + (1+s+s²/6)(1+s+s²/2)/(1−s)² ≤ 4.443503 < **9/2**,
   where the last comparison is an exact rational evaluation at S (every term is increasing).
   Also |B| ≤ (s+s²+s³/6)/(1−s) < 1/16 < e⁻¹. Thus at 0 the small-|w| branch of P_final
   applies; at ∞, |w|⁻¹ = |B|/2 < e⁻¹ gives the large-|w| branch. These inequalities also establish
   D > 0. The normalization is valid because 2π/125 < S.
3. **Control the infinity profile.** The rational bounds 157/50 < π < 22/7 follow from
   π = 16 arctan(1/5)−4 arctan(1/239) and alternating series through degrees 5 and 7 (I2g): the tangent
   of 4 arctan(1/5)−arctan(1/239) is 1 and that angle lies between 0 and π/2.
   I2h checks the positive Taylor bounds e < 3, e² < 8, e⁴ < 55, e⁶ < 500,
   e^{11/2} < 250, e^{224/15} > 3000000, and e^{13/5} > 88/7.
   On the nonconstant part of k, u = (log(2e⁻²/ρ)+4)/8 ≥ 15/16; put v = 1−u ∈ (0,1/16].
   With L = log(K₀/K_∞) = 2 log(4π)−2 ∈ (0,16/5), β = exp(−L s(v)), where s is the smooth step.
   Since s(v) ≤ exp(16/15−1/v) ≤ exp(−224/15) < 1/3000000 and
   s′(v) ≤ 258 exp(−224/15), we get
   **|1−β| < 1/500000**, **|dβ/d log ρ| < 1/25000**, hence
   ‖∂K−Id‖ ≤ δ = 21/500000 and |K(v_q)| ≤ |v_q|.
   For the derivative estimate, e^{−1/v}/v² is increasing on (0,1/16], with derivative
   e^{−1/v}(1−2v)/v⁴ > 0, and (1−v)⁻² < 2; multiply by e^{16/15}.
   Moreover K = Id when rad ≤ 1/250, since e⁶ < 500.
4. **Scale and blend.** At 0, ‖DΛ₀−M₀‖ ≤ 18π² rad < 178 rad. At ∞ the chain rule gives
   ‖DΛ_∞−M_∞‖ ≤ 9π² rad(1+δ)+πδ on the profile's support, which has rad > 1/250;
   therefore ‖DΛ_∞−M_∞‖/rad ≤ 9π²(1+δ)+250πδ < 89. Off that support the same bound holds
   with δ = 0. Integrating each derivative bound from 0 along a radial segment yields the respective
   quadratic remainders 89 rad² and (89/2) rad². Finally |∇c| ≤ 500 by L1, so throughout the shells
   σ_min(DΛ_mod,0) ≥ 2π−178/125−500·89/125² > **2**,
   σ_min(DΛ_mod,∞) ≥ π−89/125−500·(89/2)/125² > **1**.
   The same bounds hold on the inner balls, where the map is exactly linear. In fact Λ_mod is injective
   on each ball: after multiplication by Q₀ᵀ or Q_∞ᵀ, its derivative has symmetric part bounded below
   by 2Id or Id, respectively; integrate along the segment between two points. This certifies the local
   nonsingularity and uniqueness needed for I2, without assuming I3–I5 or modifying the construction.

**I2 margins (`sym/i2_margins.py`, 4000 points per shell, rad ∈ [0.004, 0.008]):**

| cusp | max \|Λ − Mv_q\|/\|v_q\|² | max \|DΛ − M\|₂/\|v_q\| | scale | bound from N4: scale − max\|DΛ − M\| − 500·max\|Λ − Mv_q\| | measured min σ_min(DΛ_mod) |
|---|---|---|---|---|---|
| 0 | 39.5 | 79.0 | 2π | **4.42** | 5.64 |
| ∞ | 19.7 | 39.5 | π | **2.20** | 2.82 |

These measurements remain numerical evidence; the certificate uses the larger bounds proved above.
The measured minima 5.6446 / 2.8197 exceed the proven lower bounds 2 / 1. The infinity shell is not
entirely on the constant-k plateau: at its outer equator |s| = 0.004e² and |t| = 250e⁻² < e⁴,
because Σ_{k=0}^7 6^k/k! > 250 (type (c), checked in `sym/i2_margins.py`). The same script evaluates
k/K_∞−1 ≈ 4.87·10⁻⁷ at this one point (numerical evidence). An infinity bound must include this profile
transition; it cannot be obtained by simply halving the cusp-0 constants. Step 3 includes this term.
**I2 is certified**, with the Jacobian identity N4 and the new local identities checked symbolically.

## 6. Status (2026-10-07) and resume checklist

### 6.1 The lemma list of the proof scheme (notes §13.1), with its status

This table concerns the original fitted-lift architecture of §§3–5. Section 8 gives a
separate global construction by quadratures, removing its dependence on lemma L.

| # | statement | status |
|---|---|---|
| A | Y's charts, gluings, lattice coordinates, Φ, Θ, ψ̂, T³-actions: identities | symbolic (the paper's `verify_construction.py`; here C1–C3, C6) |
| B | lattice system nonsingular: det = (\|t\|²+1)(log\|t\|(\|t\|²−1) + 2\|t\|²) > 0 | **manifest** (C1d) |
| C | μ̃ strictly monotone along fibres (dμ̃/dρ = −2Var), range ℝ, equivariant, closed on the central fibre | **manifest** (C4, C5) |
| D | λ-blend a₁ = (1−λ)μ̃ + λ(x₄^slice − ½) monotone, no wrap needed | **manifest** (§3.3; P < 0 by §4.1) |
| E | node correction monotone (I1) | **corrected rational certificate** (§4.2; −a₁′/(2rad) ≥ 7/20 at both cusps) |
| F | \|T\| strictly increasing (needs s′ ≤ 2) | **proven** (L1, L4) |
| G | ν strictly increasing (P_final) | **manifest** (L3) |
| H | R a diffeomorphism (fibre twist, det 1), moebius, P_final zones: closed inverses | symbolic (L5, L6) |
| I | exact linear models M₀ = 2πQ₀, M_∞ = πQ_∞ | symbolic (`symbolic_identities.py`) |
| J | linear-blend shells nonsingular (I2) | **certified** (§4.3, I2a–n): derivative constants 178 / 89, singular-value margins > 2 / > 1, including the infinity profile |
| K | lifts land in S⁴, cover Λ, equivariant; linear in the balls | symbolic (L7, §3.4) |
| L | E ≠ 0 off the tube (I3); G_A ≠ 0 on the tube (I4); \|Arg g\| < π on the shell (I5) | **open global inequalities** (§5); single-lift alternative additionally requires smooth division (L8), not just agreement of first derivatives |
| M | φ: Y → T² equivariant, invariant, closed inverse | symbolic (C6, §3.2) |
| N | F∘G = id on each piece; pieces of G agree on overlaps | **C8 verifies the closed linear-ball pieces; C9 gives the full tensor there**; other pieces and overlaps remain to be assembled, subject to L; numerical inversion outside the caps has a central-fibre limitation (T14) |
| O | DG invertible ⇒ J = G*J_Y smooth, J² = −1, N = 0 ⇒ Newlander–Nirenberg | **C10 verifies the generic six-dimensional N=0 identity**; its local-diffeomorphism hypothesis remains essential |

**Generic integrability certificate (type (a), C10; T9).** Let F be any smooth local diffeomorphism
between real six-dimensional coordinate domains, A=DF, and J=A J_std A⁻¹. Then J²=−Id by matrix
multiplication. At any chosen point a constant real linear change of target coordinates makes A=Id;
this leaves an arbitrary symmetric second jet H^a_{ij}=H^a_{ji}. Differentiating gives
∂_k J=H_k J_std−J_std H_k, where (H_k)^a_i=H^a_{ik}.
Substitute into
N^a_{ij}=Σ_p[J^p_i ∂_p J^a_j−J^p_j ∂_p J^a_i−J^a_p(∂_i J^p_j−∂_j J^p_i)].
`sym/chart_identities.py` C10 expands all 216 components with 126 independent Hessian entries and
verifies that every component vanishes identically, as does the derivative of J²+Id. Since N is a
tensor and the point was arbitrary, this proves the identity for general A. Equivalently, the frame
F_*∂_i commutes, and J has constant coefficients in that frame. Applying this to F=G⁻¹ proves
integrability for a smooth local inverse G, including the cap formula of C9. It does not establish
the smoothness, existence or global injectivity of a proposed F or G; I3–I5 and the remaining
domain/overlap checks are separate obligations. A control with A=Id+q₂E₀₀ (an invertible frame matrix
near 0, not a coordinate Jacobian) gives N¹_{0,3}=−1 at 0, although J²=−Id (C10c, zero-based indices).

### 6.2 What this session established

* The symbolic specification of F and G is complete at the level of blocks (§3): every block is closed-form or a
  certified-implicit one-variable inversion, except the two Λ_mod shells (3D implicit, certificate J) and the
  lift's dependence on E, G_A (L).
* Three lemmas of notes §13.3 are now proven or manifest (F, G, D); one is unnecessary (the wrap bound).
* New exact facts: φ(t) = 2|t|²/(1−|t|²) at both cusps; the lattice determinant; the node-variable form of the
  theta series with degree j² terms; the two-term truncation accurate to O(rad⁴).
* I1 and I2 now have explicit rational certificates (§4.2–4.3); I2 uses a common elementary local formula
  and includes the infinity profile transition. Measured margins remain separately labelled.
* The honest status of the global conditions: **I3–I5 remain inequalities over a continuum**; a purely symbolic
  proof of the whole theorem is not available with this architecture, and §5 explains why. They are certifiable by
  validated numerics (notes §12.12) or by analytic bounds in the structural regions.

### 6.3 Not done / next steps

1. I2 is now certified (§4.3); use the reduced local formula when assembling the shell inverse and its domains.
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
* **Architecture (C), single compatible lift — an open design:** one weight-1 function E with a
  manifest zero set and a smooth nonvanishing quotient (v₁+iv₂)/E near the whole arc (L8);
  matching linearisations alone is insufficient. With that stronger condition, lift B alone is smooth off the
  poles (it is smooth across Γ′ since its only denominator vanishes on Γ_E), the gauge, ψ, the shell, `fit_psi`
  disappear, and `fit_arc` disappears because the arc is explicit. Ingredients: (i) in the middle, E =
  e^{2πix₁}·ϑ with an r-dependent modulus τ(r) ∈ ℍ and explicit zero (one simple zero per cell, classical);
  (ii) at the cusps, the single-term local model y₀ + x̄₀ near the arc (exact for |t| ≤ 0.3) with zero at
  {|x₀| = |y₀|, arg t = π} and the rigid induced structure z_E = log|y₀| + (i/2)arg t (independent of the model's
  coefficients α(t), β(t)); (iii) gluing cusp ↔ middle by a log-ratio on an annulus where their quotient is smooth
  and nonvanishing (matching τ(r) and the zero alone does not establish this), with the ratio's windings computed from the explicit
  local model; (iv) compatibility of Λ with E along the arc by an explicit anti-conformal twist
  (v₁+iv₂) ↦ A(λ)(v₁+iv₂) + B(λ)·conj(v₁+iv₂) in the w-chart (closed-form A, B from derivatives of closed-form
  maps; |B| < |A| certifies the linear normal map, but smooth division still needs a separate argument by L8). Residual inequalities:
  I1, I2 (explicit constants, certified), I3-cusp = E ≠ 0 off the arc on |t| ≤ 0.3 (the aligned-phase
  argument; reducible to explicit-constant series bounds by narrowing the χ-supports so that the overlaps have
  |y₀|² manifestly > |t|), and the 1D sign condition of (iv). Smooth division and a global extension of the
  local twist are additional open requirements; the single-lift construction is not yet certified.
* Fallback for any residual inequality that resists: validated numerics (notes §12.12), stated as such.

**What the symbolic J would be:** per piece, J(q) = DF(G(q))·J₀·DG(q), with the implicit 1D inverses as symbols
carrying their defining identities and IFT derivative rules; checks: F∘G = id per piece, DF(G)·DG = 1, overlap
agreement through the chart transitions, domain covering (manifest unions), the generic naturality lemma, and the
lemma list of §6.1 with explicit constants.

## 8. Global construction by exact quadratures

This section specifies a different diffeomorphism **G_Q: S⁴×T² → Y**, using the same Y, heights,
node corrections and T²-trivialization. It replaces the equivariant lift, and takes the optional twist
R to be the identity. Neither E, G_A, an arc fit, a phase fit nor I3–I5 occurs in its definition.
The resulting tensor is **J_Q = G_Q*J_Y**; it need not equal the earlier J_round pointwise, since
the parametrization has changed. They describe the same complex manifold if the earlier map is
also a diffeomorphism. The formulas here do not assume that unresolved assertion.

**Expression class.** In addition to elementary functions, the specified convergent height series
and certified scalar inverses, the expression uses ordinary definite integrals of explicitly given
smooth real functions of one variable, and one contraction limit with ratio <17/25. These are
exact quadratures, not finite elementary closed forms. The formulas below specify the integrands,
paths, chart changes, limits and inverse; there is no unspecified ODE solution or arbitrary bundle
isomorphism. Numerical truncations are separate. All algebraic identities used below are checked
in C1–C15, L1–L9 and I2a–p; signs use manifest positivity (b) or the explicit bounds (c) in §4.
The differentiation and gluing arguments are given explicitly. This extension of the expression
class, and the global conclusion, require collaborator review under AGENTS.md §3.1.2.

### 8.1 The base map and a globally specified inverse

Write ℋ(z₁,z₂)=(2 Re(z₁z̄₂),2 Im(z₁z̄₂),|z₁|²−|z₂|²). Let
Λ_* = P_final ∘ moebius ∘ h_B, with every block as in §3.3, **R=Id**, and put
Λ_Q equal to its node-linear modification. The two exact local formulas for the v-component
of Λ_* are B₀(p)=B(2πp) and B_∞(p)=O B(2π I K(p))/2 from §4.3. Set

  R_b=1/125, r_b=1/250, c_b(r)=s((r−r_b)/(R_b−r_b)),
  f_j(p)=M_j p+c_b(|p|)[B_j(p)−M_j p],  |p|≤R_b.

Here p=ℋ(√a₀ x̄_s,y_s); X=−√(1−|f₀(p)|) or +√(1−|f_∞(p)|).
Outside these balls use Λ_*. The local expressions agree to every order at |p|=R_b,
because the smooth step is identically 1 on the outer side. At p=0 use the linear expression.
The integer choices in the node correction introduce no seam: on rad≤1/20,
I1 gives |μ̃|≤1/20+(27/10)(1/20)²=227/4000<1/2. Thus only one nearest-height
representative can support the correction. Where the nearest integer changes, the
correction is identically zero on a neighbourhood. The smaller base-blend balls
satisfy the same condition.

**The raw inverse, including the central fibres.** L3–L6 and §3.3 specify the inverse away
from u₁u₂=0: recover w, u∈S³, |t|, t and the circle-valued height, then invert the strictly
monotone fibre coordinate. On u₁=0 set t=0 and a=arg(u₂)/(2π) mod 1 in [−1/2,1/2];
on u₂=0 use the infinity chart and a=−arg(u₁)/(2π)+1/2 mod 1. Set

  d(r)=c_n(r)r+(1−c_n(r))r/(1+r),
  c_n(r)=1−s((r−1/50)/(3/100)),   r=d⁻¹(|a|), 0≤r≤1.

For a>0 take x_s=√(r/a₀), y_s=0; for a<0 take x_s=0, y_s=√r; for a=0 take
the node. The two choices a=±1/2 represent the same orbit by ψ̂, with **d(1)=1/2**.
Strict monotonicity of d on its blend
is I1 applied on the axes; elsewhere d′=1 or 1/(1+r)²>0. Thus no polar division is
needed on a central fibre. The radius is closed outside the blend: r=|a| for |a|≤1/50,
and r=|a|/(1−|a|) whenever the resulting r≥1/20.

For t≠0 a convenient certified bracket for the height root is
ρ∈[(a+5/2)P,(a−3/2)P], P<0. C7 gives |μ̃−(ρ/P−1/2)|<1/2, and I1 gives
|μ̃_new−μ̃|≤(27/10)(1/20)²=27/4000 wherever the correction is active.
The blended coordinate therefore differs from ρ/P−1/2 by less than 2027/4000<1.
Its derivative is strictly negative, so the bracket contains exactly one root. Bisection specifies
the inverse without a Newton initial guess. Local smoothness follows by differentiating this
nonzero derivative. At a nonnodal central point, (Re t,Im t,height) are local orbit coordinates:
dt has rank two and the scaling flow fixes t while changing the height strictly (I1/C5).
At the nodes the coordinates p and the formulas B_j supply the regular description.
The clutching identity a₂′−a₁=arg(t)/(2π) mod 1 (C15a), together with the positive derivative of |T|,
proves that h_B is one-to-one and onto S³, including its two coordinate circles. Hence this
raw inverse has no omitted orbit or duplicate orbit. No degree measurement is used.

**The shell inverse is a contraction (c).** Let Π(p)=p min(1,R_b/|p|), Π(0)=0, and
Δ_j=f_j−M_j on the ball. The bounds in §4.3 give

  Lip(Δ₀)≤178R_b+500·89R_b²=534/125,
  Lip(Δ_∞)≤89R_b+500·(89/2)R_b²=267/125.

Projection Π is nonexpanding: C14d–e express the difference of squared distances as
(r−R_b)(r+R_b−2s z)≥0 when r≥R_b≥s and |z|≤1, or
(r−s)²+2(rs−R_b²)(1−z)≥0 when r,s≥R_b. Thus for every v∈ℝ³ the iteration

  p₀=M_j⁻¹v,   p_{k+1}=M_j⁻¹[v−Δ_j(Π(p_k))]                         (8.1)

has contraction ratio κ<17/25 at either cusp (I2o, using the rational Machin bounds).
Indeed successive distances decrease by κ; summing their geometric tail gives convergence
and |p−p_k|≤(17/25)^k |p₁−p₀|/(1−17/25). Its limit is the unique solution of
v=M_jp+Δ_j(Πp). If |p|≤R_b, this is the required shell/ball inverse. If |p|>R_b,
use the raw global inverse above. Choose j=0 for X<0 and j=∞ for X>0; at X=0
neither ball is possible and use the raw inverse.

**Why this domain test is exact.** Apply the same argument to δ_j=B_j−M_j;
its Lipschitz bounds are 178R_b and 89R_b, smaller than those just used. The two
extended maps M_jp+Δ_j(Πp) and M_jp+δ_j(Πp) are bijections of ℝ³ by their
convergent inverses, and coincide outside the ball, including its boundary. Therefore their
images of the ball coincide: an exterior inverse for either is also the exterior inverse for
the other, uniquely. Replacing the raw map on a ball changes neither its image nor its
boundary. The two ball images are disjoint, one in each sign of X; their |v| is <1/16
by §4.3 and convexity. Consequently Λ_Q is globally bijective. Its inverse is smooth off
the two poles by I1, L3–L4 and I2; on each pole ball its v-coordinate is exactly linear.
This supplies both a global injectivity argument and an explicit inverse-domain test.

`tlift.transport_base_local/chart/middle` and `tinverse.transport_base_inverse` evaluate this
new base independently. The latter uses 96 contraction steps and 64 bisections for numerical
checks; (8.1) and the infinite height series define the exact map. The direct stereographic
formula (P₁,P₂,P₃)/(1+P₀)−(0,1,0) is used away from P₀=−1 (the infinity node).
It also covers P₀=+1, where the legacy rationalization has a removable 0/0.

### 8.2 A connection with a positive normalization

The circle action is free away from the two nodes: a nonzero deck power changes
H by a nonzero integer, and on a fixed representative either x_s or y_s is nonzero,
so (λ⁻¹x_s,λy_s)=(x_s,y_s) forces λ=1. The middle coordinate action has the same
property (C2–C3). Thus the phase charts below are principal circle charts, with no
finite stabilizer omitted.

On a cusp slice away from the node, define the real one-form

  α_loc=Im(−a₀ x̄_s dx_s+ȳ_s dy_s)/(a₀|x_s|²+|y_s|²),
  S(H)=Σ_{n∈ℤ} χ(H+n),
  α_j=Σ_{n∈ℤ} χ(H+n)(ψ̂ⁿ)*α_loc / S(H),                       (8.2)

where H=μ̃ is the **uncorrected** height and χ(h)=s((9/10−|h|)/(1/2)). The sum is
locally finite, with at most two positive terms. For an integer n with |H+n|≤1/2,
χ(H+n)≥s(4/5)>1/2, hence S(H)>1/2 (b; C11e). Every local denominator is a sum
of squares and is positive off the nodes. The circle generator is K=(−ix_s,iy_s),
so α_loc(K)=1 (C11a); each summand is circle invariant. Reindexing n under H→H+1
shows that (8.2) descends through ψ̂. On the slab |H|<3/5 the only possible indices
are −1,0,1; the n=1 term requires H<−1/10 and therefore y_s≠0, and the n=−1
term requires H>1/10 and therefore x_s≠0. Indeed on y_s=0 the height is nonnegative,
and on x_s=0 it is nonpositive (C4). The explicit ψ̂ and ψ̂⁻¹ of C3 therefore have
no zero denominator on their coefficient's support, including the central axes.

An alternative expression on the torus is especially short. Put ρ=log|ζ|,
P=log|m(t)t|, ρ_n=ρ+nP, and

  b_n=a₀|t|² e^(−2ρ_n)/(a₀|t|² e^(−2ρ_n)+e^(2ρ_n)),
  A_j=Σχ(H+n)(n−b_n)/S(H),
  α_j=d arg ζ+A_j d arg t.                                      (8.3)

Here m(t)=ν exp[−2|t|²/(1−|t|²)] is positive real at both cusps (C9a).
The argument of ψ̂ⁿζ is arg ζ+n arg t, giving (8.3) from (8.2) (C11b).
Under ψ̂, A_j→A_j−1 and d arg ζ→d arg ζ+d arg t, so the form is unchanged
(C11c). Near a pole, rad≤1/250 gives |H|≤rad+(27/10)rad²<1/10 by I1.
Only n=0 contributes there, with χ(H)=1. Thus α_j=α_loc **exactly** throughout
both linear balls. Also |t_c|≤e²/500<8/500<3/10, so the global cutoff (8.4)
equals its cusp value throughout these balls. No estimate of a complex phase is needed.

On the middle part set α_m=2π dx₁+arg(t) dx₄. The branch change
(arg t,x₁)→(arg t+2π,x₁−x₄) preserves it (C11d), as do x₄→x₄+1 and x₁→x₁+1.
Define β₀=1−s((|t|−3/10)/(1/5)) and β_∞=1−s((|1/t|−3/10)/(1/5)), and set

  α=β₀α₀+β_∞α_∞+(1−β₀−β_∞)α_m.                              (8.4)

The two cusp supports are disjoint (|t|<1/2 and |t|>2). All weights are nonnegative
and sum to one. Thus α is a smooth circle-invariant connection with α(K)=1 on
W minus its two fixed points. It is exactly the Hopf connection on the linear balls:
unitarity and constant scaling in the linear models preserve Im(Z* dZ)/|Z|².
Smoothness of the height series follows locally from its explicit terms: on |xy|≤r₀<1
and |x|²,|y|²≤M, its positive/negative tails are bounded by
M^k r₀^{k(k−1)}, k≥1, since g≤0. Derivatives on a smaller coordinate neighbourhood
have the same quadratic exponential decay times a polynomial in k and fixed powers of
M and r₀⁻¹. This justifies termwise differentiation and includes either axis by the
nonnegative exponents in the chart series. Cutoffs are constant near their absolute-value
origins, so those absolute values introduce no nonsmooth point.

### 8.3 Transport is a scalar integral with explicit chart changes

Use these local phase coordinates and references on the free part of W:

| chart | phase θ | reference at the same orbit |
|---|---|---|
| middle, a branch of log t and a lift of x₄ | 2πx₁ | x₁=0 |
| cusp, y_s≠0, on a chosen height slab | arg y_s | y_s>0, x_s=t/y_s |
| cusp, x_s≠0, on a chosen height slab | −arg x_s | x_s>0, y_s=t/x_s |

The reference's modulus is the root of §8.1. These charts cover: t∈ℂ* uses a middle
chart; at either central fibre every nonnode has x_s≠0 or y_s≠0; height shifts by
integers move it into |H|<3/5. In a phase chart a write α=dθ_a+A_a, where
**A_a=r_a*α**, the pullback along its explicit reference r_a. Its coefficients are
obtained directly from (8.2)–(8.4), or by differentiating the reference root. For a
root f(b,r)=0, the rule is dr=−(∂_r f)⁻¹d_b f; its denominator was certified in §8.1.

On an overlap write θ_b=θ_a+h_ba(base), with h defined modulo 2π. All transitions
are explicit: from the y-phase to the x-phase, h=−arg t; under ψ̂ⁿ in the y-phase,
h=n arg t; for a log branch change in the middle, h=−2πx₄. For middle-to-cusp
transitions, use §3.2's Φ⁻¹ or Θ, then take the T²-slice and read arg y_s (or −arg x_s).
In particular, for a middle reference z=−i x₄ log t, w=2tx₄, set at cusp 0
z_c=z−πs(s−1), w_c=w−2ts, s=log t/(2πi); at ∞ set z_c=z−iw/t,
w_c=w/t², t_c=1/t. Solve b₂ē₂+b₃ē₃=−w_c, and put
ζ_s=exp(i z_c+b₂α₂+b₃α₃). Its argument is the y-phase transition.
These are closed formulas (C1–C3); no lift E or G_A enters them.

For a base path γ on a segment [a,b] contained in chart a, transport is

  θ_a(b)=θ_a(a)−∫_a^b A_a(γ(s))[γ′(s)] ds.                     (8.5)

Across an endpoint add h_ba and continue. A compact path admits a finite partition
into such segments; any such partition gives the same point. To check this, equality
of α gives A_b=A_a−dh_ba, so the two integrals differ by
−h_ba(γ(b))+h_ba(γ(a)), exactly cancelled by the endpoint changes (C11f).
Refinement, concatenation and reversal hold by addition and negation of integrals
(C11g). Thus (8.5) defines an explicit invertible map T_γ of circle fibres, commuting
with their circle action. Smooth dependence follows locally by choosing one finite
partition valid in a neighbourhood of a given path and differentiating its smooth
integrands on compact intervals. There is no branch requirement on the size of an
integral. Changing a real argument lift by 2π changes only a phase representative.

### 8.4 The north mismatch has an explicit global real logarithm

Fix ε=1/100 and x_−=−√(1−ε), x_+=√(1−ε). I2p gives ε<π/250, so both
levels are strictly inside their respective linear balls. For n∈S² define the base path

  γ_n(s)=Λ_Q⁻¹((1−s²)n,s),  −1<s<1.

Write L₀ and L_∞ for the **slice** linear maps of §3.4 (omit the T²-coordinate).
For (Z,X)∈S⁴ with Z≠0 set r=|Z|²=1−X², n=ℋ(Z)/r, Z_−=√(ε/r)Z, and

  G_H(Z,X)=T_{γ_n|[x_−,X]}(L₀⁻¹(Z_−,x_−)).                   (8.6)

When X<x_− the interval is oriented backwards. Formula (8.5) makes (8.6) circle
equivariant and smooth on the free part. In a linear cap, γ_n is radial in the local
Hopf coordinates and α_loc vanishes on radial vectors. Thus transport in that cap
is precisely positive radial scaling, an identity, not an asymptotic assertion.

At X=x_+ define g:S²→S¹ by

  g(n)=Z₀* pr_Z L_∞(G_H(Z₀,x_+))/ε,   ℋ(Z₀)=εn, |Z₀|²=ε.       (8.7)

Both vectors have the same Hopf image and length, hence differ by a unique unit
phase: C15b writes their rank-one projectors entirely in terms of their Hopf data,
so their complex lines coincide; equality of lengths gives modulus one. The expression
is independent of Z₀ because (8.6) is equivariant. Explicit
references Z₀=√ε z(n) are

  z(n)=(√((1+n₃)/2),(n₁−in₂)/√(2(1+n₃)))       if n₃≥0,
  z(n)=((n₁+in₂)/√(2(1−n₃)),√((1−n₃)/2))       if n₃≤0.

Their denominators are at least √2 in their stated regions; their Hopf identities
are L7/C14a. The two definitions of g agree on overlap by equivariance, so g is
globally smooth. Define η=Im(ḡ dg). Writing g=u+iv gives η=u dv−v du and
u du+v dv=0. Consequently dη=2du∧dv=0 (C12f), and dg=i g η (C13c).

Here is a real logarithm without an unproved Arg bound. For a∈ℝ² set

  S_N(a)=(2a₁,2a₂,1−|a|²)/(1+|a|²),
  S_S(a)=(2a₁,2a₂,|a|²−1)/(1+|a|²),
  P_N(a)=∫₀¹ (S_N*η)_{ta}[a] dt,   P_S(a)=∫₀¹ (S_S*η)_{ta}[a] dt.

Choose θ₀=Arg(g(0,0,1)) once in (−π,π], and e=(1,0). On the two punctured spheres define

  Θ(S_N(a))=θ₀+P_N(a),
  Θ(S_S(a))=θ₀+P_N(e)−P_S(e)+P_S(a).                           (8.8)

These are mutually compatible real functions. Indeed for a closed one-form A,
∂_j∫₀¹ A_i(ta)a_i dt=∫₀¹ ∂_t[t A_j(ta)]dt=A_j(a); C13b checks the derivative
identity for an arbitrary symmetric first jet. Therefore both differentials in (8.8)
equal η. Their difference is constant on the overlap S² minus its two poles, which
is connected (latitude and longitude paths connect any two points), and is zero at
(1,0,0). This proves agreement and smoothness, including both poles using the other
chart. Finally d(e^(−iΘ)g)=0 and its value at the north pole is 1, so

  **g=e^(iΘ) everywhere on S².**                                (8.9)

No winding computation, small-phase estimate or logarithm of a vanishing section
is hidden in (8.8). To evaluate it one may use S_N for n₃≥0 and S_S for n₃≤0;
then |a|≤1, keeping the quadrature domains bounded.

### 8.5 The diffeomorphism, its inverse and every domain

Set χ_X=s((X−x_−)/(x_+−x_−)). On the free part define

  **G_W(Z,X)=e^(−iχ_X Θ(n)) · G_H(Z,X).**                       (8.10)

For X<x_− this is exactly L₀⁻¹, by radial transport and χ_X=0. For X>x_+
radial transport preserves the mismatch g(n), so (8.9) gives exactly L_∞⁻¹.
Define G_W at the two poles by these same linear maps. Thus the three open domains
{X<x_−}, {−1<X<1}, {X>x_+} cover S⁴, and their formulas agree on their overlaps.
No n or angle is evaluated at Z=0. This proves smoothness at both fixed points.

The inverse on W minus its fixed points is equally explicit. Given w, compute
(v,X)=Λ_Q([w]), r=1−X² and n=v/r, and set

  w_−=T_{γ_n|[X,x_−]}(e^(iχ_X Θ(n))·w),
  (Z_−,x_−)=L₀(w_−),
  **F_W(w)=(√(r/ε)Z_−,X).**                                   (8.11)

Use L₀ or L_∞ directly at the fixed points and in the matching caps. Equivariance,
transport reversal, cancellation of the two phases and cancellation of the positive
radial factors give **F_W G_W=id and G_W F_W=id**, including caps (C8, C11g–h,
C14b–c). Smoothness follows from the same local integral formulas. In local base
and phase coordinates G_W has the form (b,θ)↦(Λ_Q⁻¹b,θ+f(b)), with Jacobian
det DΛ_Q⁻¹≠0. There is no fibre monotonicity or phase-size condition left to prove.

Restore the T²-coordinate using the existing φ. Given (w,φ)∈W×T², its unique
point in Y is specified as follows:

* On a cusp core |t_c|<3/10, choose a representative (x_s,y_s) of w, put
  H=μ̃(x_s,y_s), C=φ+H c̃(t_c), A=C·α(t_c), and take
  (x,y,w_c)=(x_s e^(−A),y_s e^A,C·ē(t_c)).
* In the middle, take a branch representation (t,x₁,x₄) of w and put
  (x₂,x₃)=φ+β₀D₀+β_∞D_∞, then (z,w_c)=Σx_j e_j(t), ζ=e^(iz), using §3.1.

The middle expression is used wherever t≠0,∞, so these domains cover Y.
Here D₀=ξ₂₃+(H−x₄+ξ₄)c̃₀ and D_∞=(H+x₄)c̃_∞ on their annuli.
C12b–e verify the corrections and their invariance under both deck and log-branch
changes. On the cusp cores the middle formula gives exactly C=φ+Hc̃;
the source coordinates are related by Φ or Θ (C1–C3). Adding an integer to φ
is the holomorphic lattice action, so the resulting Y-point is independent of the
angle representatives. Thus these formulas define the inverse τ⁻¹ of
τ:Y→W×T², τ(y)=([y],φ(y)). Define

  **G_Q(Z,X,φ)=τ⁻¹(G_W(Z,X),φ),**
  **F_Q(y)=(F_W([y]),φ(y)).**                                  (8.12)

The two compositions are identities, using (8.11) and C6/C12. Locally on the free
part, in the bundle coordinates (base, circle phase, φ), the Jacobian is block
triangular with diagonal blocks DΛ_Q⁻¹,1,I₂; C11i
checks its determinant. On either cap its nonzero determinant is the one in C9.
Hence these are smooth mutually inverse maps on the entire stated manifolds.

### 8.6 The symbolic tensor on the round embedded manifold

The following is a full coefficient prescription, including the derivative rules;
it does not leave the choice of a diffeomorphism to an existence theorem. For
q=(q₀,…,q₄,e^(2πiφ₂),e^(2πiφ₃)) choose k with |q_k|>1/√6 and a sign for q_k.
Use coordinates u=(the other four q_j,φ₂,φ₃), so
q_k=±√(1−Σ_{j≠k}q_j²). Such patches cover S⁴ because max|q_k|≥1/√5>1/√6;
two ordinary angle charts on each circle complete the cover of S⁴×T².
Evaluate (8.1)–(8.12) in any source holomorphic chart, and form a complex coframe η:

  middle: η=(dz,dw_c,dt),  dz=Σ e_j^z dx_j+Σ x_j d e_j^z,
                            dw_c=Σ e_j^w dx_j+Σ x_j d e_j^w;
  cusp:   η=(dx_s−x_s dA, dy_s+y_s dA, d(C·ē)),
          dA=Σ α_j dC_j+Σ C_j dα_j,
          dC=dφ+c̃ dH+H dc̃,  dt_c=y_s dx_s+x_s dy_s.             (8.13)

The cusp coframe is (e^A dx,e^(−A)dy,dw_c), a nonsingular complex multiple of the
holomorphic coframe, so it gives the same tensor. The explicit dH is the quotient
derivative of the series in C9f–g, with zero-exponent derivative terms omitted on
the axes. The elementary lattice derivatives come from §3.1. Derivatives of scalar
roots use the rule in §8.3; derivatives of the shell root use
dp=(Df_j(p))⁻¹dv, with I2 certifying its denominator. Derivatives of an integral use
the endpoint rule and differentiation under its integral sign. Repeating these
rules specifies all orders; one never differentiates a discontinuous chart selector.
Where the selector changes, use either overlapping smooth chart formula.

Let M_ab be the coefficient matrix in
(Re η₁,Im η₁,Re η₂,Im η₂,Re η₃,Im η₃)^T=M du,
and let J₀=diag([[0,−1],[1,0]],[[0,−1],[1,0]],[[0,−1],[1,0]]). The answer in
these six real coordinates is

  **J_Q(u)=adj(M(u)) J₀ M(u) / det M(u).**                      (8.14)

All entries of M are specified by (8.13) and the rules above; det M≠0 by (8.12)
and the nonzero coframe factors. For the ambient 9×9 matrix, let E=D_u q and
R=du∘P, where P is the orthogonal tangent projector, and use

  **J_Q^ambient(q)=E adj(M) J₀ M R / det M.**                   (8.15)

Explicitly the sphere part of E inserts the four coordinates and has the missing
row −q_{j≠k}^T/q_k; its circle columns are 2π(−sin2πφ_j,cos2πφ_j).
R extracts the other four sphere coordinates after I₅−q q^T and has circle rows
(−sin2πφ_j,cos2πφ_j)/(2π). Therefore RE=I₆ and ER=P (C9i–k, permuting the
missing sphere coordinate). It follows algebraically that (J_Q^ambient)²=−P,
J_Q^ambient P=P J_Q^ambient=J_Q^ambient, and J_Q²=−Id on the tangent bundle.

On overlaps, G_Q gives the same Y-point by §§8.3–8.5, and the holomorphic source
transition has a complex-linear differential commuting with J₀. Thus (8.14)–(8.15)
are chart independent. They are smooth at poles by exact agreement with C9 on
|Z|²≤ε, across
central fibres by the toric charts and (8.2), across shell boundaries by the flat
cutoff, and across torus/log branches by C2/C12. Pulling back Y's holomorphic charts
by the proven diffeomorphism G_Q is directly a complex atlas on the whole round
S⁴×T². Equivalently, C10 checks all 216 Nijenhuis components for an arbitrary
symmetric second jet of its coordinate map; hence **N_{J_Q}=0 identically**.
Here C10 is applied to the coordinate Jacobian of G_Q. The coframe matrix M in
(8.13) may differ by a complex-linear factor, which commutes with J₀ and cancels
from (8.14); it is not being assumed to be an arbitrary integrable frame.
The construction is T³-equivariant, and the source actions are holomorphic, so the
tensor is invariant under the prescribed circle and torus actions.
In the six angular target coordinates, T²-translations have identity differential.
Consequently M may be evaluated at φ=0 after retaining its dφ columns, as in C9;
the ambient matrices E and R still use the actual circle angles.

### 8.7 Verification scope

The global assertion in §§8.1–8.6 is an analytic construction by exact quadratures.
Its quantitative inputs are precisely I1, I2 and L1–L4, with their recorded margins;
the new normalization is >1/2 and the new inverse contracts by <17/25. The remaining
steps are the displayed identities and the elementary calculus/gluing arguments.
There is no inference from a small sampled Nijenhuis tensor to a global theorem.
The source input is Y's specified holomorphic atlas and lattice actions, checked in
the manuscript's `verify_construction.py` and C1–C3. No preprint-dependent sphere
classification input or nonconstructive diffeomorphism theorem is used in §8.

The implementation additions evaluate the new base, its inverse and the connection;
they **do not replace `Finv`, `Psi_W` or `J_round` by a complete numerical quadrature
evaluator**. Thus their old numerical results must not be quoted as tests of J_Q.
The exact tensor (8.14) is specified by quadratures and derivative rules, while C9
also supplies a direct numerical evaluator for J_Q on the reduced caps |Z|²≤ε.
Outside those reduced caps it evaluates the earlier construction's tensor.
The recorded numerical
checks for the new blocks, their sample sizes and errors are appended below and in
the ledger. A finite elementary expression without quadratures remains unestablished.

**New numerical observations (not proofs).** `check_inverse.check_transport_base` uses 2000
uniform S⁴ targets, 1024 local targets per cusp, and 256 axis targets per cusp (4560 total):
maximum shell-root error 3.6·10⁻¹⁷; forward/inverse residual 1.1·10⁻¹⁴. `check_tlift`
uses 4000 torus points and 4000 chart comparisons per cusp: connection deck error
≤1.4·10⁻¹⁵; chart agreement ≤2.2·10⁻¹⁶. It also checks 128 exact axis tangents per
cusp, 512 middle/cusp tangent comparisons per cusp, 3000 middle branch/deck comparisons,
and four smooth paths crossing central axes. For the latter, 96-node Gauss–Legendre
quadrature gives composition/reversal/gauge residual ≤4.4·10⁻¹⁶. An independent
test of (8.8) on 256 sphere points with a known phase whose range exceeds 2π has
unit-phase error ≤9.0·10⁻¹⁵, includes both poles, and checks agreement of the real
primitives on 64 equator overlaps to 1.2·10⁻¹⁴. These test the new blocks and quadrature identities;
they are not an evaluation of the complete G_Q or J_Q at global random points.
The base check also compares all 16 second-derivative directions at each exact pole
against the real quadratic Hopf model. The numerical base uses the inner polynomial
radius and an inactive-cutoff guard there; differentiating a norm at zero is avoided.

## 9. Closed-form reductions and the limits of the present result

### 9.1 Expression class and audit

Here a *finite elementary formula* permits arithmetic, algebraic roots, exp, log,
trigonometric functions and their local inverses, on specified finite chart/piece domains.
The piecewise smooth step s is permitted. Unevaluated quadratures, infinite sums or iterations,
and naming the inverse of a newly defined function are excluded. Matrix adjugates and
determinants are finite polynomial operations. This convention is stronger than §8's exact
expression class; it does not change the campaign's original goal.

In particular (8.14) is an exact global prescription, **not** a global elementary coefficient
formula. Even the earlier cap coframe C9 still uses H and dH. Section 9.3 removes those
dependencies on two complete pole neighbourhoods. No general impossibility theorem for
elementary complex structures on S⁴×T² is claimed.

### 9.2 Eliminate both scalar profile inverses

These replacements define a new parametrization; they are not identities with the old
profiles. Their identities are L10–L11 in `sym/lemmas_1d.py`; inequalities below are (b)/(c),
with arithmetic checked in I2q–s. The old implementation remains the default.

**Radial block.** For w≠0, put R=|w| and O=ρ_O. Replace P_final by

  v=2Ow/(1+R²),  X=(R−1)/√(1+R²).                              (9.1)

Then |v|+X²=1, dX/dR=(R+1)/(1+R²)^(3/2)>0, and −1<X<1. Its inverse is

  a=X√(2−X²),  w=O⁻¹v/(1−a).                                  (9.2)

The denominator is positive: a²=1−(1−X²)²<1. For numerical stability when X>0,
use 1/(1−a)=(1+a)/|v|². L10 verifies both compositions and the denominator identity.
Thus ν and ν⁻¹ are removed, including at the equator.

**Base profile.** Write A=(4π)², B_e=e⁴, and replace k by

  k_A(r)=√((A+B_e r²)/(1+A r²)),  |T|²=(A y+B_e y²)/(1+A y), y=r².

The derivative in y is (A+2B_e y+AB_e y²)/(1+Ay)²>0, with range [0,∞). The endpoint
scales are k_A(0)=4π and k_A(∞)=e²/(4π). For u=|T|² the inverse is

  y=[√(A²(1−u)²+4B_e u)−A(1−u)]/(2B_e),  r=√y.               (9.3)

For u≤1 use the equal expression y=2u/[A(1−u)+√(A²(1−u)²+4B_e u)] to avoid
subtraction. Both k_A(|t|) near t=0 and 1/k_A(1/|s|) near s=0 are smooth functions
of the respective squared modulus. The h_B phases and their clutching identity C15a
are unchanged. Radial monotonicity, I1/C5 and that clutching identity therefore prove
that the revised raw h_B is still bijective, including its two coordinate circles.

**Recheck the shells.** The old I2 constants are not silently reused. Define

  β(ρ)=√[(1+e⁴ρ²/(64π²))/(1+4π²ρ²)],
  K_A(p)=(βp₁,βp₂,p₃),  ρ²=p₁²+p₂².

Before the final radial block, the new raw node maps are

  B₀,A(p)=B(2πK_A(p)),  B_∞,A(p)=O_out B(2π I K_A(p))/2.       (9.4)

B, I and O_out are §4.3's maps. L11d proves that the same β occurs at both cusps:
k_A(ρ/2)/(4π)=K_∞/k_A(2e⁻²/ρ). On these balls the source height is still p₃.
Together with I2b/d this derives (9.4) from the raw composition, with a regular node limit.

For R_b=1/125, set a=e⁴/(64π²), b=4π². The previous rational bounds on π,e⁴ give
0<a<b<40. Consequently

  40/41 < β ≤ 1,
  0 ≤ 1−β=(b−a)ρ²/[(1+bρ²)(1+β)] ≤ 40ρ²,
  |ρβ′|=(b−a)ρ²/[(1+bρ²)²β] ≤ 41ρ².

The first bound follows from β²≥1/(1+40R_b²)>(40/41)². The transverse eigenvalues
of DK_A are β and (1+2aρ²+abρ⁴)/[(1+bρ²)²β]=β+ρβ′, both in (0,1]; the third
is 1. Thus |K_A(p)|≤|p|, ‖DK_A‖=1 and ‖DK_A−1‖≤81|p|². Apply §4.3's bound
‖DB(z)−Q‖≤(9/2)|z| on |z|≤51/1000 to obtain

  ‖DB₀,A−M₀‖ ≤ 18π²|p|+162π|p|² ≤ 182|p|,
  ‖DB_∞,A−M_∞‖ ≤ 91|p|.

Integration along the segment from 0 to p gives value defects ≤91|p|² and
≤(91/2)|p|². Use the same blend f_j,A=M_jp+c_b(|p|)(B_j,A−M_jp). Since
|c_b′|≤500, the defect Lipschitz bounds are

  L₀=182/125+500·91/125²=546/125=4.368,  L_∞=273/125.

Hence σ_min(Df₀,A)>19/10, σ_min(Df_∞,A)>19/20, and the projected-defect
iteration (8.1) has ratio **<7/10** at both cusps. I2q–s verify the rational comparisons
using the same Machin enclosure. Applying the same argument to the smaller raw defect
proves equality of raw/modified ball images, as in §8.1. This proves global injectivity
of the shell replacement. Numerical margins (4000 points per cusp, I2 script): minimum
singular value 5.636589 / 2.815100 versus proved 1.9 / 0.95; maximum defect-derivative/scale
0.183538 / 0.185246 versus proved <0.7.

**Combine the two replacements.** On each shell put v₀=f_j,A(p), h₀=|v₀|, and

  v=2v₀/(1+h₀²),  X=±(1−h₀)/√(1+h₀²),                       (9.5)

with − at cusp 0 and + at ∞. Outside use the revised raw h_B, stereographic
projection and (9.1). They agree to all orders on each shell boundary: in the south
raw v₀=Ow, and in the north raw v₀=Ow/|w|². Both give (9.5). The raw bounds give
h₀<1/16<1, so this change is invertible there. For a target (v,X) first form

  v₀=v/[1+|X|√(2−X²)].                                       (9.6)

Use its shell contraction and domain test. Outside that image use (9.2), inverse
stereographic projection, (9.3), and the existing height inverse. The equality of raw
and modified ball images is preserved by (9.5). Thus the resulting base Λ_A is globally
bijective, with smooth inverse on the free part. In the inner balls h₀=scale·rad is
polynomial in the four slice coordinates, and (9.5) has the following smooth lifted model.

**Global assembly.** Keep §8.2's connection and replace Λ_Q by Λ_A in §8.4's paths.
Put ε=1/100, ε_A=2ε/(1+ε²)=200/10001 and x_±=±99/√10001. Replace L_j by

  Z=√[2/(1+h₀²)] Z₀,  X=±(1−h₀)/√(1+h₀²),
  Z₀=√scale·U_j(√a₀ x̄_s,y_s),  h₀=|Z₀|².                   (9.7)

In (8.6)–(8.7) replace ε by ε_A and use these cap maps and their inverses. Set
χ_X=s((X−x_−)/(x_+−x_−)); compute the new mismatch and logarithm using the new paths.
No old mismatch is assumed to equal the new one. Positive radial scaling preserves
the Hopf connection, so transport is exactly radial in these caps. The overlap,
inverse and smoothness proofs of §§8.3–8.6 apply with these substitutions, giving
global maps G_A,F_A. Both scalar profile inversions have been removed from their
actual prescription, rather than just from an isolated candidate block.

Implementation: `tlift.transport_base_chart/middle` and `tinverse.transport_base_inverse`
with **both** `algebraic_radial=True, algebraic_profile=True`; defaults retain §8.
`transport_base_local` and `transport_shell_inverse` take the latter flag.
Height inversion and the convergent shell iteration are still required.

### 9.3 An algebraic complex-structure tensor on both pole neighbourhoods

Change the T² trivialization by a global translation cancelling H near the nodes.
In a cusp let n be the nearest integer to H, reduce by ψ̂⁻ⁿ, and denote the reduced
point by (x₀,y₀), H₀=H−n and r₀=a₀|x₀|²+|y₀|². Define, with zero extension,

  η(r₀)=1−s((r₀−1/250)/(1/250)),
  b_j=η(r₀)H₀ c̃_j(t_c),  b=b₀+b_∞,  φ_E=φ+b.                (9.8)

These are real two-vectors before reduction modulo ℤ². The supports are disjoint:
|t_c|≤e²R_b/2<4/125 there, which cannot hold simultaneously for t and 1/t.
I1 gives |H₀|≤R_b+(27/10)R_b²<1/100 on the support, so the reduced representative
is unique on a neighbourhood of it; at an integer-choice seam η is zero nearby.
A deck change increases H,n together and preserves H₀ and the reduced point (C17a).
Thus b is globally smooth on W, including its nodes and support boundaries.
Translation φ↦φ+b(w) has inverse φ=φ_E−b(w) and determinant one. These are identities
and explicit support bounds (C17a/b/e), with smoothness from the displayed flat cutoff.

Define **G_E(q,φ_E)=G_A(q,φ_E−b(G_{A,W}(q)))**; its inverse has torus component
φ_E. This is a global change of parametrization of Y. Its tensor J_E=G_E*J_Y has
the same global pullback/overlap proof as §8, and still uses quadratures off the caps.

On either cap h₀≤ε, r₀=h₀/scale<1/250, so η=1,n=0 and C=φ_E **exactly**.
The complete elementary cap map, including the pole, is as follows. Put r=|Z|²,

  s=√(1−r²),  g=(1+s)^(-1/2),  Z₀=gZ,
  (q₁,q₂)=U_j* Z₀/√scale,  x_s=q̄₁/√a₀,  y_s=q₂,  t=x_s y_s,
  W=φ_E,2−iσφ_E,3,
  (x,y,w)=(x_s e^(−tW), y_s e^(tW), σ(W̄−Wt²)),              (9.9)

where σ=1 at cusp 0 and −1 at ∞, and U_j,scale,a₀ are as in §3.4.
The identity h₀=r/(1+s) follows from (9.7); r≤ε_A<1/50. No height, series,
inverse function, limit or integral occurs here. Integer changes of φ_E act by
the given holomorphic T²-action, so there is no torus-angle seam (C9d with H=0).

The tensor itself is algebraic in the ambient coordinates and the constants π,e:
torus invariance permits W=0 when computing its coordinate coefficients, retaining dW.
Write Z_R for the four real components of Z. First take

  dZ₀=D_g dZ,  D_g=g I₄+2g′ Z_R Z_Rᵀ,  g′=r/[2s(1+s)^(3/2)],
  D_g⁻¹=g⁻¹[I₄−r Z_R Z_Rᵀ/(1+s)].                            (9.10)

Apply (9.9)'s constant unitary/conjugation map to obtain dx_s,dy_s, and set
dW=dφ_E,2−iσdφ_E,3. Form the complex coframe

  η₁=dx_s−x_s t dW,  η₂=dy_s+y_s t dW,  η₃=σ(dW̄−t²dW).     (9.11)

Let M_E be its real coefficient matrix in (Re Z₁,Im Z₁,Re Z₂,Im Z₂,φ_E,2,φ_E,3).
With E,R as in (8.15), the complete formula is

  **J_E^amb=E adj(M_E) J₀ M_E R / det(M_E),**                  (9.12)
  det(M_E)=−σ(1−|t|⁴)g⁴/(scale²a₀s) ≠ 0.                    (9.13)

L10e/f prove the radial eigenvalues g,g,g,g/s and (9.10). C9l gives the constant
slice determinant and C17d the triangular coframe determinant σ(1−|t|⁴).
Also |t|²≤ε²/(4scale²a₀)<1/1000 by e⁴<55 and π>157/50. Thus
1−|t|⁴>999999/1000000, s>√(1−1/2500)>0 and g⁴≥1/4. These explicit (b)/(c)
bounds include Z=0. The coframe is a complex-linear rescaling of the differential
of the actual map (9.9) (C17c), so C10 proves N_J=0. The global inverse and overlap
arguments above make this the restriction of J_E, not an unrelated local field.

`jround.elementary_cap_map`, `elementary_cap_coframe` and `J_elementary_cap`
implement (9.9)–(9.12); `tinverse.algebraic_cap_slice` implements the slice inverse.
`tlift.elementary_torus_shift_chart` requires a reduced cusp chart;
`elementary_torus_shift_middle` evaluates the global correction on the regular part.

### 9.4 Remaining strategies and restricted impossibility results

**Special-function compression works for the height value only.** On t≠0 set
P<0, ρ=log|ζ|, q=e^P and z=i(ρ−P/2). Then

  H(ρ,P)=θ₃′(z,q)/(2i θ₃(z,q)),                               (9.14)

in the convention of [NIST DLMF 20.2.3](https://dlmf.nist.gov/20.2.E3).
C16a/b check equality term by term with the Gaussian sum; normal convergence
justifies differentiation. On this imaginary argument the denominator is a sum
of positive real terms. This is a standard special-function expression, not an
elementary one, and does not evaluate H⁻¹. Naming that inverse would not solve
the task. At t=0 use the exact axis formulas, not (9.14) with P=−∞.

**No algebraic equivariant height.** Here is a restricted impossibility proof (c),
with explicit constants. Fix t≠0, set Y=|ζ|² and q=e^(2P)∈(0,1). Equivariance gives

  h(qY)=h(Y)+1,  h(q^nY₀)=h(Y₀)+n.                            (9.15)

Suppose a nonzero polynomial relation Σ_{j=0}^m a_j(h(Y))Y^j=0 holds on Y>0.
Divide out the least power of Y, so a₀ is a nonzero polynomial. Fix Y₀>0, write
a_j(h(Y₀)+n)=Σ_l b_jl n^l, and let D be the maximum degree. Substitution yields

  |a₀(h(Y₀)+n)|≤C n^D q^n,  C=Σ_{j=1}^m Σ_l |b_jl|Y₀^j,     (9.16)

for integers n≥1, by the triangle inequality and q^(jn)≤q^n. Put c=−log q>0.
For D≥1, x^D e^(−cx/2) has maximum K=(2D/(ec))^D at x=2D/c, by its displayed
derivative in L12a. For D=0 take K=1. The right side of (9.16) is ≤CK e^(−cn/2).
If d=deg a₀, then b₀d≠0, and for
n≥max(1,2Σ_{l<d}|b₀l|/|b₀d|) the left side is ≥|b₀d|n^d/2≥|b₀d|/2.
For C>0 choose also n>(2/c)max(0,log(2CK/|b₀d|)); this is a contradiction.
For C=0 the contradiction is immediate. **No such algebraic height exists.**

This excludes a fixed finite Gaussian moment sum, which is rational in Y when
t is fixed. It also excludes an algebraic formula for Y as a function of the
original height, since that would give the same polynomial relation. It does not
by itself exclude an elementary inverse for ρ=(log Y)/2. The original H is analytic on Y>0
by normal convergence on complex annuli; a polynomial relation on any open interval
would therefore extend to all Y>0 and give the same contradiction. This does not
exclude elementary heights: log Y/log q solves (9.15). That particular candidate
has the known central-fibre collapse of notes §13.3 and cannot replace H here.

**An endpoint potential cannot reproduce all transports of this connection.**
On the middle region α=2πdx₁+θdx₄, θ=arg t, and its curvature is dθ∧dx₄.
The base-form integral around a rectangle is (θ₁−θ₀)(x₄,1−x₄,0), by the four
line integrals (L12b). Side lengths 1/10 give 1/100, not 0 modulo 2π. Thus a single
endpoint potential cannot reproduce this transport for all paths. This does not
exclude an elementary formula along the particular meridians of §8.4 or a different
connection/bundle isomorphism.

**Principal Arg is not automatically a smooth global logarithm.** For example,
g(n)=exp(4in₃) has smooth real lift 4n₃ but crosses the principal cut at n₃=π/4.
Simple connectedness gives a real lift, not a bound on its range. For the actual
g_A neither a certified principal-branch bound nor a closed replacement for (8.8)
has been established.

**Analytic-only smoothing cannot keep an exact cap and change it elsewhere.**
A real-analytic function on a connected interval that equals a fixed affine function
on a subinterval equals it throughout, by power series and continuation. This
restricted observation does not exclude the allowed piecewise elementary smoothing.
For example, e^(−1/u), extended by 0 for u≤0, is flat at 0: e^v≥v^N/N! for
arbitrary N gives superpolynomial decay; every derivative is a finite sum of powers
of 1/u times that exponential. L12c checks low orders; the arbitrary-order proof
is this positive-series bound, not an extrapolation from those checks.

### 9.5 Completion status and verification

| dependency | outcome |
|---|---|
| ν⁻¹ in P_final | eliminated globally by (9.1)–(9.2) |
| inverse of r k(r) | eliminated globally by (9.3); revised shells certified |
| H,dH in the pole tensor | eliminated by (9.8); algebraic tensor (9.12) |
| general height value | standard θ₃ expression (9.14), still not elementary |
| general height inverse | unresolved; an algebraic inverse for Y excluded by (9.15)–(9.16) |
| shell inverse | contraction with ratio <0.7; no finite formula obtained |
| meridian transport and phase logarithm | quadratures remain; unrestricted endpoint-potential/Arg shortcuts fail |

These reductions and restricted negative results do **not** exhaust possible constructions.
Nor do they exclude cancellations of difficult intermediate functions in the tensor itself.
A finite height/gauge redesign, closed evaluation of the specified integrals, or another
global parametrization remains possible. **The full elementary global expression is open;
its impossibility has not been proved.** The exact quadrature construction is retained.

New numerical checks (d), separate from the proofs: C16 compares θ₃ and the independent
height sum at 25 points and 60-digit precision, discrepancy ≤9.3·10⁻⁶¹. Each algebraic
scalar block has 4097 inverse tests over 24 orders of magnitude, relative errors
≤2.4·10⁻¹⁵ / 5.8·10⁻¹⁶. The revised base has 4560 random/shell/pole/axis inverse tests,
residual ≤1.3·10⁻¹⁴, and all 32 pole second-jet comparisons. Each new shell has 4000
independent composition/derivative comparisons, errors ≤8.6·10⁻¹⁵ / 2.4·10⁻¹³.
The new tensor has 256 tests per cap against the Jacobian of the elementary map at
nonzero torus coordinates, plus integer torus changes, J²=−P, tangency, the determinant,
exact poles and axes. Coefficient discrepancy ≤8.9·10⁻¹⁶; independent finite-difference
Nijenhuis tests on 44 points per cap give residual <1.3·10⁻¹⁰ after division by
max(1,cancellation scale); eight-point nonintegrable controls have median residual
0.26 / 0.22 under the same normalization. All 14 central-axis derivatives through
order three per cusp agree with independent Taylor polynomials to <10⁻¹⁰.
There are also 1024 torus-correction chart/middle/deck/branch comparisons per cusp,
exact cancellation on the inner caps, 256 exact pole/axis values per cusp, and
1024 revised-base chart/middle/deck/branch comparisons per cusp. These cap checks
are not a complete global numerical evaluation of J_E.

## 10. Literature-led continuation: two height replacements and an algebraic target

This section records further reductions, **not completion of the global closed-form
tensor**. There are two different height designs: §10.2 is finite elementary in the
forward direction; §10.5 has a closed inverse in standard elliptic functions. They
must not be conflated. The latter is an extension of the elementary expression class.
The original Gaussian construction and all its defaults remain available.

### 10.1 What the checked literature supplies

The following are source-specific conclusions, not additional assumptions about Y.

| primary source inspected | relevant result | consequence for this campaign |
|---|---|---|
| Di Scala–Kasuya–Zuddas, [*Non-Kähler complex structures on R⁴*, §2–3](https://msp.org/gt/2017/21-4/gt-v21-n4-p13-p.pdf) | a smooth Hopf suspension, its composition with Hopf, and an analytic nodal-torus model; the final smooth identification uses gluing/isotopy | construct a specified algebraic target (§10.4), but do not claim the paper supplies a coordinate formula for our diffeomorphism |
| The same authors, [*On embeddings of almost complex manifolds in almost complex Euclidean spaces*, §4](https://arts.units.it/retrieve/e2913fdd-2206-f688-e053-3705fe0a67e0/PHE25.pdf) | separate the two critical values by an orientation-preserving diffeomorphism of S³ | replace that unspecified map by the explicit boost below |
| Zung, [*A note on focus-focus singularities, II*, §4](https://www.math.univ-toulouse.fr/~tienzung/Maths/Focus2_lmp2002.pdf) | circle reduction and topological gluing; the remark corrects an earlier smooth-classification claim for multiple pinches | topology alone does not certify the required smooth lift |
| Bolsinov–Izosimov, [*Smooth invariants of focus-focus singularities and obstructions to product decomposition*, Theorem 3.3 and §3.2](https://arxiv.org/pdf/1706.07456) | liftable plane germs are divisible smoothly by z or by its conjugate; gluing carries smooth information | supports checking exact factorization, as in L8; this theorem concerns their focus-focus model and is not automatically a theorem about our Hopf lift |
| Massamba–Thompson, [*The Universal Connection and Metrics on Moduli Spaces*, §3](https://arxiv.org/pdf/math/0311198) | finite-dimensional bundle embeddings and universal connection matrices | suggests a finite projection construction; the required concrete embedding, nonvanishing estimates and pole compatibility still have to be supplied |
| [DLMF 22.2](https://dlmf.nist.gov/22.2), [22.11](https://dlmf.nist.gov/22.11), [22.16](https://dlmf.nist.gov/22.16) | theta/modulus relations, the Fourier series of dn, and amplitude/incomplete-integral inversion | closes the redesigned scalar height and its collar inverse in §10.5 |

In particular, neither matching the number/signs of the singularities nor invoking
an isotopy supplies its finite formula. The multi-pinched smooth invariants in the
fourth source do not themselves obstruct this example with separate nodal fibres.
No new dependency on the manuscript's 2026 inputs is introduced.

### 10.2 A finite elementary cusp height

Put c=1/16 and let s be the existing flat smooth step. Define K on u≥0 by

  K(u)=u,                                                     0≤u≤1/16,
  K(u)=u+s(log(16u)/log 8)(u/(1+u)−u),                         1/16<u<1/2,
  K(u)=u/(1+u),                                               1/2≤u≤2,
  K(u)=1−K(1/u),                                             u>2.       (10.1)

The last line uses only the preceding lines and is not an infinite recursion.
It gives K(u)+K(1/u)=1 for u>0. Flatness of s gives agreement of every derivative at
1/16,1/2,2,16. On the lower transition, writing λ=s(log(16u)/log 8),

  K′=(1−λ)+λ/(1+u)²−λ′u²/(1+u),  0≤λ′≤2/(u log 8).

Thus K′≥4/9−2/(3 log 8)>1/9. Here log 8>2 follows from
log 2=2 atanh(1/3)>2/3, using the positive power series. On the reciprocal
transition K′(u)=K′(1/u)/u²>0; all other pieces have displayed positive derivatives.
These are identities (a) and manifest estimates (b), checked in L13.

In a cusp chart put X=|x|², Y=|y|², v=XY, a₀=ν², and

  e=exp(−4v/(1−v)),  A=a₀eX,  q=a₀ev=AY.

For |t|≤1/2, 0≤q<1/8. Indeed e≤1/(1+4v/(1−v)) and a₀≤1, so
q≤v(1−v)/(1+3v); the inequality with 1/8 has numerator
8v²−5v+1=8(v−5/16)²+7/32>0 (L13d).

For the proof only, set S_q(u)=Σ_{j≥0}K(q^j u). If q^N u≤1/16, then exactly

  S_q(u)=Σ_{j=0}^{N−1}K(q^j u)+q^N u/(1−q).                  (10.2)

Consequently this is a finite evaluation, with an exact geometric tail, not a
truncated approximation. Define H_f=S_q(A)−S_q(Y). On the regular part let
P=(log q)/2, ρ=log|y|, n=floor(ρ/P), Y₀=exp(2ρ−2nP), A₀=q/Y₀. Then
q≤Y₀,A₀≤1, and the entire formula is just

  H_f=n+K(A₀)+K(qA₀)−K(Y₀)−K(qY₀)+q²(A₀−Y₀)/(1−q).         (10.3)

This has four seed evaluations. The floor chooses a representative, as in the
existing atlas. It does not create a discontinuity: S_q(u)−S_q(qu)=K(u), and
K(Y)+K(1/Y)=1 give H_f(qY)=H_f(Y)+1 exactly. Formula (10.2), with a locally
fixed N, proves smoothness across all artificial representative seams and q=0.
On a compact chart, choose N with (1/8)^N max(A,Y)<1/16; then all terms and the
rational tail are smooth on a neighbourhood. L13 verifies the recurrence and
geometric-tail identities, independently of the numerical evaluator.

At fixed nonzero t the derivative is

  ∂ρH_f=−2Σ_{j≥0}[q^j A K′(q^j A)+q^j Y K′(q^j Y)]<0.       (10.4)

The positive tail is geometric when its arguments reach the linear piece; thus
termwise differentiation is justified directly by (10.2). Equivariance makes the
range all of R. Reflection A↔Y gives H_f(P−ρ)=−H_f(ρ), and H_f(0)=−1/2.
On the central fibre H_f(x,0)=K(a₀X), H_f(0,y)=−K(Y); no axis is collapsed.

### 10.3 Certificate and scope of the finite-height substitution

On R=a₀X+Y≤1/20 both A,Y≤1/16, so

  H_f=(A−Y)/(1−q),  −∂ρH_f/2=(A+Y)/(1−q).                   (10.5)

Here is a uniform explicit-constant estimate (c). Write u=a₀X, v=Y just for
this paragraph, m₀=u−v, r=u+v, and z=uv/a₀. The earlier e⁴<55 gives
z≤55r²/4≤11/320, q=euv≤r²/4≤1/1600, and 1−e≤4z/(1−z). Hence

  H_f−m₀=[−u(1−e)+q(u−v)]/(1−q),
  |H_f−m₀|≤[16·55/(27(1−11/320))+1/4]/(1−1/1600) · r³
           <35r³.                                           (10.6)

We used u²v≤4r³/27 and uv≤r²/4, whose square-factor identities are already
checked. Also (A+Y)/((1−q)r)≥e≥1−4(11/320)/(1−11/320)=265/309.
For exactly the existing node cutoff c_n=1−s((r−.02)/.03), the corrected height
H_f+c_n(m₀−H_f) therefore satisfies

  −∂ρH_f,new/(2r)≥265/309−(2·35/.03)(.05)³>14/25.            (10.7)

Furthermore |H_f|≤r+35r³<.055 throughout the support. Thus the reduced node
integer is unique, including across nearest-integer seams. N6 checks the constants.
Measured margins on 200000 samples per cusp are .998628 / .903375; measured
cubic remainder constants are .528244 / 32.280107, against 14/25 and 35 respectively.

The corrected height is still exactly m₀ on r≤.02. Therefore the local base
formulas, the I2 certificates and both elementary pole tensors of §9 are unchanged.
Replace H consistently in the cusp height, the unwrapped collar, the torus
trivialization and the connection. Equivariance, positivity (10.7), the same
middle matching, and the unchanged local shells supply the hypotheses of §8's
assembly proof. This defines a further **quadrature construction**, G_f. Its
complete tensor is still not a finite closed form.

Implementation: `T.FINITE_CUSP0/INF`, `finite_height_seed`, and `cusp_angle`;
the independent base/inverse, connection and torus correction accept
`finite_height=True`. Reduced chart functions instead take the finite cusp object.
These are separate from the original `Psi_*`, `Finv`, `tfull` and `J_round`.
In particular, there is no complete numerical API for G_f/J_f. The implemented
inverse still bisects (10.3) and the scalar collar. The finite seed's transition
does not have a proved elementary inverse; the scalar problem has not disappeared.

### 10.4 A specified algebraic round-S⁴ target

The literature's suspension can be written, with |z₁|²+|z₂|²+X²=1, as

  Σh=(2 Re(z₁ conj z₂),2 Im(z₁ conj z₂),|z₁|²−|z₂|²,X√(2−X²)) ∈ S³.

For p=(a,b,c,d)∈S³ define

  k(p)=((√2 a+1),b,c,d)/(√2+a),
  k⁻¹(p)=((√2 a−1),b,c,d)/(√2−a).                            (10.8)

Both denominators are ≥√2−1. Direct substitution gives the sphere and inverse
identities; Dk has conformal factor 1/(√2+a) on tangent vectors (C18a–d).
Define f_alg=h∘k∘Σh, using h(a,b,c,d)=(2(ac+bd),2(bc−ad),a²+b²−c²−d²).
This is a global algebraic map S⁴→S² with f_alg(0,0,0,0,ε)=(0,−ε,0).

It has no other critical points: writing R=√(1−X²), the three nonzero singular
values of DΣh are 2R,2R,2R/√(2−X²). The Hopf differential has its two nonzero
singular values equal to 2. Surjectivity of these differentials proves rank two
off the poles, quantitatively

  σ₂(Df_alg)≥4R/((√2+a)√(2−X²))≥4R/(2+√2)>R.                (10.9)

At each pole all first derivatives vanish. In the two target tangent coordinates
(first and third components), the quadratic part is
(p₃+εp₂,√2 p₁), where p=Hopf(z₁,z₂). The determinant of the Hessian of the
real pencil λQ₁+μQ₂ is 64(λ²+μ²)² (C18f/g). These are identities (a) and
manifest inequalities (b). This paragraph does not infer an exact local normal
form or a fibrewise equivalence with W from the Hessian alone.

`suspended_hopf`, `sphere3_boost`, and `algebraic_torus_fibration` implement this
target. The route still needs an explicit smooth matching of the source fibres,
its global gluing and its circle lift. Thus f_alg is **not a complex structure on
S⁴**, nor a completed construction of J on S⁴×T².

### 10.5 A closed height inverse in standard special functions

This alternative really evaluates an inverse, rather than renaming the old one.
It replaces the Gaussian height, not merely its notation. At t≠0 set

  q=e^(2P), B=θ₃(0,q)², κ=θ₂(0,q)²/B, k=√(1−κ²), K(k)=−PB,
  h_e(ρ,P)=am(−2Bρ,k)/π−1/2.                                (10.10)

Here am is the continuous real Jacobi amplitude, not principal arcsin(sn) on
every period. These conventions and K(k)=−PB follow from the complementary
nome identities in DLMF 22.2.1–2. DLMF 22.16.2 and 22.16.10–11 give

  h_e(ρ+P,P)=h_e(ρ,P)+1,  ∂ρh_e=−2B dn(−2Bρ,k)/π<0,
  ρ=−F(π(a+1/2),k)/(2B),  when h_e(ρ,P)=a.                  (10.11)

F is the standard real incomplete elliptic integral of the first kind with its
continuous amplitude convention. Numerically reduce a=n+δ, −1/2≤δ<1/2,
and use nP−F(π(δ+1/2),k)/(2B), so every branch is explicit.

For the nodal normalization use the degree-one real lift

  R_c(h)=n+atan2(c sin(πδ),cos(πδ))/π,
  n=floor(h+1/2), δ=h−n, c>0.                               (10.12)

Its inverse is R_(1/c), and R_c′=c/(cos²(πh)+c²sin²(πh))>0.
Set H_e=R_(π/2)(h_e). In (10.11) replace a by R_(2/π)(a) to invert H_e.

**Central extension and verification of the identity.** With A=ν²e^(−2f)X,
Y=|y|² and q=AY as above, an equivalent expression is

  h_e=(2/π)Σ_{j≥0}[atan(q^j A)−atan(q^j Y)].                 (10.13)

Here the series proves smoothness; the closed expression and inverse remain
(10.10)–(10.12). To verify (10.13), differentiate with respect to ρ. Its
derivative is −(2/π)Σ_{n∈Z}sech(2ρ+n log q). Put L=−P>0.
The Fourier transform of sech(2x) is π/(2 cosh(πω/4)): substitute u=e^(4x)
and use the [beta integral](https://dlmf.nist.gov/5.12.E3) and
[gamma reflection formula](https://dlmf.nist.gov/5.5.E3). Thus the periodized
derivative has constant coefficient −1/L and cosine coefficients
−(2/L)sech(π²n/(2L)). The derivative of (10.10) has the same coefficients by
DLMF 22.11.3, K=LB and K′=πB/2; C19a/b check the coefficient simplifications.
Both heights equal −1/2 at ρ=0. This proves equality, rather than relying on
sampled agreement. Exponential decay justifies the Fourier summation.

Near the node, (10.13) converges normally as a holomorphic series in the independent
variables (A,Y,q): on |A|,|Y|,|q|≤1/8 its tail starting at j=N is bounded by
(64/63)(|A|+|Y|)|q|^N/(1−|q|), from the arctangent power series. This proves a
smooth real extension after substitution of the smooth chart variables. At either
non-nodal axis finitely many initial terms have positive finite arguments; the
remaining tail has the same local argument. In particular

  H_e(x,0)=R_(π/2)((2/π)atan(a₀X)),
  H_e(0,y)=−R_(π/2)((2/π)atan(Y)),
  rad=tan((π/2)R_(2/π)(|a|))                                (10.14)

on the corresponding central axis. Thus its central-axis inverse is elementary.
No logarithm of zero is evaluated at the node or on an axis.

**A certificate near the node.** Let L_a be the sum of arctangents in (10.13),
without 2/π, and let H_lin=(A−Y)/(1−q). On r≤.05,

  |L_a−H_lin|≤r³/[3(1−q³)],  |L_a|≤r/(1−q),
  H_e=T(L_a),  T(s)=atan((π/2)tan(2s))/π,
  T′(s)=1/[cos²(2s)+(π²/4)sin²(2s)],  |T(s)−s|≤2|s|³.

The last inequality follows from 0≤1−T′≤6s² (π²<10) and integration of
the polynomial bound. Combining with (10.6) gives |H_e−m₀|<38r³. Its fibre
slope satisfies

  −∂ρH_e/(2r)≥[265/309]/(1+.05²) · [1−6(.05/(1−1/1600))²].

If one applies the old node correction, the lower bound after subtracting
(2·38/.03)(.05)³ is >13/25. C19c/d and N7 check these explicit constants.
The 200000 samples per cusp in N7 give cubic remainder constants
2.289646 / 32.441924 and corrected slopes ≥.982202 / .899635, respectively,
against the proved constants 38 and 13/25. These samples evaluate eight terms
of the independent arctangent sum; on this ball the omitted tail is <2·10⁻²⁷.
That correction would reintroduce an implicit inverse in its transition; it
has **not** been silently appended to the closed height block.

**The whole cusp collar can also be inverted.** A convex blend of two invertible
heights need not retain a closed inverse. Instead put b(r)=1−s((r−.3)/.2),
k_b=b(r)k, K_b=K(k_b), c_b=1+(π/2−1)b(r), and define

  a_e=R_(c_b)(am(2K_bρ/P,k_b)/π−1/2),
  ρ=P F(π(R_(1/c_b)(a_e)+1/2),k_b)/(2K_b).                  (10.15)

Use the same period reduction as in (10.11). At r≤.3 this is H_e; at r≥.5,
k_b=0, K_b=π/2 and c_b=1, so it is exactly ρ/P−1/2 with inverse P(a_e+1/2).
All jets match at both boundaries because the cutoff is flat. At fixed t,
strict monotonicity follows from dn>0 and R_c′>0, with no blend-derivative term.

Scalar reference functions in `cusp.py` implement the forward value, inverse,
normalization, exact axes and complete collar. They use mpmath's standard sn,
cn, theta and elliptic F; mpmath's argument is the parameter k², not DLMF's
modulus k. Extra working precision resolves 1−k² when q is tiny. C19 compares
120 forward/inverse cases with independent arctangent sums at 80-digit precision,
including P=−100, translated cells, exact seams and near-seams. It also checks
20 axis values/inverses, the node, 72 complete collar inverses, and 36 seam jets
through order three against the independent sum. Maximum discrepancies are
1.5·10⁻⁸⁰, 1.1·10⁻⁸¹ and 1.1·10⁻⁸¹ for the 120 value/inverse, collar and jet
groups, respectively. The evaluator retains its guard digits in the returned
mpmath value. Rounding that value before inversion can lose accuracy near a
nearly flat fibre; the inverse now converts its input inside the guarded context.

**Scope:** (10.15) closes the scalar height/collar block in standard special
functions. It does not invert the existing node-correction blend or the I2 shell,
and it has not yet been assembled into a replacement global base/lift/tensor.
The finite elementary design of §10.2 and this elliptic design are alternatives.

### 10.6 Finite bundle operations and remaining strategies

**A finite real logarithm with an explicit derivative bound.** Suppose a specified
smooth g:S²→S¹ satisfies ||Dg||≤L in the round metric. This is an input requiring
a proof, not a sampled bound. Let N=floor(8L)+1. In either stereographic chart
S_N,S_S of (8.8), on |a|<2 define the finite expression

  P_±(a)=Σ_{j=1}^N Arg(g(S_±(ja/N)) conj(g(S_±((j−1)a/N)))).  (10.16)

Every argument lies in the right half-plane with real part >7/8. Indeed
DS_±^T DS_±=4I/(1+|a|²)² (C20a), so each path segment has length at most
2|a|/N<4/N. Its images differ in C by ≤4L/N<1/2, by integrating the stated
derivative bound along that segment. For unit complex numbers z,w,
Re(z conj w)=1−|z−w|²/2>7/8 (C20b). Thus each Arg in (10.16) means the single
smooth elementary expression atan(Im/Re); no winding or branch choice is hidden.

Use the same constants as (8.8): Θ_N=Arg(g(north))+P_N(a) and
Θ_S=Arg(g(north))+P_N(e)−P_S(e)+P_S(a), e=(1,0). Telescoping the finite phase
ratios (C20c) proves e^(iΘ_±)=g. Both are smooth real functions. Their difference
on the connected band −3/5<n₃<3/5 is a continuous multiple of 2π and vanishes
at (1,0,0), so it vanishes throughout the overlap. This gives a global smooth
real logarithm. For evaluations choose north when n₃≥0, south when n₃≤0;
then |a|≤1. The endpoint derivative bounds and identities above are types (c)/(a).
This removes the phase-logarithm quadrature **provided a finite formula for g
and a certified L have been supplied**. It does not assert those inputs for the
campaign's unfinished finite lift. No bound on the total phase range is needed.

**A finite circle-bundle isomorphism.** Let P₀,…,P_M be smooth rank-one Hermitian
projections with ||P_j−P_(j−1)||<1. If P_(j−1)z=z and ||z||=1, put

  z_j=P_j z_(j−1)/||P_j z_(j−1)||.                            (10.17)

The denominator is ≥1−||P_j−P_(j−1)||>0 by the triangle inequality. This is
smooth and S¹-equivariant. If u,v are unit frames of two consecutive lines and
a=v* u≠0, the forward map multiplies the frame by a/|a| and the reverse by
conj(a)/|a|; their product is 1 (C20d). Hence reversing the finite chain is its
exact inverse. For a projector path with ||dP/ds||≤C on an interval of length D,
M=floor(2CD)+1 equal steps ensure consecutive distances <1/2. This bound follows
by integrating the derivative estimate. The finite operation is a bundle map in
its own right; it does not approximate or claim equality with connection transport.

`tlift.finite_phase_log` and `finite_projection_lift` implement (10.16)–(10.17).
The former requires the certified derivative bound as an explicit argument.
Tests use g(n)=exp(i(12n₁+9n₂n₃+5n₃)), whose derivative is bounded by
12+9+5=26 since |grad(n₂n₃)|≤1 on S². Its phase takes both −12 and 12, so
the test covers more than three turns. There are 259 value/exponential tests,
including both poles, and 65 real-logarithm and tangent-derivative overlap tests.
The projection tests use 128 independent 64-step chains and check nonvanishing,
image, norm, equivariance and reversal. These are numerical observations (d),
separate from the proofs above: phase-value error ≤7.1·10⁻¹⁵, real-overlap error
≤3.6·10⁻¹⁵, tangent-overlap error ≤3.4·10⁻¹⁴, and projection error ≤1.2·10⁻¹⁵.

The remaining applications are precise:

1. **Build the inverse base first.** With (10.15), the unmodified inverse chain is
   closed in standard functions. Blend that inverse to the linear pole model in
   target coordinates, and certify its derivative there. This could remove the
   shell solve from G and hence from J=G*J_Y, even if F's coordinate expression
   remains implicit. The new shell certificate and its smooth overlaps still need
   to be supplied; the old forward-shell certificate does not prove them.
2. **Apply the finite bundle maps.** Supply an explicit finite bundle embedding,
   the constants C,D in (10.17) along every meridian, pole agreement, and the
   resulting mismatch's constant L in (10.16). An existence theorem for embeddings
   supplies none of those concrete inputs. The two finite operations are complete;
   their application to W is not. No completed projection lift is claimed.
3. **Match directly to (10.8).** This avoids treating the target fibration itself
   as unspecified. It still requires a smooth fibrewise matching and global
   clutching calculation, including the source torus translation; equal singular
   data do not constitute that calculation.

The scalar advances do not close the circle lift or the modified-shell inverse.
The full finite global tensor and the campaign's requested finite inverse pair
are therefore still unresolved. No general impossibility result follows from
this investigation, and no new global numerical J evaluator is claimed.


## 11. Finite elementary global tensor and two-way inverse

The earlier implementation of this section gives a finite elementary expression
for the global tensor by an inverse-first argument. Its global statement is for
human review; the checks and practical limits are in §11.8. Sections 11.2,
11.4 and 11.6 give a two-way replacement, including F_cl. Its exact
identities and differential bounds have been checked as specified below;
the two-way global formula and its domains are certified there.

### 11.1 A concrete ten-component circle-bundle embedding

Write w(h)=1−s((|h|−1/4)/(3/20)). It is 1 on |h|≤1/4 and 0 on
|h|≥2/5. For either cusp, let H be its chosen smooth deck-equivariant height.
For c=0,1/2 choose k_c=floor(c+1/2−H) and (x_c,y_c)=ψ̂^{k_c}(x,y). Define

  V_c=w(H+k_c−c)(√a₀ conj(x_c),y_c)/√(a₀|x_c|²+|y_c|²).

Use zero for the entire pair wherever its weight vanishes. A deck translation
changes k_c by the opposite integer, so these pairs descend exactly. A
representative seam is surrounded by an open set on which that pair vanishes.
Every height class is within 1/4 of 0 or 1/2; consequently one of the two weights
is exactly 1. On |H|≤.6, the first active pair uses the original chart. The second
uses that chart for H>0, and ψ̂ once for H<0. The latter's active support has
H<−.1, so y≠0, including on the central fibre. Thus these are smooth through
either central axis, away from the node. This is a finite chart argument, not an
infinite deck sum. C21 checks the reindexing and phase identities.

On the middle use θ₀=Arg(t), θ₁=Arg(−t)+π and

  d_j=1−s((|θ_j−jπ|−π/2)/(π/4)),
  V_m,j=d_j exp(i[2πx₁−(θ_j−Im log t)x₄]),  j=0,1.

Each component vanishes near its own argument cut. The change log t→log t+2πi,
x₁→x₁−x₄ cancels in the exponent. Changing x₄ by an integer changes it by
an integer multiple of 2π. At least one d_j equals 1.

Let β₀,β_∞ be (8.4)'s cutoffs, and β_m=1−β₀−β_∞. Concatenate, rather than
add, the separate components:

  V=(β₀V₀,0, β₀V₀,1/2, β_∞V_∞,0, β_∞V_∞,1/2, β_mV_m,0, β_mV_m,1),
  e=V/√(V*V) ∈ S¹⁹⊂C¹⁰,             P=e e*.                 (11.1)

The four cusp entries in this display are pairs. At most two radial weights are
nonzero. Since each constituent squared norm is between 1 and 2,

  1/2≤β²+(1−β)²≤V*V≤2,   β²+(1−β)²−1/2=2(β−1/2)².        (11.2)

These are manifest inequalities (b), with the square identity in C21. Circle
rotation multiplies every component by the same phase. Therefore P is a smooth
rank-one Hermitian projector on the free orbit space, and e identifies each
source circle with the unit circle in im P. Relative to any reference e_ref,
the inverse on that circle is multiplication of the source reference by e_ref*e.
On rad≤.004, the height estimates give |H|<.1, and |t_c|<.3. Thus e is exactly
(√a₀ conj x,y)/√rad in the first pair of its own cusp block. All other components
are zero. This gives the required exact homogeneous pole models.

`tlift.bundle_embedding_middle/chart` implement this for both previous heights.
`cusp.elliptic_bundle_embedding` is the scalar closed-special-function version.
Its previously implemented elementary inverse-first option uses the known inverse
height coordinate; the two-way replacement below supplies its inverse. In middle
data that coordinate refers to the principal-log cusp representative. The scalar
implementation first replaces log t by its principal representative and replaces
x₁ by x₁+kx₄, where log t_old=log t_principal+2πik. This preserves (z,w) and
keeps the supplied height attached to the correct cusp representative.

### 11.2 A two-way elementary height and complete collar

This replacement of the inverse-first height is implemented in the scalar
evaluator and global forward F_cl evaluator. Its
elementary inequalities below are type (b)/(c); C29a–e,g check the specified
algebraic identities, while sampled height compositions and seam jets are
recorded in §11.8. Put
δ=1/10, c=1/5, η=1/100. The allowed parameter set is

  D = {0≤q≤1/8, λ=0} ∪ {1/2500≤q≤1/8, 0≤λ≤1}.

On its first piece define κ=2, γ=0 and d=1, without evaluating log q. On the
second put k=−log q, κ=2+λ(k−2), γ=λk and
d=(1−λ)+λ·2√q k. Thus κ∈[2,8], 0≤γ≤8, .3<d<1.5 and q≤e^(−κ).
For a≥0 write

  g(a)=exp(κ(a−1/2))−q exp(−κ(a−1/2)),
  G_*(a)=d sinh(γa)/γ,      G_*(a)=da when γ=0.             (11.3)

Both functions increase strictly; G_* maps [0,∞) onto [0,∞). Set
z₀=g(0)≥0, z_c=g(2/5), z₁=g(1/2)=1−q. Here z_c>.4485,
z_c−z₀>.43 and z_c−g(c)>.26. For example, using q≤e^(−κ), the
first two lower bounds reduce respectively to the minima on κ∈[2,8] of
e^(−κ/10)−e^(−9κ/10) and e^(−κ/10)−e^(−κ/2).
Each derivative has one zero and gives an interior maximum; the smaller
endpoint values are e^(−4/5)−e^(−36/5)>.4485 and
e^(−4/5)−e^(−4)>.43. The same argument for
e^(−κ/10)−e^(−3κ/10) gives z_c−g(c)>.26. These strict margins
ensure that every denominator used next is positive.

For z<z_c put E_t(z)=exp(η/(z_c−z)), A₀=E_t(z₀), E₀=E_t(0),
Λ=exp(A₀−E₀) and

  μ=d E_t′(0)/(g′(0)E_t′(z₀)),
  F(u)=(Λu²+μu)/(1+u),
  M(U)=log(1+F(exp(U)−1)),
  T(z)=z_c−η/log[E₀+M(E_t(z)−A₀)].

For z≥z_c set T(z)=z. On [z₀,z_c), u=exp(U)−1≥0 and
F′(u)=Λ+(μ−Λ)/(1+u)²>0. Therefore T increases from 0 to z_c;
it maps [z₀,z₁] bijectively onto [0,z₁]. Its inverse is elementary:

  F⁻¹(v)=2v/[√((μ−v)²+4Λv)+μ−v],
  M⁻¹(V)=log(1+F⁻¹(exp(V)−1)),
  T⁻¹(w)=z_c−η/log[A₀+M⁻¹(E_t(w)−E₀)]               (0≤w<z_c),

with T⁻¹(w)=w for w≥z_c. All square roots are positive and every log
argument is positive. At z₀, M′(0)=μ gives
T′(z₀)g′(0)=d=G_*′(0). At λ=1, q>0 and κ=k, one has z₀=0,
Λ=μ=1 and G_*=g, so T is identically the identity. The correction is
also flat at z_c: with U=E_t(z)−A₀ and v=exp(−U),

  exp(M(U)−U−logΛ)
    =(1−v)²+(μ/Λ)v(1−v)+v/Λ.                         (11.3a)

The three coefficients multiply nonnegative weights whose sum is one
(type (b)); C29b checks the exact identity. C29a checks F⁻¹ and F′,
and C29c checks the tuned endpoint derivative.
Since U grows like exp(η/(z_c−z)), (11.3a) gives T(z)−z and all its
mixed parameter derivatives as polynomials in U and (z_c−z)⁻¹ times
exp(−U); they vanish at the terminal identity branch. The exact flat-tail
identity is C29b; the vanishing derivatives follow from the displayed
exponential estimate.

One flat radial conjugation makes the inner core exact. For x>0 set
S(x)=log(1+exp(x−1/x)) and S(0)=0. For δ<a<c put

  E_h(a)=exp(c/(c−a)),  φ(a)=E_h(a)−e²,  ψ(a)=E_h(a)−e,
  R(a)=ψ⁻¹(S(φ(a))),  ψ⁻¹(v)=c−c/log(e+v),

and set R(a)=a for a≥c. It is an increasing bijection
(δ,∞)→(0,∞), with R(a)<a on (δ,c), R(a)→0 as a↓δ and
R−id flat at c. Its inverse is the identity above c; below c it is

  R⁻¹(y)=φ⁻¹(S⁻¹(ψ(y))),
  φ⁻¹(v)=c−c/log(e²+v),
  S⁻¹(y)=[w+√(w²+4)]/2,  w=log(exp(y)−1).               (11.3b)

The positive square-root branch is forced by x−1/x=w and x>0. C29g
checks this positive inverse and both logarithmic radial-coordinate
inverses; C29e checks the b_R′ identity below. At zero
the inverse R⁻¹(0)=δ is used only as a one-sided limit in the composite.
Put f₀=T∘g and h=G_*⁻¹∘f₀. The tuning above gives h(0)=0,
h′(0)=1. Define on 0≤a≤1/2

  f(a,q,λ)=G_*(a)                                      (a≤δ),
  f(a,q,λ)=G_*∘R⁻¹∘h∘R(a)                            (a>δ),     (11.4)

and extend it oddly. It is strictly increasing. The two-way formula is

  f⁻¹(y)=G_*⁻¹(y)                                    (0≤y≤G_*(δ)),
  f⁻¹(y)=R⁻¹∘f₀⁻¹∘G_*∘R∘G_*⁻¹(y)                   (y>G_*(δ)),
  f₀⁻¹=g⁻¹∘T⁻¹,
  g⁻¹(z)=1/2+κ⁻¹log[(z+√(z²+4q))/2],
  G_*⁻¹(y)=asinh(γy/d)/γ,  with value y/d at γ=0.

The inverse for negative y is its odd reflection. To verify both
compositions on the outer branch, apply successively G_*⁻¹, R, G_*,
f₀⁻¹ and R⁻¹ to (11.4); the adjacent inverse pairs cancel in exactly
that order. Reversing the order gives the other composition. C29a,d,g
check the nontrivial scalar inverse identities used here. Its arguments stay in the
displayed domains: R(a)≤a≤1/2, while R(u)≤u for u=G_*⁻¹(y), so
G_*(R(u))≤y≤1−q. The γ=0 quotients have smooth Taylor limits;
asinh(v)=log(v+√(1+v²)) makes the inverse finite elementary.
For the inner seam write a=δ+s. Then φ(a)=20e²s+O(s²), whence
R(δ+s)=A(s)exp[−1/(20e²s)] with A smooth and positive. Since
h(x)=x+x²B(x,q,λ), log[h(R(a))/R(a)] is flat in s. The derivative of
log R(δ+s) is 1/(20e²s²)+O(1); the mean-value theorem and repeated
implicit differentiation make R⁻¹hR(a)−a flat with every mixed parameter
derivative. Thus (11.4) glues to G_* to all orders. At c, R−id and
R⁻¹−id are flat by the E_h exponential. These are type-(b)/(c) smoothness
arguments, not numerical seam tests.

The outer height is exactly the exponential on an open interval. Indeed
G_*(c)≤.3 sinh(1.6)/1.6<.446, whereas z_c>.4485; hence
h(2/5)>c. For a≥2/5 both R(a)=a and R⁻¹(h(a))=h(a), while
T(g(a))=g(a). Consequently f(a)=g(a) on [2/5,1/2],
f(1/2)=1−q. At λ=1, h=id and f=g on the entire cell.

For 0<|t_c|≤1/2 put r=|t_c|,
P=log ν+log r−2r²/(1−r²), q=exp(2P), λ=s((r−.3)/.2).
Reduce a=n+b with n=floor(a+1/2), −1/2≤b<1/2, and set

  D=√(f(b,q,λ)²+4q),  Y=(D−f(b,q,λ))/2=2q/(D+f(b,q,λ)),
  ρ₀=(log Y)/2,       ρ(a,t_c)=nP+ρ₀(b,q,λ).             (11.5)

Use the rationalized Y when f≥0. The old outer convex blend is absent:
it would again require an implicit forward inversion. Since f_a>0,
ρ₀,a=−f_a/(2D)<0, and f(±1/2)=±(1−q) gives q≤Y≤1.
On |b−1/2|<.1, Y=q exp[−κ(b−1/2)]; on |b+1/2|<.1,
Y=exp[−κ(b+1/2)]. Thus ρ and every derivative agree across height-cell
seams, and ρ(a+1)=ρ(a)+P. At r=.5, λ=1 and κ=k=−2P, so
ρ=P(a+1/2) exactly; the flat endpoint of λ gives all mixed collar jets.
At r≤.3, λ=0 on an open neighbourhood of the central fibre; no log q is
formed at q=0. Both pieces of the parameter domain agree where they meet.

At t_c=0 the nonnode axis reference is

  a>0: (x,y)=(√f(a,0,0)/ν,0),
  a<0: (x,y)=(0,√f(−a,0,0)).                           (11.6)

The quadratic equation Y²+fY−q=0 gives smooth open axis charts, with
e=exp(−4|t_c|²/(1−|t_c|²)) and D as above:

  a>0:  x=√[(D+f)/(2e)]/ν, y=t_c/x;
  a<0:  y=√[(D−f)/2],      x=t_c/y.                       (11.6a)

For .5≤a<.6 put b=a−1, f_b=f(b,q,λ), D_b=√(f_b²+4q),
and use x=√[2/(e(D_b−f_b))]/ν, y=t_c/x. For −.6<a≤−.5
put b=a+1 and use y=√[2/(D_b+f_b)], x=t_c/y. The exact open
seam formulas, valid on |a−.5|<.1 and |a+.5|<.1 respectively, are

  x=exp[κ(a−.5)/2]/(ν√e), y=t_c/x;
  y=exp[−κ(a+.5)/2],       x=t_c/y.

They also supply transverse derivatives at t_c=0, where κ=2. Near either
node λ=0 and f=a on |a|≤δ, so the forward source height is exactly

  a=e u−v,  u=ν²|x|², v=|y|², e=exp(−4uv/(ν²−uv)).        (11.7)

Away from the node the forward formula is the displayed finite f⁻¹ applied
to eu−v, with q=euv and λ from |t_c|; no quadrature or root solve occurs.
C29d checks the positive quadratic root and the displayed open seam identity.

### 11.3 Inverse-first pole matching

For the elliptic height, and also for (11.7), write the forward local height as
H=m+R_H(u,v), m=u−v, r=u+v. On r≤1/250 the explicit estimates are

  |R_H|≤38r³,  |(R_H)_u|,|(R_H)_v|<300r²,
  |grad_p R_H|≤700r²,  p=Hopf(ν conj x,y).                    (11.8)

Here is the derivative certificate, type (c). With q=euv, 1−e≤221uv,
|e_u|≤222v, |e_v|≤222u, |A_u|≤1.001, |A_v|<.004 for A=eu,
and |q_u|≤1.001v, |q_v|≤1.001u. Differentiate
H_lin−m=[−u(1−e)+q(u−v)]/(1−q). The two derivative bounds, divided by r²,
are at most 112.001/(1−q)+35.035r²/(1−q)² and
223.251/(1−q)+35.035r²/(1−q)², respectively; both are <225.
The arctangent-minus-linear sum contributes <2r² to either derivative, since
Σq^(3j)=1/(1−q³) and Σj q^(3j−1)=q²/(1−q³)². The normalization contributes
<14r², using |1−T′|≤7r² and |dL/du|,|dL/dv|<2. Thus 300 is a common bound.
For (11.7) alone these extra arctangent terms are absent and the derivative bound
follows directly from −u(1−e). Finally u=(|p|+p₃)/2 and v=(|p|−p₃)/2 have
Euclidean gradient norms at most 1, giving the last bound. I2t checks every rational
endpoint comparison.

Let M_j be the previous linear base model and B_j,A its unmodified algebraic
profile map (9.4). The normalized map

  C_j(p)=M_j⁻¹ B_j,A(p+(0,0,R_H(p)))

satisfies |C_j(p)−p|≤20r² and ||DC_j−I||≤40r on r≤1/250: use §9's
bounds 15r² and 30r for M_j⁻¹B_j,A−I, |p+(0,0,R_H)|≤1.001r, and (11.8).
The projected-defect contraction proves existence and uniqueness of C_j⁻¹ on
|p_L|≤10⁻⁴; its solution satisfies

  |p_S−p_L|≤80r_L²,  ||Dp_S−I||≤100r_L,  r_L=|p_L|.          (11.9)

This is a certificate for the explicit inverse chain: inverse stereographic
projection, the quadratic profile inverse (9.3), and (11.5) (or §10.5 in the
elliptic variant). The contraction is used to prove the bounds, not to evaluate
that inverse. `elliptic_pole_raw_inverse(...,elementary=True)` evaluates this chain.
C22 checks its stable north stereographic formula and its Hopf projector identities.
In the elementary cap it simplifies further: f=a, so set R=√(a²+4q),
u=(R+a)/(2e), Y=(R−a)/2 and p_S=(2ν Re t_c,−2ν Im t_c,u−Y).
Use the rationalized expression for the small numerator. This is smooth on both
nonnode axes and avoids dividing |t_c|² by Y. No logarithmic height is evaluated
by this cap implementation.

For z∈C², p_L=Hopf(z), r_L=|z|², put s=|p_S| and Q=Π(p_S/s), where
Π(n)=(I+[[n₃,n₁+in₂],[n₁−in₂,−n₃]])/2. The raw local lift is

  g_j(z)=√s Qz/|Qz|.                                        (11.10)

The denominator never vanishes here. In fact (11.9) gives
|n_S−n_L|≤200r_L, ||Q−Π(n_L)||≤100r_L, and
h=|Qz|²/r_L=1−|n_S−n_L|²/4≥1−10000r_L². In
 g_j=c(I+D)z, D=Q−Π(n_L), c=√(s/r_L)/√h,
the following bounds hold on r_L≤10⁻⁴:

  |c−1|≤42r_L, c≤1.005, ||D_p D||≤300, |D_p c|≤210.

For completeness, differentiating the normalized directions gives
||D_p(n_S−n_L)||≤(100+400+80)/(1−80r_L)<600;
|D_p(s/r_L)|≤400 and |D_p h|≤60000r_L. These imply the displayed bounds on c.
Since ||D_z Hopf||=2√r_L, they give

  |g_j(z)−z|≤200r_L^(3/2),   ||Dg_j−I||≤5000r_L.            (11.11)

The latter estimate has slack: the sum of the three derivative contributions is
<1500r_L. I2u verifies these constants.

An earlier closed pole patch, retained here with its checks and obstruction,
used ε=10⁻⁸ and the map

  ĝ_j(z)=z+s((|z|²−ε)/ε)(g_j(z)−z).                        (11.12)

It is exactly linear for |z|²≤ε and exactly the raw lift for |z|²≥2ε. In its
transition ||Dĝ_j−I||≤13200ε<1/1000. On the entire ball |z|²≤10⁻⁴ it is
≤1/2. Extending the defect by radial projection gives a contraction inverse for
both local maps, proving that the replacement preserves the image of that ball
(the maps agree on its boundary). This proves a smooth injective pole replacement.
The forward inverse of (11.12) is not asserted to be elementary. The two-sided
elementary replacement used below is (11.12h).

There is a finite elementary forward formula on the **raw** part of each cap
(identity (a), C25a; bound (c)). Given a nonzero source spinor
\(w=(\nu_j\overline{x_s},y_s)\), set

  p_S=Hopf(w),  s=|w|²,  p_L=C_j(p_S),  r_L=|p_L|,
  P_L=Π(p_L/r_L),  z_L=√r_L P_Lw/|P_Lw|.                (11.12a)

The formula inverts the raw \(g_j\) whenever \(0<s\leq1/250\) and
\(0<r_L\leq10^{-4}\); the additional condition \(r_L\geq2\varepsilon\)
is needed only for the patches to agree with \(g_j\).
All entries are finite elementary: (11.7) gives the height in these caps,
and (9.4) gives \(C_j\). The bound preceding (11.9) gives
\(|p_L-p_S|\leq20s^2\), hence
\(|p_L/r_L-p_S/s|\leq40s/(1-20s)\leq4/23<2\).
If \(Q_S=\Pi(p_S/s)\), the rank-one identity
\(Q_SP_Lw=h w\), \(h=(1+(p_S/s)\cdot(p_L/r_L))/2\geq525/529\), gives
\(|P_Lw|^2=sh>0\) and shows directly that
\(g_j(z_L)=w\). Since both patches equal \(g_j\) at this \(r_L\), (11.12a)
also inverts them there. On \(|w|^2\leq\varepsilon\), the inverse of either
patch is simply \(w\).
In (9.7) take \(Z_0=\sqrt{s_j}U_jz_L\) and \(h_0=s_jr_L\) to obtain the
round coordinates. On \(1-|X|<\alpha/2\) they agree with the S⁴ component
of the global forward map, by (11.17). The intervening part of the earlier
blend (11.12) has no finite elementary forward inverse supplied here.

The blend cannot be replaced by a smooth diffeomorphic lift obtained merely
by changing the phase over the **unchanged south raw base for (11.7)**.
The identities are type (a), C25b–d, and the bound below is type (c).
Indeed, put
\(p=(p_1,p_2,p_3)=\operatorname{Hopf}(w)\), \(r=|p|\),
\(\rho^2=p_1^2+p_2^2\), \(\kappa=2\pi\),
\(u=(r+p_3)/2\), \(v=(r-p_3)/2\),
\(\beta_2=e^4/(128\pi^2)-2\pi^2\), and \(R_3=-4u^2v\).
Here \(\beta=1+\beta_2\rho^2+O(r^4)\) and
\(R_H=R_3+O(r^5)\). From
\(C_0(p)=\kappa^{-1}Q_0^{-1}B(\kappa\beta p_1,\kappa\beta p_2,
\kappa(p_3+R_H))\), with \(Q_0\) orthogonal, and
\(|B(a,b,c)|^2=2(d-\cos c)/(d+a\cos c+b\sin c)\),
\(d=\sqrt{1+a^2+b^2}\), one obtains

  |C_0(p)|²=r²−κp₁r²+F₄+O(r⁵),
  F₄=2β₂ρ⁴+2p₃R₃+κ²[r²(p₁²−ρ²/2−p₂p₃)−ρ⁴/4−p₃⁴/12],
  |C_0(p)|=r−κp₁r/2+N/(2r)+O(r⁴),  N=F₄−κ²p₁²r²/4. (11.12b)

On the real spinor plane \(w=(x,y)\), we have
\(p_1=2xy\), \(p_2=0\), \(p_3=x^2-y^2\), and \(r=x^2+y^2\).
As a polynomial-divisibility test over \(\mathbb C\), substituting \(x=iy\)
in \(N\) gives \(16y^8\eta\), where
\(\eta=1+e^4/(64\pi^2)-16\pi^2/3<1+55/576-48<-46\), using
\(e^4<55\) and \(\pi>3\) from §4.3. Thus \(x^2+y^2\) does not divide \(N\), so the
homogeneous degree-six term \(N/(2r)\) is not a polynomial. The remainder in
(11.12b) becomes \(O(|w|^8)\); a \(C^6\) function of \(w\) would have a
polynomial sixth-order Taylor term. Consequently
\(w\mapsto|C_0(\operatorname{Hopf}(w))|\) is not \(C^6\) at zero. Any
\(C^6\) lift \(L\) preserving this raw base would satisfy
\(|L(w)|^2=|C_0(\operatorname{Hopf}(w))|\), an impossibility. This only
obstructs the unchanged raw base; it places no restriction on a redesigned
base or on an elementary inverse of the existing blend. In particular, a
smooth diffeomorphic inverse-direction lift preserving the raw base would have
a smooth local forward inverse \(L\), which has just been excluded.

The earlier blended shell (11.12) has an exact **three-real-variable forward reduction**
(identities (a), C26a–c). Given \(w=\widehat g_j(z_L)\), let
\(b=\operatorname{Hopf}(w)\). The unknown raw source orbit is
\(q=\operatorname{Hopf}(g_j(z_L))\); put \(p=C_j(q)\),
\(r_S=|q|\), \(r_L=|p|\), \(n_S=q/r_S\), \(n_L=p/r_L\),
\(h=(1+n_L\cdot n_S)/2\), and
\(\vartheta=s((r_L-\varepsilon)/\varepsilon)\). On
\(\varepsilon<r_L<2\varepsilon\), the sole orbit equation is

  b=Φ_j(q):=(1−ϑ)²p+ϑ²q
    +ϑ(1−ϑ)√(r_L r_S/h)(n_L+n_S).                   (11.12c)

Indeed, with \(P_L=\Pi(n_L)\), \(Q_S=\Pi(n_S)\), and
\(g_j(z_L)=\sqrt{r_S/(r_Lh)}Q_Sz_L\), expand the Hopf matrix of
\((1-\vartheta)z_L+\vartheta g_j(z_L)\) and use
\(P_LQ_S+Q_SP_L=hI+(\Sigma(n_L)+\Sigma(n_S))/2\), where
\(\Sigma(n)=2\Pi(n)-I\). In particular,

  |b|=(1−ϑ)²r_L+ϑ²r_S+2ϑ(1−ϑ)√(r_L r_S h).       (11.12c′)

The phase requires no further solve. Projection gives the positive real multiple
\(P_Lw=[(1-\vartheta)+\vartheta\sqrt{r_Sh/r_L}]z_L\), hence, once \(q\)
is known,

  z_L=√r_L P_Lw/|P_Lw|.                              (11.12d)

Conversely, if \(q\) solves (11.12c), the spinor in (11.12d) has Hopf
orbit \(p=C_j(q)\); its patched image has Hopf vector \(b\) and the same
positive \(P_L\)-projected phase as \(w\), so it equals \(w\).
The margin \(h\geq525/529\) above keeps every denominator nonzero.
The contraction proof for \(\widehat g_j\) gives uniqueness of the orbit
solution in the certified shell cap: two solutions there would give two
preimages of \(w\) by (11.12d). Equations (11.12c)–(11.12d) isolate that earlier
forward problem as inversion of the finite elementary three-variable map
\(\Phi_j\); they do not provide a finite elementary formula for that inverse.
The flat step itself has an elementary inverse (L2a–b); the unresolved part is
its coupling to \(C_j\) and the Hopf projector.

A separate obstruction rules out smoothing the **south canonical raw lift**
by a real scalar taper, even if that taper may change the raw orbit radius
(quadratic-jet identities (a), C27a–c). Put \(\kappa=2\pi\),
\(p=\operatorname{Hopf}(z)\), \(r=|p|=|z|^2\), and write
\(\Sigma(d)=\left(\begin{smallmatrix}d_3&d_1+id_2\\
d_1-id_2&-d_3\end{smallmatrix}\right)\).
The second jet of \(B\) in §4.3, together with
\(\beta-1=O(r^2)\) and \(R_H=O(r^3)\), gives

  C_0(p)=p+Δ₂(p)+O(r³),
  Δ₂(p)=κ((−p₁²+p₂²+p₃²)/2+p₂p₃, −p₁(p₂+p₃), −p₁p₃). (11.12e)

Consequently \(q=C_0^{-1}(p)=p-\Delta_2(p)+O(r^3)\).
In (11.10), \(|q|=r-p\cdot\Delta_2(p)/r+O(r^3)\) and
\(h=1+O(r^2)\). Expanding its normalized projector, the radial and
directional terms cancel as follows, uniformly in the spinor direction:

  g₀(z)=z−Σ(Δ₂(p))z/(2r)+O(|z|⁵).                       (11.12f)

On the real spinor plane \(z=(x,y)\), the first component has

  Im g₀,₁(x,y)=κ x y²(x²−y²)/(x²+y²)+O((x²+y²)^(5/2)).

The homogeneous cubic term is not a polynomial: its numerator at \(x=iy\)
is \(-2i\kappa y^5\ne0\), so \(x^2+y^2\) does not divide it.
A \(C^3\) Taylor term would be a polynomial, hence \(g_0\) is not \(C^3\)
at the pole. For any real scalar \(a(z)\to1\), the imaginary first
component of \(a(z)g_0(z)\) has the same cubic term plus \(o(|z|^3)\);
such a rescaling cannot be \(C^3\). This concerns scalar rescalings of the
canonical raw lift. The preceding C25 obstruction separately concerns
all smooth lifts preserving the unchanged south raw base.

The pole patch can instead have a **finite elementary inverse in both
directions**. The following identities are type (a), C28a–d; the inequalities
and margins are type (b)/(c). Put \(a=\varepsilon=10^{-8}\),
\(c=3\varepsilon/2\), \(b=2\varepsilon\), and \(R_*=10^{-4}\). For \(t<c\)
let \(E(t)=\exp(c/(c-t))\), and set

  φ(r)=E(r)−e³  (a<r<c),   ψ(s)=E(s)−e  (0<s<c),
  T(x)=log(1+exp(x−1/x))  (x>0).                         (11.12g)

Each is strictly increasing onto \((0,\infty)\):
\(T'(x)=\exp(x-1/x)(1+x^{-2})/[1+\exp(x-1/x)]>0\).
Their finite elementary inverses are

  φ⁻¹(x)=c(1−1/log(x+e³)),   ψ⁻¹(y)=c(1−1/log(y+e)),
  T⁻¹(y)=(L+√(L²+4))/2,   L=log(exp(y)−1).

Define \(S(r)=\psi^{-1}(T(\phi(r)))\) for \(a<r<c\), and \(S(r)=r\) for
\(r\geq c\). Its inverse on \(0<s<c\) is
\(S^{-1}(s)=\phi^{-1}(T^{-1}(\psi(s)))\), and it is the identity for
\(s\geq c\). Thus \(S:(a,\infty)\to(0,\infty)\) is an increasing
elementary diffeomorphism. In fact \(S(r)<r\) on \((a,c)\):
\(\phi(r)<\psi(r)\), while \(T(x)<x\) follows from
\(e^{-x}<1/(1+x)\) and \(e^{-1/x}<x/(1+x)\).
For \(a<|z|^2<R_*\), write

  D(z)=√[S(|z|²)/|z|²] z,   D⁻¹(w)=√[S⁻¹(|w|²)/|w|²] w,
  H_j(z)=D⁻¹(g_j(D(z))),  H_j⁻¹(w)=D⁻¹(g_j⁻¹(D(w))).       (11.12h)

Set \(H_j(z)=z\) on \(|z|^2\leq a\), and \(H_j^{-1}(w)=w\) on
\(|w|^2\leq a\). The exterior inverse is asserted on \(g_j(B_{R_*})\),
where \(B_{R_*}=\{z:|z|^2<R_*\}\). Here \(g_j^{-1}\) is exactly the
projector formula (11.12a), without its restriction to the raw collar.
All branches of (11.12h) have a fixed number of elementary operations.

Here are the seam checks. At \(r=a+t\), \(\phi(r)\) has a simple zero
with \(\phi'(a)=6e^3/\varepsilon\),
\(T(x)\) is flat at \(x=0\), and \(\psi^{-1}\) is smooth at zero; hence
\(S(a+t)\) is flat in \(t\). At \(r=c-t\), let \(q=e^{-c/t}\),
\(A=e^3\), \(B=e\), and
\(D_0(q)=B-A-q/(1-Aq)+\log(1+\exp(-q^{-1}+A+q/(1-Aq)))\).
Then \(D_0\) extends smoothly to \(q=0\), and, with
\(\ell=\log(1+qD_0(q))\), direct substitution gives

  S(c−t)−(c−t)=t²ℓ/(c+tℓ).                            (11.12i)

Since \(q\) is flat in \(t\), this joins \(S\) and its inverse to the
identity to all orders at \(c\). To check the inner spinor seam, write
\(z=\sqrt u\,\xi\). On the spherical blowup of the origin, (11.7) gives
\(R_H(un)=u^3\) times a smooth function of \((u,n)\), and (9.4) is analytic.
Thus \(C_j(un)=u\widetilde C_j(u,n)\), with
\(\widetilde C_j(0,n)=n\). The parameterized inverse function theorem gives
\(C_j^{-1}(un)=uM_j(u,n)\), smooth with \(M_j(0,n)=n\). Formula (11.10)
therefore has \(g_j(\sqrt u\,\xi)=\sqrt u\,G_j(u,\xi)\), with \(G_j\)
smooth and \(G_j(0,\xi)=\xi\). Its squared output radius and unit direction
have the forms \(v=u(1+uA_j(u,\xi))\) and \(\eta=\xi+uB_j(u,\xi)\),
with smooth \(A_j,B_j\). Put \(x=\phi(r)\), \(u=S(r)\),
\(x'=T^{-1}(\psi(v))\). Since
\(h_0(t)=(\exp(\psi(t))-1)/t\) is smooth and positive at zero,
\((\exp\psi(v)-1)/(\exp\psi(u)-1)
=(1+uA_j)h_0(v)/h_0(u)=1+u\widetilde D_1(u,\xi)\),
with smooth \(\widetilde D_1\). Taking its logarithm yields
\(F(x')-F(x)=uD_1(u,\xi)\), where \(F(x)=x-1/x\) and \(D_1\) is smooth.
The exact identity
\(x'-x=uD_1xx'/(1+xx')\) gives
\(x'-x=x^2uW(x,u,\xi)\): the resulting equation
\(W[1+x^2(1+xuW)]=D_1(1+xuW)\) has derivative one in \(W\) at
\((x,u)=(0,0)\), so the implicit function theorem gives smooth \(W\).
Consequently both the output-radius change and
\(\eta-\xi\) are flat in \(r-a\). Thus \(H_j\) joins the identity smoothly,
with derivative \(I\), at \(|z|^2=a\). Away from that sphere it is a
composition of diffeomorphisms. It is \(S^1\)-equivariant because \(D\) is
radial and (11.10) is equivariant.

The image and overlap have explicit margins. By (11.11), on
\(|z|\leq\sqrt{R_*}=.01\), the raw defect \(g_j(z)-z\) has Lipschitz
constant at most \(5000R_*=1/2\) and size at most
\(200R_*^{3/2}=.0002\). The projected-defect contraction shows that
\(g_j(B_{R_*})\) contains every \(|w|<.01-.0002=.0098\), in particular
\(B_c\). Since \(D\) maps the annulus \(a<|z|^2<R_*\) onto the punctured
\(B_{R_*}\) and \(D^{-1}\) rearranges only \(B_c\), adding its identity
inner ball proves \(H_j(B_{R_*})=g_j(B_{R_*})\). The inner and outer images
have disjoint squared radii \(\leq a\) and \(>a\), so (11.12h) is the
inverse on that image. Its exterior formula is defined throughout this image:
for \(w\in g_j(B_{R_*})\), either \(|w|^2\geq c\) and \(D(w)=w\), or
\(|w|^2<c\) and \(D(w)\in B_c\subset g_j(B_{R_*})\).
For \(r=|z|^2\geq b\), (11.9) gives

  |g_j(z)|²≥r−80r²≥b−80b²=1.9999968·10⁻⁸>c,         (11.12j)

with margin \(4.999968\cdot10^{-9}\). Both \(D(z)=z\) and
\(D^{-1}(g_j(z))=g_j(z)\) there, so \(H_j=g_j\) exactly for
\(2\varepsilon\leq r<R_*\); the middle comparison uses
\((r-80r^2)'=1-160r>0\) on \(r\leq R_*\). For the modified inner cap,
\(u=S(r)\leq r\) and \(v=|g_j(D(z))|^2\leq u+80u^2\leq r+80r^2\).
If \(v\geq c\), the output radius is \(v\); if \(v<c\), it is
\(S^{-1}(v)<c\). The identity branch satisfies the same estimate, so
\(|H_j(z)|^2\leq\max(c,r+80r^2)\). This replacement changes the raw base inside its
small ball, as required by the C25 obstruction; it does not resolve the
separate inverse of the old global height transition; the candidate (11.4)
supplies that inverse in finite elementary form.

### 11.4 A fixed global projector bound

Here and below the base height is **(11.3)–(11.7)**. Its
two-way scalar implementation and type-(c) constants have been checked as
described below; the end-to-end global inverse is given in §11.6. The old
Gaussian, finite forward and elliptic heights are not used to define this
base. The height's inner cutoff is 1/10; δ in this section denotes only
the target cap thickness. Put

  δ=10⁻⁹, α=10⁻⁵, X₋=−1+δ, X₊=1−δ, L=10⁸⁰⁰⁰.

For n∈S² and X∈[X₋,X₊], evaluate the unmodified inverse base at
(v,X)=((1−X²)n,X). Apply (11.1) to any reference point over that base and set
P(n,X)=e_ref e_ref*. The reference cancels in this formula. In detail, first

  w=O⁻¹v/[1−X√(2−X²)], y=w+(0,1,0),
  u₁=(1−|y|²+2iy₁)/(1+|y|²), u₂=2(y₂+iy₃)/(1+|y|²).

Use (9.3) to invert u₁/u₂=t k_A(|t|), then (11.5) for the cusp radius.
The height is Arg D/(2π) modulo 1, where D=u₂ at zero and D=−conj u₁ at
infinity. Use the primary argument on Re D≥−.9|D| and
a=1/2+Arg(−D)/(2π) otherwise. Their open domains Re D>−.95|D| and
Re D<−.85|D| overlap; the latter gives .4<a<.6 since cos(π/5)<.85.
Thus a central-circle cell seam is interior to a smooth argument chart, with the
neighbouring-cell formulas following (11.6a). On .5≤|t|≤1 use
x₄=Arg(u₂)/(2π)+ξ₄+1/2. On
1≤|t|≤2 use x₄=Arg(u₁)/(2π)−λ_∞δ(t)−(1−λ_∞)/2−1/2, with
λ_∞,δ(t),ξ₄ as in §3. On the central circles use (11.6). In a cusp choose y>0
when y≠0; across its opposite axis choose x>0. These are local reference gauges,
not a global section. For central-fibre derivatives use the open formulas (11.6a),
not a pointwise axis value. The profile inverse is also evaluated without a polar
division: with T equal to the ratio to the dominant S³ coordinate and z=|T|²,
the local coordinate is T√[c/(A(1−z)+√(A²(1−z)²+4Bz))],
where c=2 for t and c=2B for 1/t. This is smooth at T=0.
In a collar convert (x,y,w_c=0) to the middle with Φ/Θ
and the real Cramer formulas. All operations here are elementary. The code is
`elliptic_raw_base_reference(...,elementary=True)` and
`elliptic_meridian_projector(...,elementary=True)`; the retained function names
also support the older alternative.

The following is an explicit, deliberately coarse **type (c) derivative
certificate**:

  ||D_(n,X) P|| < 10⁶⁰⁰⁰ < L.                               (11.13)

The norm on the domain is the product of the round metric and dX², and the
matrix norm is the operator norm. The details are included because a compactness
argument with an unspecified bound would not give a fixed finite expression.
Sixteen numerical derivative samples reached Frobenius norm 32.8697,
leaving a measured margin greater than 10⁶⁰⁰⁰−33 to the displayed
proved bound; this is evidence (d), not part of the certificate.

1. **Target coordinates.** Let d=1−X√(2−X²). Since
   d(1+X√(2−X²))=(1−X²)², d≥δ²/2. The derivatives of
   (1−X²)n and X√(2−X²) have norms at most 2. Thus
   ||Dw||≤12δ⁻⁴ and ||Du||≤24δ⁻⁴<100δ⁻⁴; inverse stereographic
   projection has derivative norm ≤2.

   In fact the raw-base derivative in this first step has a uniform
   bound, independent of δ. Write n′=O⁻¹n and
   β=(1+X√(2−X²))/(1−X²), so w=βn′ and y=βn′+(0,1,0). C52a
   checks β′=(1+β²)/√(2−X²). Inverse stereographic projection is
   conformal with derivative factor 2/D, where
   D=1+|y|²=β²+2+2βn′₂. The X and sphere directions in Dy are
   orthogonal, with lengths β′ and β; also β′>β because
   1+β²≥2β and √(2−X²)≤√2<2. Thus

     ||Du||=2β′/D<6.

   Indeed D≥d_min:=β²−2β+2 and
   3d_min−(1+β²)=2(β−3/2)²+1/2>0, while √(2−X²)≥1.
   This type-(a)/(b) improvement is not used in the much coarser
   final bound (11.13). C52b checks the two sum-of-squares identities,
   and C52c checks conformality of the inverse stereographic map.
   In a 20001-point one-variable sample of the angular worst case
   n′₂=−1, the largest ||Du|| was 3.814212, a measured margin
   greater than 2.185788 below 6 (C52d, type (d)).

   The next, radial-profile inverse also admits moderate uniform
   constants. Put A=(4π)² and B=e⁴, so 144<A<160 and
   50<B<55; the new lower bound follows already from
   e⁴>Σ_(k=0)^7 4^k/k!=16319/315>50. For r≥0, y=r², the
   direct and reciprocal ratio radii are

     R₀(r)=r√[(A+By)/(1+Ay)],
     R_∞(r)=r√[(y+A)/(Ay+B)].

   Their derivatives are the two positive fractions displayed below.
   The dominant-coordinate choice gives |T|≤1. Since R₀(1)>1 and
   R_∞(√8)>1 (the latter is equivalent to 64>B), the inverse
   arguments satisfy y<1 and y<8, respectively. On these intervals
   monotonicity in A and B bounds the derivative numerators below by
   N₀=144+100y+7200y² and N_∞=144y²+100y+7200, and the squared
   denominators above by D₀=(1+160y)³(160+55y) and
   D_∞=(160y+55)³(y+160). Both inequalities

     144N₀²−D₀>0  (0≤y≤1),
     1600N_∞²−D_∞>0  (0≤y≤8)

   are manifest after substituting y=t/(1+t) and y=8t/(1+t):
   after multiplication by (1+t)⁴ the respective coefficient lists,
   in ascending powers of t, are

     (2985824,16013641,317849979,147378331,7082236169),
     (82917380000,348241629000,717234559000,407299703000,74507626600).

   C53a–b check the derivative and polynomial identities. Thus
   R₀′>1/12 and R_∞′>1/40. Their angular scales are R₀/r>1
   (because A+By−1−Ay=(B−1)+(A−B)(1−y)>0 for y≤1) and
   R_∞/r>1/4 (because (y+A)/(Ay+B)>144/1335>1/16).
   For either inverse complex radial map, the differential norm is
   therefore <40. If T is the ratio of the two S³ coordinates, the
   dominant denominator has modulus at least 1/√2 and Cauchy--Schwarz
   gives |dT|≤2|du|. Hence the inverse ratio coordinate `small` in
   the source has ||D small||<40·2·6=480. Where the source uses
   t=1/small instead, |small|≥1/2, so ||Dt||<1920.
   On a cusp branch t_c=small, or t_c=1/small with
   |small|>2, hence ||Dt_c||<480 in either case. The height a is
   Arg of a dominant S³ coordinate (with a constant half-cell shift
   in the other argument chart), whose modulus is at least 1/√2.
   Therefore ||Da||<6√2/(2π)<2. Since ν≤1, the cusp coordinates
   b=(2ν Re t_c,−2ν Im t_c,a) satisfy
   ||Db||<√(960²+2²)<961. This uses π>3 and
   960²+4<961², checked in C53c.
   These are type-(b)/(c) bounds on the ratio/profile stage, not a
   replacement for the height or projector bounds. C53d's separate
   20001-point samples measured minimum radial derivatives
   0.2145580812 (direct) and 0.0318077198 (reciprocal), margins
   greater than 0.1312 and 0.0068 over the proved 1/12 and 1/40.

   On the zero cusp the dominant |u₂|>1/2,
   and on the infinity cusp |u₁|>1/2. On .3≤|t|≤10/3 both are >1/4.
   To check this last bound, use 144<A<160, 1<B<55 and
   144r²/(1+160r²)≤|t k_A(r)|²≤(160+55r²)/144 for .3≤r≤4.
   The cusp bounds follow by the same calculation at r≤.5, or its reciprocal.
   Ratio derivatives cost at most 128||Du||. The radial profile derivative is
   (A+2By+ABy²)/[(1+Ay)^(3/2)√(A+By)]>10⁻⁶ for 0≤y≤16.
   The reciprocal profile derivative is
   (Ay²+2By+AB)/[(Ay+B)^(3/2)√(y+A)]>1/100 for 0≤y≤1/4.
   Angular inverse derivatives are ≤13. Argument derivatives cost at most
   4||Du||. Consequently the local coordinates
   b=(2ν Re t_c,−2ν Im t_c,a) have ||Db||<961, in particular
   ||Db||<10¹²δ⁻⁴ as used below.

2. **Height and source coordinates.** Reduce the representative to |a|≤1/2;
   (11.1) is deck invariant. Write p=q or λ for a parameter derivative on the
   second part of D; on its q=0 part take p=q and λ=0. The following is a
   type-(c) certificate for (11.4). On the collar,
   ν²≥e⁻⁴>1/55, |t_c|²≥9/100 and
   exp[−4|t_c|²/(1−|t_c|²)]≥e^(−4/3)>1/4; hence
   q>9/22000>1/2500, k<8 and .3<d<1.5. Differentiating (11.3),
   with |κ_q|≤2500, |κ_λ|<6, |d_q|<300 and |d_λ|<1, gives
   |g_p|<1800, |g_ap|<20000 and |(z_c)_p|<260.
   Moreover .14<g′(0)<1.5, E_t′(0)/E_t′(z₀)>.29,
   so .05<μ<11; 1≤Λ<1.02. The ratio E_t′ follows from
   E_t′(0)/E_t′(z₀)=[E₀/A₀][(z_c−z₀)/z_c]²,
   E₀/A₀>e^(−.01/.43) and
   (z_c−z₀)/z_c≥1−e^(−.4κ)>.55.

   For u≥0, both F′(u) and F(u)/u lie in [.05,11]; therefore
   1/220<M′(U)<220. Put A=E_t(g(a)), B=E_t(f₀(a)).
   By (11.3a), .74<exp(B−A)<4, or −.31<B−A<1.4.
   To see the lower constant,
   discard the positive μ-term in (11.3a), use 1/Λ>50/51, and minimize
   (1−v)²+(50/51)v; its minimum is 1925/2601>.74.
   The upper bound is 1+10v(1−v)≤3.5<4.
   Since E_t′(z)=E_t(z)(log E_t(z))²/η,
   f₀,a=g′(a)M′(U)[A/B][log A/log B]².
   For A≤2, A,B≥1, log A,log B>.012 and B<3.4, so the bracket
   lies in (2.8·10⁻⁵,7000). For A≥2, the preceding bounds on B−A
   put it in (.18,3). With .14<g′(a)<12 this gives

     10⁻⁸<f₀,a<10⁸.                                    (11.13a)

   Here are parameter bounds, including the moving terminal cut. Directly
   differentiating E_t gives |(A₀)_p|<120, |(E₀)_p|<14,
   |Λ_p|<140 and |μ_p|<2·10⁶. In the last estimate, divide the
   logarithmic derivative of μ into those of d, g′(0), E_t′(0),
   E_t′(z₀); their respective bounds are 1000, 20000/.14,
   1200 and 10000. At A≤20,
   |A_p|≤210000A(log A)²<4·10⁷, |M_p|<5·10⁷,
   |B_p|<9·10⁹, and η/[B(log B)²]<70, hence |f₀,p|<10¹².
   At A≥20, write B=A+Δ, with
   Δ=log(1+A₁e⁻ᵁ+A₂e⁻²ᵁ),
   A₁=(μ+1−2Λ)/Λ and A₂=1−μ/Λ.
   Here |A_i|<20, |(A_i)_p|<3·10⁶, e⁻ᵁ<3e⁻ᴬ,
   |Δ|<80e⁻ᵁ and |Δ_p|<e⁻ᵁ(120|A_p|+1.3·10⁷).
   Since |A_p|≤210000A³ and sup_(A>0)A³e⁻ᴬ=27e⁻³<1.35
   (use e³>20 from Σ_(j=0)^8 3^j/j!=89641/4480>20),
   |Δ_p|<1.5·10⁸ and |A_pΔ|<7·10⁷. Differentiate
   f₀−g=η(1/log A−1/log B), using
   |[1/(B(log B)²)]′|<.001 and 1/[B(log B)²]<.01;
   this gives |f₀,p−g_p|<2·10⁴. Thus |f₀,p|<10¹² on the whole
   cell. On 0≤a≤c, z_c−g(a)>.26 gives E_t′(g)<.16,
   E_t″(g)<1.3, U<.04 and exp(U)−1<.05. Direct product and
   quotient differentiation yields
   |U_a|<2, |U_p|<450, |U_ap|<36000,
   |M_UU|<5·10⁴, |M_p|<1.1·10⁵, |M_Up|<4·10⁶.
   Therefore |B_a|<450, |B_p|<2.2·10⁵,
   |B_ap|<6.2·10⁷. With W(B)=1/[B(log B)²],
   W<7000 and |W′|<1.5·10⁶, so
   |f₀,ap|=η|W B_ap+W′B_aB_p|<2·10¹².
   Since f₀(0,p)=0, |f₀,p(a)|<2·10¹²a on 0≤a≤c.

   The one-soft radial distortion has an independent type-(c) bound.
   On 1/10<a<c write E=E_h(a), L=log E, x=E−e²,
   y=S(x), D=S′(x), F=e+y, M=log F. Then

     R′=DE L²/(F M²),
     b_R=R/R′=cF M(M−1)/(DE L²),
     b_R′=M²+M−1−[F M(M−1)/(ED)](1+2/L)
                     −F M(M−1)D′/D².                 (11.13b)

   C29e checks this exact derivative identity.

   The identities S′=σ(1+x⁻²) and
   S″=σ(1−σ)(1+x⁻²)²−2σx⁻³, σ=(1+e^(−x+1/x))⁻¹,
   give 0<S(x)<x, S/S′≤2x and |SS″/(S′)²|<6. Indeed
   exp(−x)≤1/(1+x) and exp(−1/x)≤x/(1+x) give S<x.
   If x≤1, w=x−1/x≤0, then S≤exp(w), σ≥exp(w)/2,
   so S/S′≤2x²≤2x. If x≥1, σ≥1/2 and S<x give
   S/S′≤2x. Finally
   SS″/(S′)²=S exp(−w)−2(S/σ)/[x³(1+x⁻²)²].
   Its first term is ≤1: for x≤1 use S≤exp(w), and for
   x≥1 use S<x and x exp(1−x)≤1. The second is <4 by
   the preceding S/σ bounds. Thus its absolute value is <5<6.
   For E≥20, x>.6E≥12, E−6<F<E
   (use S>x−1/x and e<11/4 from e⁴<55),
   1≤D≤1+x⁻², |D′|<3/x² and
   0≤L−M<−log(.7)<.4. With v=F/E, rearrange the first two
   terms of b_R′ as

     1+(1−v/D)(M²+M−2)+2(v/D)(M−1)(1−M/L).

   Here D≥1 follows from x−1/x≥2log x for x≥1,
   so exp(−x+1/x)≤x⁻². Their upper bounds are 1, 4 and .8: use
   1−v/D≤6/E+3/E² and the decreasing function
   (log²E+log E)/E on E≥20. The remaining term is
   <(25/3)(log E)²/E<4, also decreasing there. Thus
   |b_R′|<10 on this outer part. The following sharper type-(c)
   estimate closes E<20. C54 checks the rational comparisons and
   finite exponential series used here. For x≤1, S<log 2<3/4,
   F<7/2 and M<13/10, the last inequality following from
   exp(13/10)>7/2. Also S/D≤2x²≤2 and, in the displayed formula
   for SD′/D² above, its first nonnegative term is ≤1 while the
   second is at most 4x/(1+x²)²≤3√3/4<13/10. Their difference
   therefore has absolute value <13/10. Since E>64/9 and L≥2,
   the three terms of (11.13b) have absolute values respectively
   less than

     199/100,   (273/160)·4·9/64,   (273/160)·13/10.

   Their sum is <6. For 1≤x and E≤20, D≥1, F<E, and
   1≤M≤L<3. Set v=F/E. The first two terms of the rearranged
   expression are nonnegative; its second summand is at most
   (L−1)²/(2L)<2/3, by
   4(M−1)(L−M)≤(L−1)². For 1≤x≤3,
   F<6, M<9/5 and the first summand is <101/25.
   The two nonnegative terms in SD′/D² are each ≤1: for the
   second, use S/σ≤2x and 4x²/(1+x²)²≤1. Their difference
   has modulus ≤1, so the last term of (11.13b) has modulus
   <(6·9/5)/(8/3)=81/20. The sum is <10.

   For 3≤x and E≤20, x<13 and
   1−v/D<6/E+1/x²: use x−S<1/x,
   e²−e+1<6 and 1≤D≤1+x⁻². On E≥10,
   (log² E+log E−2)/E<3/5. Indeed
   log 10<7/3, and the derivative of
   (3/5)E−(log² E+log E−2) is
   3/5−(2log E+1)/E>0 for E≥10: the quotient decreases
   since its derivative is (1−2log E)/E²<0, and at E=10
   it is <17/30<3/5.
   Consequently the first summand is <1+18/5+10/x².
   The modulus of SD′/D² is at most
   x exp(−x+1/x)+4x²/(1+x²)²; both terms decrease for x≥3.
   On [3,4], [4,6], [6,13], respectively, the last term
   of (11.13b) is bounded by

     (21/4)·(61/100),   (297/40)·(7/20),   18·(7/50).

   Here F<7,9,16 and M<2,11/5,3 respectively. The exponential
   parts are <1/4, 1/8, 1/30 at the left endpoints; the rational
   parts are at most 9/25, 64/289, 144/1369, with
   1/8+64/289<7/20 and 1/30+144/1369<7/50.
   Adding the first and second summands gives a total <10 on
   all three intervals. The exponential comparisons follow from
   the positive partial sums listed in C54. Hence |b_R′|<10
   globally. C29i's 1999 fixed interior samples have maximum
   |b_R′|=3.301394 at a=.1300000, measured margin 6.698606
   to the proved bound 10; this type-(d) sample is not the proof.
   Since d log b_R/d log R=b_R′, if x=R(a), H=R⁻¹(h(x)),
   then |log[b_R(H)/b_R(a)]|≤10|log[h(x)/x]|, including
   the identity branch above c. Also b_R(1/10)=0 and in fact
   0<b_R(H)<2 when H<c. For E≤20, M−1≤S/e,
   F/E<1, M/L²<3/4, S/D≤2x<26 and e>2 give
   b_R<c·(3/4)·13<2. For E≥20, F/E<1,
   D≥1 and M≤L give b_R<c. Above c, b_R(H)=H<4.

   The identity G_*′(h)=√[d²+(γ f₀)²] gives
   .3<G_*′(h)<9. From (11.13a), h_x and h/x both lie in
   [10⁻⁹,4·10⁸]. Consequently the distortion ratio lies between
   10⁻⁹⁰ and 10⁹⁰. Since
   H_a=h_x(x/h)b_R(H)/b_R(a) and G_*′(H)∈(.3,9),

     10⁻¹⁰⁹<f_a<10¹⁰⁹.                                (11.13c)

   The C23d derivative samples range from .466903 to 6.682791;
   their measured lower and upper margins are .466903−10⁻¹⁰⁹ and
   10¹⁰⁹−6.682791, respectively. These are evidence (d).

   Differentiating G_* by its integral sinh(γs)/γ shows
   |G_*,p(s)|/s<10⁷ for s=h or H. Here s<4 and
   G_*(s)≤1 for both: G_*(h)=f₀≤1−q and
   G_*(H)=f≤1−q. Thus
   cosh(γs)=√[1+(γG_*(s)/d)²]<√[1+(8/.3)²]<30.
   Use |d_p|<300 and |γ_p|≤2500 in that integral derivative.
   The local |f₀,p(x)|<2·10¹²x and global |f₀,p|<10¹²
   give |h_p|/h<10²³. Parameter differentiation of the radial
   composite gives f_p=G_*,p(H)+G_*′(H)b_R(H)h_p/h;
   hence |f_q|,|f_λ|<10²⁶, also on the direct core branch.

   Distinguish Y₀=e^(2ρ₀) from the physical source modulus
   Y_*=e^(2ρ)=|y|². On the reduced representative they now agree
   throughout the cusp. Since q≤Y_*≤1, eu=q/Y_*≤1,
   u<4 and the source radius r=u+Y_*≤5; also e>1/4.
   At t_c=0 use (11.6a), avoiding a quotient q/Y_* on an axis.
   In the core, differentiating f(a,q(u,Y_*))=eu−Y_* with
   |(eu)_u|<1600, |(eu)_(Y_*)|<6300 and
   |q_u|,|q_(Y_*)|<10000 gives |a_u|,|a_(Y_*)|<10¹⁸⁵⁰.
   At fixed t_c, q is fixed and
   a_m=(eu+Y_*)/(r f_a)>10⁻¹⁸²⁰. Thus the core triangular
   inverse has norm <10³⁶⁷². In the collar,
   2/3<|P|<4, |P_r|<7, q>1/2500, |q_r|<2 and |λ_r|≤10.
   For D=√(f²+4q), .04<D<2 and

     10⁻¹⁸²⁰<−ρ_a=f_a/(2D)<10¹⁸²¹,
     |ρ_q|=|(1/Y_*−f_q)/(2D)|<10²⁸,
     |ρ_λ|=|f_λ|/(2D)<10²⁸.

   Hence |ρ_r|<10³⁰, |a_ρ|<10¹⁸²⁰ and |a_r|<10¹⁸⁵⁰.
   In orbit coordinates p=(2ν Re t_c,−2ν Im t_c,m),
   D_transverseρ=p_transverse/(4Y_*r), ρ_m=−1/(2r),
   |D_transverseρ|<1000 and ||D_transverse t_c||<4.
   Thus |D_transverse a|<10¹⁸⁵¹ and
   a_m>10⁻¹⁸²²; the collar triangular inverse is <10³⁶⁷⁵.
   Together the two regions give ||D_p a||<10¹⁸⁵².
   In source coordinates z=(ν conj x,y), ||D_z Hopf||=2√r,
   so ||D_z a||<10¹⁸⁵³√r.

3. **Embedding derivatives.** The first normalized cusp frame has derivative
   norm ≤2/√r. An active shifted frame has a<−.1. In the core,
   f(1/10)=1/10 and f is strictly increasing, so on the active shifted
   support a<−1/10 one has Y_*−eu=|f(a)|>1/10 and
   |y|>1/√10>1/4.
   In the collar |y|>1/55.
   Its unnormalized coordinates are
   (exp(2|t_c|²/(1−|t_c|²))/conj y,
    exp(−2|t_c|²/(1−|t_c|²))conj z₁ y²).
   Their norm is ≥1. Since |D_z |t_c|²|<500, differentiating these
   formulas bounds the normalized derivative by 10⁸. (The largest first
   coordinate term is <2·55·2000+2·55²<10⁶.) Each height weight has
   derivative <14·10¹⁸⁵³√r by step 2. The four-component cusp block
   therefore has derivative <10¹⁸⁵⁶/√r, using r≤5.

   In a collar, |t|∈[.3,10/3], |log t|<6, and a local branch with
   |Arg y|≤π gives |log y|<8. The explicit Φ/Θ and lattice formulas
   have |z_m|,|w_m| and all lattice entries ≤20, with derivatives
   ≤1000√r in source z coordinates; here |D_z t|≤100√r and
   |D_z log y|≤55. These inequalities follow by differentiating
   z_m=−i log y+πs(s−1), w_m=2ts, s=log t/(2πi), or
   z_m=−i log y,w_m=0 at infinity. The lattice determinant is
   (1+|t|²)[log|t|(|t|²−1)+2|t|²]>.1.
   A 3×3 determinant has magnitude ≤6·20³<50000 and derivative
   ≤6·3·20²·1000√r<10⁷√r. Cramer's rule then gives lattice
   coordinates <5·10⁵ and derivatives <6·10¹³√r; x₁ has derivative
   <10¹⁵√r. In an active argument chart θ_j−Im log t is locally
   0 or 2π, so the middle phase derivative is <10¹⁶√r.
   Radial and angular weights cost <10⁴√r. Concatenation and
   V*V≥1/2 therefore give the common bound ||D_z e||<10¹⁸⁵⁶/√r.
   Branch changes are exact phase/atlas identities; these estimates use
   local representatives and do not differentiate a discontinuous selector.

4. **Combine.** The projector kills the vertical circle direction. The
   horizontal right inverse of the Hopf derivative has norm 1/(2√r).
   Combining it with the larger of the two triangular inverse bounds
   in step 2 gives ||D_b P||<10⁵⁵³²/r.
   If r≤1/250, (11.8)–(11.9)'s forward estimate gives |v|<14r;
   since |v|=1−X²≥δ, r≥δ/14. Otherwise r>1/250≥δ/1000.
   Thus always r≥δ/1000. The target-coordinate bound in step 1 is
   <10¹²δ⁻⁴=10⁴⁸, so ||D_(n,X)P||<10⁵⁵⁹²<10⁶⁰⁰⁰<L.
   In the pure middle (.5≤|t|≤2), T=u₁/u₂=t k_A(|t|): the ratio
   bound and the radial/angular inverse-profile bounds in step 1 give
   ||Dt||<128·100·10⁶δ⁻⁴<2·10¹⁰δ⁻⁴. The same Cramer bound applied to
   ξ₄, together with the independent argument coordinate, then gives
   ||Dx₄||<10¹⁵||Dt||+100δ⁻⁴<10²⁷δ⁻⁴, and hence
   ||DP||<10³⁰δ⁻⁴<L. The radial cutoffs are flat at the chart boundaries.

I2v checks the profile lower bounds from step 1, the final powers of ten,
the counts (11.14), and agreement with the implementation: ALL PASS.
C29a–h check the specified exact identities; the intervening differential
inequalities are the type-(c) certificate above, independently audited
line by line. None of its bounds is inferred from the samples in §11.8.

### 11.5 Fixed finite projection and phase formulas

All counts in this section are **fixed integers**, independent of the evaluation
point:

  N=4L+1,              L_g=32(N+1)L,              K=8L_g+1.  (11.14)

They have 8001, 16003 and 16004 decimal digits, respectively. Finite products and sums
below are abbreviations for finite compositions of elementary functions with
these exact counts. There is no limiting operation or unspecified mesh choice.

For Z≠0 put ξ=Z/√(1−X²), n=Hopf(ξ). Let I₀ξ be U₀*ξ in the first pair of
C¹⁰ and I_∞ξ be U_∞*ξ in its fifth and sixth entries, all others zero, where

  U₀=2^(-1/2)[[1,1],[i,−i]],    U_∞=2^(-1/2)[[1,1],[1,−1]].

At h_j=X₋+j(X−X₋)/N, let P_j=P(n,h_j). Define

  e₀=P₀ I₀ξ/|P₀ I₀ξ|,
  e_j=P_j e_(j−1)/|P_j e_(j−1)|, j=1,…,N,   E=e_N.           (11.15)

The initial denominator is positive by the local cap estimate:
||P₀−I₀Π(n)I₀*||≤100r_L, r_L≤2δ/3. Consecutive projectors differ
by <2L/N<1/2, so every subsequent denominator is >1/2. Reversal gives the
exact inverse circle map, as proved in §10.6. The formula is smooth and
S¹-equivariant on the whole slab.

A useful derivative estimate avoids an exponential loss. In a smooth unit
reference ξ over either stereographic disk |a|≤2, ||Dξ||≤10. Hence
||De₀||≤4L. Write, for any unit domain direction,
D e_j=(D P_j)e_j+i A_j e_j with A_j real. Differentiating the normalized
projection and taking its imaginary inner product gives

  A_j=A_(j−1)+Im(e_(j−1)*P_j[(D P_j)+(D P_(j−1))]e_(j−1))
                         /|P_j e_(j−1)|².

Therefore |A_j−A_(j−1)|≤8L and ||DE||≤(5+8N)L.
This follows from P(DP)P=0, P e=e and D|e|²=0; C23c checks the scalar
coefficient cancellation and C23c2 checks the full two-line tangent-jet identity.
Components perpendicular to both lines have zero inner product in this formula.
The phase direction propagates with coefficient one.

For the south band X∈[−1+δ,−1+α] or the north band
X∈[1−α,1−δ], put

  d_j=P(n,X) I_jξ/|P(n,X) I_jξ|,       g_j=d_j*E ∈ S¹.

Here r_L≤2α/3<10⁻⁴, so the cap estimate again makes the denominator nonzero.
The mismatch is independent of the phase of ξ; it is a function of (n,X).
The preceding derivative estimate gives
||Dg_j||≤(9+8N)L<L_g. Thus a real logarithm Θ_j(n,X) is given by **finite**
elementary arguments as follows. Take anchor a₀=−1+α or a_∞=1−α. Start with
the single constant Arg g_j(north,a_j), add

  Σ_(l=1)^K atan(Im R_l/Re R_l),
  R_l=g_j(north,a_j+l(X−a_j)/K)
                       conj(g_j(north,a_j+(l−1)(X−a_j)/K)),

and then add the S² primitive (10.16), using count K at this X. On the south
stereographic chart add its same equator correction P_N(e)−P_S(e).
Every variable ratio has real part >7/8 because its endpoint images differ
by <4L_g/K<1/2. Telescoping proves exp(iΘ_j)=g_j; the connected-overlap
argument of §10.6 proves exact real agreement of its two sphere charts.
The anchor argument is a constant, so its branch cannot create a variable seam.

There is a reduced-count numerical warning (d), not a certificate for the
fixed count (11.14). C55a evaluates g∞ at two directions on
the same latitude, separated by azimuth 1/20, at each of
N=400, 800, 1600. Both real parts are below −.9, while the
imaginary parts have opposite signs and magnitudes above .17;
the sampled projection norms exceed .59. Thus those nearby endpoint
principal arguments lie on opposite sides of their cut. C55b gives
one N=64 meridian with adjacent pure-cusp and pure-middle
embedding supports, whose inner product is exactly zero in the
evaluated branches. These observations concern only the stated
samples and reduced counts; they neither prove a cut crossing nor
nonvanishing for the certified fold. A north-band one-Arg replacement
requires its own global proof.

On the **south band** the phase walk in this formula has an equivalent
single-argument expression:

  Θ₀(n,X)=Arg g₀(n,X),       |Θ₀(n,X)|<10⁻⁶.             (11.15a)

This does not change the defined map. Here is a type-(c) certificate, with
the algebraic identities checked by C51. Put R=1/150000 and ε=100R=1/1500.
For every south-band mesh point, the cap estimate above gives r=|p_L|≤R.
The source is in the unshifted south cusp cell: (11.9) gives
|p_S|≤r+80r²<10⁻⁵; (11.7) gives |a|≤u+v=|p_S|<1/4,
and, since ν₀=1, p_S,1+ip_S,2=2conj(t_c) gives
|t_c|≤|p_S|/2<.3;
the other cusp-frame weight is zero. Thus
P(r)=I₀Π(p_S(r)/|p_S(r)|)I₀* and
Q=I₀Π(p_L/r)I₀* is independent of r. By (11.9),
||P−Q||≤100r≤ε. With p_L=rn and p_S=rn+η, the same bound gives
|η|≤80r², |η′|≤100r and |p_S|≥r(1−80R). Projecting
p_S′=n+η′ perpendicular to p_S/|p_S|, the contribution of n is
at most |η|/r≤80r. The Pauli-matrix identity therefore yields

  ||P′(r)||≤90/(1−80R)<91.

Fix q=I₀ξ, a unit vector spanning Q, and set
f(r)=P(r)q/|P(r)q|. Its denominator is at least √(1−ε²),
|q*f(r)| is that positive real denominator, and
||f(r)−q||²=2(1−|P(r)q|)≤2ε²; in particular
||f(r)−q||<2ε. Differentiating the normalization gives
||f′(r)||≤||P′(r)||/√(1−ε²)<92. The radii of the mesh
increase with X: (11.17a) has
dX/dh=(1+h)/(1+h²)^(3/2)>0 for h=2πr. Hence
Σ_j||f(r_j)−f(r_(j−1))||≤V:=92R=23/37500.

Write e_j=c_j f(r_j) in (11.15), so c₀=1 and
c_j/c_(j−1)=〈f(r_j),f(r_(j−1))〉/
|〈f(r_j),f(r_(j−1))〉|. Since 〈q,f(r)〉 is real,
the imaginary part of each numerator equals that of
〈f(r_j)−q,f(r_(j−1))−f(r_j)〉. Its magnitude is at most
2ε||f(r_j)−f(r_(j−1))||, while the real part is at least
1−V²/2>0. The sum of the principal increments is consequently
bounded by

  2εV/(1−V²/2)=2300/2812499471<10⁻⁶,

with proved margin greater than 1.82·10⁻⁷. This sum is a continuous
real lift of g₀. Its principal argument is the sum itself. The lift
defined above has the same exponential and agrees at the north-point
anchor, so connectedness of S² times the band proves (11.15a).
For comparison only, C47f's two south samples at reduced counts have
|Θ₀|≤3.54·10⁻¹⁸ (type (d)); they are not used in the bound.

Set χ₀=1−s(((X+1)−α/2)/(α/2)) and
χ_∞=1−s(((1−X)−α/2)/(α/2)), using each only in its band, and put

  E_c=exp(−i[χ₀Θ₀+χ_∞Θ_∞])E.                               (11.16)

Where a weight is zero the corresponding phase need not be evaluated. On the
inner half of either band, E_c=d_j exactly. To obtain the W-point, take the
reference from §11.4 and multiply it by the unit scalar e_ref*E_c: in a cusp,
(x,y)↦(conj(λ)x,λy); in the middle, x₁↦x₁+Arg λ/(2π).
Integer changes of the latter are the original circle lattice, so this is an
exactly specified point of W.

### 11.6 Close the poles and obtain the tensor

On the cap of sign j, define the **elementary** local formula

  z_L=U_j*Z /√[s_j(1+|X|√(2−X²))],
  z_S=H_j(z_L),  (x_s,y_s)=(conj(z_S,1)/ν_j,z_S,2),          (11.17)

where (s₀,s_∞)=(2π,π), (ν₀,ν_∞)=(1,e⁻²), and H_j is (11.12h)
with the elementary raw lift (11.10) and inverse height (11.5). The formula includes the exact pole: there
it uses its identity piece and never evaluates n or a vanishing Hopf denominator.
Use (11.17) for 1−|X|<α/2; use (11.15)–(11.16) otherwise.
The cap satisfies |Z|²=1−X²<α and
1+|X|√(2−X²)>2−α/2, so (11.17) gives
r_L<α/[π(2−α/2)]<1.7·10⁻⁶. The patched bound after (11.12j) gives
|z_S|²≤max(c,r_L+80r_L²)<2·10⁻⁶<(.0098)². Thus every source point in this
cap lies in the certified image of H_j, and the elementary H_j⁻¹ of
(11.12h) applies there.

Conversely, for a source spinor w=(ν_j conj(x_s),y_s) with
|w|²<2·10⁻⁶, the inradius in §11.3 puts w in that image. Set
z_L=H_j⁻¹(w), h=s_j|z_L|², and σ₀=−1, σ_∞=+1. The local forward formula is

  X=σ_j(1−h)/√(1+h²),   Z=√[2s_j/(1+h²)] U_j z_L.       (11.17a)

The sphere and denominator identities in C28i show that (11.17a) reverses
(11.17) on its cap. It also gives a source-side selector by computing 1−|X|.

There is an open overlap: take α/3<1−|X|<α/2. If r_L=|z_L|²,
then r_L≤2(1−|X|)/3 and 1−|X|≤7r_L. On this overlap
r_L>α/21>2ε, so H_j=g_j. The homogeneous embedding (11.1) and minimal
Hopf projection (11.10) identify d_j with exactly the same raw g_j lift;
(11.16) equals d_j there. Thus the two W-formulas agree, including all jets.
They cover S⁴. The raw base is globally bijective by the positive scalar
height derivative, the positive radial profile derivative and the clutching
identity C15a. On the free part (11.15) and (11.16) are circle isomorphisms.
Before pole smoothing this is a homeomorphism, smooth with invertible derivative
away from the poles, and equals g_j near them. Replace it on each ball
r_L≤10⁻⁷ by H_j. At this smaller radius ρ=10⁻⁷, (11.11) gives raw-defect
Lipschitz constant at most 5000ρ=0.0005 and size at most
200ρ³ᐟ²<6.33·10⁻⁹. Thus the same contraction proves that
g_j(B_ρ) contains B_c; D only rearranges B_c. The raw and patched maps
agree on the boundary and have the same ball image; (11.12h)–(11.12j)
give an elementary inverse and injectivity. Their
images are disjoint cusp neighbourhoods. The replacement is smooth with
invertible differential, including at z_L=0. Given the new height
certificate, this argument establishes a diffeomorphism G_W:S⁴→W,
without a numerical degree or Nijenhuis assumption.

Restore the torus by §8.5's **explicit τ⁻¹**, using **H_f from (10.3)** for
that trivialization. This is independent of the inverse-first height a used to
parametrize the orbit space. Both heights are deck-equivariant; τ requires only
this property and its smooth central extension. On reduced representatives
A=eν²|x_s|²≤1 and Y=|y_s|²≤1. The additional .5<a<.6 chart has
A≤exp(κ·.1)≤exp(.8)<exp(1)<3 and Y≤1, by the exact open
endpoint expressions following (11.6a); in the core κ=2 and the
stronger A<exp(.2) holds.
On the adjacent −.6<a<−.5 chart, the same endpoint expressions give
A≤1 and Y≤exp(κ·.1)≤exp(.8)<exp(1)<3.
The small patched cap also has A,Y<1. Thus A,Y<3 everywhere used here.
Thus the four-seed formula with its exact q² tail applies everywhere it is
needed, including both axes. On a collar use D₀,D_∞ and Φ/Θ exactly as in
§8.5. C6/C12 prove the gluing for any such height. The inverse-first formula is

  G_cl(Z,X,φ)=τ_f⁻¹(G_W(Z,X),φ).                             (11.18)

The finite forward formula reverses these steps. Given a Y-point y,
first compute its orbit point w=[y] and its torus coordinate φ(y), using
H_f from (10.3) in τ_f. In a cusp use the **node-adapted** representative
before testing a cap. For t_c≠0 set ρ=log|y_s|, P<0,
k_d=floor(ρ/P), h=m(t_c)t_c, and

  (x₀,y₀)=(x_s h^(k_d),y_s h^(−k_d)),
  m(t_c)=ν exp[−2|t_c|²/(1−|t_c|²)],
  s=ν²|x₀|²+|y₀|².                                       (11.18a)

The displayed deck reduction puts q≤|y₀|²≤1. On t_c=0 use the
node-side axis representative: a y-axis point with |y_s|>1 is carried by
ψ̂ to (1/(ν y_s),0), and an x-axis point with ν|x_s|>1 is carried by
ψ̂⁻¹ to (0,1/(ν x_s)); otherwise use the given axis coordinates.
For a point first supplied in a middle chart whose orbit lies in a pole
cap, use Φ or Θ to obtain this cusp representative before the cap test.
Compute the torus coordinate from the same reduction. The four-seed
formula of §10.3 evaluates H_f(x₀,y₀) on its bounded A,Y domain;
deck equivariance gives H_f(x_s,y_s)=H_f(x₀,y₀)+k_d because
H_f(ψ̂p)=H_f(p)+1. Thus in a cusp φ=C−H_f(x_s,y_s)c̃,
with C and c̃ from §8.5. On an axis use the corresponding adjacent
deck integer in this same formula. In a middle chart, perform the same
reduction for each active cusp H_f in D₀ or D_∞ before using
φ=x₂₃−β₀D₀−β_∞D_∞. This keeps every four-seed evaluation in its
stated bounded domain, even for a deck-shifted source chart.
This reduction is essential: an arbitrary U₀ spinor does not give the
source-side cap radius. If s<2·10⁻⁶, let
z_L=H_j⁻¹(ν_j conj(x₀),y₀) and compute the candidate X,Z by (11.17a).
Use that cap formula when 1−|X|<α/2. Otherwise use the free formula
below. Every target point selected by the cap has s<2·10⁻⁶ by the bound
preceding (11.17a). If the candidate fails the cap test, then
1−|X|≥α/2. Since z_L∈B_(R_*), one has
h=s_j|z_L|²<2π·10⁻⁴<1 here, and
1−(1−h)/√(1+h²)≤2h≤2s_j|z_L|²<16|z_L|²,
one has |z_L|²>α/32>2ε; hence H_j=g_j and the raw
orbit formula applies. At a node
the cap formula uses H_j⁻¹(0)=0 and no argument or Hopf division.
If s≥2·10⁻⁶, a point with 1−|X|<δ<α/2 would lie in the
cap image and have s<2·10⁻⁶ by the preceding cap bound. In that
cap representative u=ν²|x₀|², Y=|y₀|² satisfy q=e uY<Y<1 when
t_c≠0, so P<log|y₀|<0 and floor(log|y₀|/P)=0; at t_c=0 both
axis moduli are below one. Thus the source selector uses this same cap
representative. Moreover, in that core
|eν²|x₀|²−|y₀|²|≤s<2·10⁻⁶, λ=0 and f=a on
|a|≤1/10, so the source height has |a|<2·10⁻⁶. The floor
reduction therefore selects exactly this node chart; at t_c=0 its
nonzero axis modulus is <1 and the axis rule selects the same chart.
This contradicts s≥2·10⁻⁶.
Thus every free branch has X∈[X₋,X₊], where its projector is defined.

For the free formula, compute the raw orbit coordinates (v,X) of w.
Here the cusp height is obtained by the displayed f⁻¹ in (11.4)
applied to eν²|x₀|²−|y₀|², followed by the elementary profile and
radial maps (9.1)–(9.4); in a middle chart use its existing elementary
coordinates. Put n=v/(1−X²), choose either explicit Hopf reference
ξ_ref(n) from §8.4, and evaluate E_c(n,X;ξ_ref) by (11.15)–(11.16).
Then

  λ_W=E_c(n,X;ξ_ref)* e(w) ∈S¹,
  F_cl(y)=(√(1−X²) λ_W ξ_ref(n),X,φ(y)).                (11.18b)

Changing the Hopf reference multiplies both E_c and ξ_ref by the same
unit phase, leaving (11.18b) unchanged. On the cap/slab overlap,
E_c=d_j and H_j=g_j, so (11.17a) and (11.18b) agree. Fibre
equivariance and reversal of every finite scalar block give the two
compositions F_cl∘G_cl=id and G_cl∘F_cl=id. The height checks and
derivative review are recorded in §§11.2, 11.4 and 11.8. C20d checks the reverse
projector-phase algebra, C28i checks the cap sphere inverse, and
C29f checks the unit-phase/Hopf forward recovery algebra. The source-side
piece-selector domains are proved just above; reduced-count end-to-end
compositions are a separate numerical safety check.

Every value in (11.18) uses finitely many arithmetic operations, positive square
roots, exp, log, sin, cos and arctan, together with flat-step pieces and atlas
selectors with verified smooth agreement. No theta/elliptic function, series,
quadrature, root solve or limiting transport occurs. The only principal Arg
of a variable outside the right half-plane selects a local chart representative;
its lattice transitions have already been specified.

For a standard six-coordinate sphere/torus chart u as in §8.6, let
h=(h₁,h₂,h₃) be the holomorphic coordinates of (11.18). Define the **finite
coefficient expression**

  M_(2j−1,b)=Re ∂_b h_j,  M_(2j,b)=Im ∂_b h_j,
  J_cl=adj(M) J₀ M/det M,
  J_cl^ambient=E adj(M) J₀ M R/det M.                         (11.19)

The derivative notation here abbreviates finite chain rules, not a differential
equation. Specifically use D(Pe)=(DP)e+PDe,
D(v/|v|)=(Dv−v Re(v*Dv)/|v|²)/|v|,
D atan(b/a)=(a Db−b Da)/(a²+b²), D exp f=exp(f) Df,
D log f=Df/f, D√f=Df/(2√f), the product/quotient rules and
D sin f=cos f Df, D cos f=−sin f Df. In the open flat-step interval,
s′=s(1−s)(u⁻²+(1−u)⁻²); its endpoint derivatives are zero.
For a product use ΣP_N⋯P_(j+1)(DP_j)P_(j−1)⋯P₀; sums have the
fixed counts (11.14). On an argument chart D Arg z=Im(dz/z). Differentiate
one smooth local representative; at a selector seam use the overlapping chart.
Write every squared complex modulus as (Re z)²+(Im z)² before differentiating.
At a central fibre use (11.6a) and the smooth profile inverse just displayed;
at a pole use the identity piece of (11.12h). No derivative of a polar 0/0 is needed.
For other phase representatives (including the final λ), use the same two open
argument charts centred at 0 and π, with the elementary formulas
2 atan(Im z/(|z|+Re z)) and π+2 atan(−Im z/(|z|−Re z)).
Their respective denominators are positive on the stated open domains, including
Re z=0 in the first chart. Their integer differences are holomorphic
lattice transitions. Thus the derivative rules act on an open source chart at
every point, including an argument-cut point; they never differentiate the jump
of a principal-argument selector. The numerical derivative helper is restricted
to interiors of such source charts.
These rules eliminate every derivative into the same finite elementary class.

Equivalently, (11.19) can be read **without a differentiation operator**. Replace
each scalar by a seven-tuple (v,w), where w∈R⁶, seed u_b by (u_b,e_b), and seed
constants by (c,0). Evaluate the finite expression (11.18) with these rules:

  (v,w)+(a,b)=(v+a,w+b),
  (v,w)(a,b)=(va,aw+vb),
  1/(v,w)=(1/v,−w/v²),   √(v,w)=(√v,w/(2√v)),
  exp(v,w)=(exp v,exp(v)w),  log(v,w)=(log v,w/v),
  sin(v,w)=(sin v,cos(v)w), cos(v,w)=(cos v,−sin(v)w),
  atan(v,w)=(atan v,w/(1+v²)).

Use the same rules componentwise for complex numbers, and conjugate both entries
when taking a conjugate. On an argument chart the rule is
(z,w)↦(Arg z,Im(conj(z)w)/|z|²). Use the displayed s′ rule for the flat step,
with coefficient vector zero on its constant pieces. For a nonzero real value,
the absolute-value rule is (v,w)↦(|v|,sign(v)w). At zero, simplify the surrounding
constant plateau first; use f=a on |a|≤.1. At a height-cell seam use its open
κ-dependent endpoint formula above, rather than differentiating
the floor selector. Choose a local chart using
the value components; its integers stay constant in that open chart. The six
coefficient components of the three complex outputs are precisely the rows of M
in (11.19). C24 verifies these arithmetic rules. Thus this is a finite scalar
coefficient prescription with explicit constants and finite elementary operations,
not merely a formula that assumes an unspecified coordinate Jacobian is known.

Under the new-height certificate, G_cl is a diffeomorphism, so det M≠0. Its pulled-back holomorphic charts
are a complex atlas on the whole standard round S⁴×T². Therefore J_cl is
smooth, J_cl²=−Id and N_(J_cl)=0. In the ambient convention,
(J_cl^ambient)²=−P_T and J_cl^ambient P_T=P_T J_cl^ambient=J_cl^ambient.
C9/C10 check these matrix and naturality identities generically; their essential
invertibility and gluing hypotheses are established for the new height above.
No manuscript input (I1)–(I3) or 2026 preprint is used. This is a new parametrization of the
paper's same holomorphic Y, hence defines that complex manifold on S⁴×T²;
it is not asserted to give identical coefficients to the legacy fitted map.

**Expression-class distinction.** (11.19) is a finite elementary, globally glued
closed expression for J under the new-height certificate above. It is not
a compact or feasible-to-expand expression: its fixed sums/products have the
enormous counts (11.14). Formula (11.18b), together with the cap branch,
gives a finite elementary inverse F_cl in the same expression class. The
old inverse-first implementation and its tests alone establish neither the new
height nor this stronger inverse pair; the new identities and type-(c)
bounds above supply the certificate. The statement requires human review.

### 11.7 Implementation and verification status

`cusp.elementary_closed_parameters()` returns all constants and counts as exact
rational numbers/integers. `elementary_projection_section`, `elementary_phase_family`,
`elementary_global_slice`, `elementary_pole_patch`, `elementary_pole_forward`,
`scalar_finite_height`, `elementary_restore_torus` and
`elementary_global_coordinates` now implement the terminal-cut height
(11.3)–(11.4), the pole patch (11.12h), and (11.15)–(11.18)
independently of the legacy numerical evaluator. The new
`elementary_forward_moment` and `elementary_cusp_height_chart` evaluate the
finite forward height, and `elementary_global_forward` evaluates (11.18b)
with the cap selector. The defaults use the revised counts (11.14).
The 38-case reduced-count two-way implementation check passes, including
four exact half-cell seams compared modulo the middle lattice.
An explicit
`max_work` guard prevents an accidental astronomical calculation; the parameter
specification and pole formula can be evaluated immediately. Supplying smaller
counts defines different finite maps for numerical checks. Such checks do not
certify those smaller counts globally and are not described as evaluation of the
full certified expression. The code checks actual nonzero/right-half-plane
conditions on each evaluated path. The scalar `mpmath` evaluator also takes
guarded identity shortcuts near the flat terminal and radial seams; its
numerical output is at working precision, while (11.3)–(11.4) define the
exact symbolic map. The `mpmath` pole evaluator has two precision
shortcuts, justified by the following type (c) bounds; the exact symbolic
definition remains (11.12h). Write \(a=10^{-8}\), \(c=3a/2\), and let
\(d=\texttt{dps}\geq1\). At the inner guard put \(r=|z|^2\),
\(x=\phi(r)\), \(u=S(r)\), and \(F(x)=x-x^{-1}\). Its condition is
\(F(x)\leq-(d+10)\log 10\). Since \(T(x)\leq e^{F(x)}\) and
\(\psi(u)\geq(e/c)u\), it gives \(u\leq(c/e)10^{-d-10}\). For
\(0<t<c/10\), \(\psi'(t)<5/c\), \(\psi(t)>2t/c\), and \(\psi(t)<2\):
use \(e^{10/9}<4\), \((c-t)^2>81c^2/100\), and convexity of \(E\).
Since \(1-e^{-\psi}\geq\psi/(1+\psi)\),
\(J(t)=\log(e^{\psi(t)}-1)\) has
\(J'(t)=\psi'(t)/(1-e^{-\psi(t)})<15/(2t)\). Also
\(J(S(r))=F(\phi(r))\), \((F^{-1})'<1\), and
\((\phi^{-1})'\leq c/(9e^3)<a/120\).

The squared radii \(u,v,s\) used in these derivative comparisons are
\(<c/10\): the raw forward radius is
\(v\leq u(1+80u)<2u\), while at the inverse guard
\(s=S(|w|^2)\leq(c/e)10^{-d-10}\) and its raw inverse radius is \(<2s\).
For the raw output squared radius \(v=|g_j(Dz)|^2\), (11.9) gives
\(|v-u|\leq80u^2\), hence \(|\log(v/u)|<160u\). The preceding derivatives
give \(|S^{-1}(v)-r|<10au\). Normalizing the two spinors and using
(11.11), \(v>u/2\), gives a unit-direction change less than
\((200+80)\sqrt2\,u<400u\). Since \(r\geq a\), these bounds imply
\(|H_j(z)-z|/|z|<1000u\). For the inverse shortcut, put
\(s=S(|w|^2)\) and \(u=|g_j^{-1}(Dw)|^2\). The raw preimage lies in
\(B_{R_*}\), so (11.9) and (11.11) give
\((1-200R_*)^2u\leq s\leq(1+80R_*)u\), hence \(s/2<u<2s\) and
\(|s-u|\leq80u^2\leq320s^2\). Therefore
\(|\log(u/s)|<640s\), \(|S^{-1}(u)-|w|^2|<40as\), and the unit-direction
change is less than \((400+320\sqrt2)s<1000s\). It follows that
\(|H_j^{-1}(w)-w|/|w|<3000s\). Both relative errors are less than
\(3000(c/e)10^{-d-10}<10^{-d}\).

At the outer shortcut let \(q=e^{-c/(c-r)}\), with \(r\) the argument's
squared radius. The implemented condition gives
\(q<10^{-\texttt{working dps}-10}\). Here \(r>c/2\). Since
\(20<e^3<21\), \(2<e<3\), \(q<10^{-11}\), and the exponential in
\(D_0\) is less than one, its three terms have absolute bounds \(19,1,1\);
thus \(|D_0(q)|<22\), \(|\ell|<44q\), and \(c+t\ell>c/2\). Hence
\(|S(r)-r|<88cq\), a relative error less than \(176q\). For
\(S^{-1}(r)\), apply the same bound at \(r'=S^{-1}(r)>r\):
\(q_{r'}<q_r\) and \(|r'-r|/r<176q_r\). Each is below working precision.

For source comparisons that must follow the **exact formula branches**,
evaluate under the context manager
\(\texttt{with cusp.elementary\_formula\_branches():}\).
It bypasses the six working-precision identity returns in the pole radial,
pole patch/forward, radial height and terminal height evaluators. The
piecewise identity formulas on their actual constant regions and the
algebraically equivalent stable formulas remain in force. Check (m1) in
\(\texttt{symbolic\_identities.py}\) matches all six guards to the source
syntax and counts them; (m2) verifies nested use and restoration after an
exception. Check (m3) reproduces the two strict-interior observations
below (sample size two). This opt-in mode still computes with finite-precision
\(\texttt{mpmath}\), and the production work guard still prevents literal
evaluation of the huge fixed counts. It supplies a branch-faithful numerical
comparison path for the symbolic circuits, not an exact arithmetic engine or
a global second-jet certificate. As type-(d) interior observations, at
\(\texttt{dps}=50\) the ordinary evaluator returns the identity for
\(r=(3\cdot10^{-8}/2)(199/200)\) and
\(v=(1/5)(199/200)\); formula-branch evaluation at 240 digits gives
\(r-S(r)=9.0129315017\cdot10^{-99}>0\) and
\(v-R(v)=3.23193414504\cdot10^{-92}>0\).
The strict inequalities \(S(r)<r\) on \((a,c)\) in §11.3 and
\(R(v)<v\) on \((\delta,c)\) in §11.2 already prove that these are
different functions throughout those open regions. Differentiating
the ordinary shortcut there would differentiate the identity instead.

The global tensor and inverse-pair statement require human review.
The checks include the new height, individual exact identities, and
reduced-count end-to-end F_cl compositions. Neither huge
fixed-count product has been numerically evaluated.

### 11.8 Checks, margins and scope of the result

The following table separates current checks from earlier historical
measurements. C29a–h pass for the individual new algebraic identities
specified in §§11.2, 11.4 and 11.6. They do not certify the parameter
bounds or the complete global F_cl composition. Sampled values and
derivatives are numerical observations, type (d).

| Check | Result and margin |
|---|---|
| C21–C22 | deck/branch identities, norm square, north stereographic inverse, Hopf projector and horizontal derivative identities: ALL PASS |
| C23–C24 / I2t–v | quadratic height, projection-phase and finite coefficient rules, specified endpoint arithmetic and implemented counts: ALL PASS |
| C25a–d | raw-cap forward projection, south raw-base norm expansion and sixth-order divisibility obstruction: ALL PASS; the explicit obstruction coefficient η<−46 |
| C25e, numerical | 18 independent raw-cap forward compositions on both caps; maximum relative error 6.3·10⁻⁶⁸ |
| C26a–c, earlier blend | projector anticommutator, Hopf vector of the earlier shell blend, phase projection and blended norm: ALL PASS |
| C26d, numerical, earlier blend | 30 two-cap earlier-shell endpoint/interior base, norm and phase compositions; maximum relative residual 4.1·10⁻⁶⁷ |
| C27a–c | south raw-base quadratic jet, normalized-projector cubic jet and nonpolynomial real-plane component: ALL PASS |
| C28a–d,h–i | elementary radial inverse, seam identities, raw-collar and cap-domain margins, rational north stereographic inversion and forward-cap identity: ALL PASS |
| C29a–h, new height | quadratic F inverse/derivative, convex tail, tuned endpoint derivative, positive quadratic root/open seam, one-soft b_R′ and its inverse coordinates, unit-phase/Hopf recovery, and cusp torus reversal: ALL PASS |
| C29i, numerical | 1999 fixed radial-derivative samples; max \|b_R′\|=3.301394, measured margin 6.698606 to the proved bound 10 |
| C54a–d, sharper radial distortion (a–c) | Exact one-soft algebra, monotonicity derivatives, positive finite exponential series, and rational interval sums certify \|b_R′\|<10: ALL PASS |
| C55a–b, reduced-count north observations (d) | Six fixed endpoint samples at N=400, 800, 1600 straddle the principal-Arg cut with projection norms above .59; one N=64 sample has adjacent disjoint cusp/middle supports and zero overlap. These are numerical warnings, not claims about the certified count: ALL PASS |
| C56–C58, direct middle lift (a,b,c,d) | Section 12 specifies two direct local matrices, positive norm-square bounds >1/480000 and >1/80000, and a phase-corrected blend over the whole pure middle with norm square ≥1/2 and explicit adjoint inverse. The overlap sign follows from 16 exact positive-coefficient polynomial identities. Separately labelled source and seam samples pass. This is a regional replacement, with cusp/pole extension and the full J circuit still open. |
| C59, direct cusp lifts (a,b,c,d) | Section 12.4 supplies separate circle lifts and adjoint inverses on both punctured cusps, including the nonnode central fibres. Signed Hopf components prove nonvanishing; rotating height frames give overlap real part >1/20 and blend norm square ≥1/2. The 1232 source samples and 96 height-seam pairs pass. Radial gluing, pole matching and the global tensor remain open. |
| C60, pole matching (a,c,d) | The C59 first matrix is changed only inside the existing small pole bands. The matrix norm square stays >0.98 and its inner formula is exactly the existing raw pole lift, so the explicit cap replacement of §11.6 applies. Forty source samples have minimum norm square >0.9999999985. This closes the regional pole-matching task; radial gluing and the tensor circuit remain open. |
| C61, radial overlaps and remaining phase (a,b,c,d) | The middle phase offset is extended through both source collars, with the north integer chosen to match at r=2. The existing coefficient certificate then applies on 49/100≤r≤100/49. The 2048 source tests and 32 seam pairs pass. A separate 32×32 phase grid suggests Re(δ_n/δ_s)>0; this sign is not proved and is the explicit next gluing gate in §12.6. |
| C62–C64, practical radial phase (a,b,c,d) | At the different exact boundary radii 1/2 and 2, cancel the source embedding to two rows. The written derivative bounds are \(|g_θ|<1900\), \(|g_φ|<2700\); two explicit coordinate-circle arguments prove zero winding. A 65536-term phase sum has every adjacent real part >31/32 and is smooth and periodic. C62 checks 720 extension points, 60 seams and 144 scalar/matrix comparisons; C64 evaluates the full certified count. This proves the alternative gluing in §§12.7–12.9 without assuming C61's global sign. |
| C65–C66, global value prototype and practical tensor test (a,d) | The assembled new map and inverse pass 65 round trips, including both caps and nonnode central fibres. C66a–b prove endpoint formulas for the phase's first/second derivatives. One actual tensor in the radial interpolation has square/N residuals <6·10⁻⁵⁸ and two independent finite-difference checks of dJ. This first evaluator uses numerical differentiation elsewhere; C67–C68 supply the finite circuit and source binding. |
| C67–C68, finite tensor jets and source binding (a,d) | A six-variable finite second-jet algebra executes the same value source, using (12.36) for the phase jets and polynomial modulus squares at zero. All unary rules are differentiated symbolically; actual product/unary bodies are bound to the checked formulas. Twenty complete definition digests, all 56 new geometric tests and seven vectorized selectors are pinned. Regional tensor points include both nodes, central fibres, caps and radial joins. The written global certificate and numerical limitations are in §12.11, for human review. |
| C30a–n, exact local second jets (a) | elementary, projection and phase second-jet identities; source checks that the certified counts reach the two recurrence sites; both cusp Cramer/torus-restoration signs and one open noncap branch of the four-seed height; inverse-map-convention \(J^2\) (36 entries) and \(N_J\) (216 entries): ALL PASS |
| C31a–d, one pure-middle projector cell (a,c) | source-linked reciprocal dominant branch with exact witness \(t=(-1+2i)/3\), explicit rank-one projector and complete local second-jet chain: ALL PASS; the rational witness bounds are written in the script |
| C32a–d, one two-active-pair cusp projector cell (a,b,c) | on \(1/10<H<1/4\), \(|xy|<3/10\), both cusp signs: the source frame and nonconstant flat weight give an explicit rank-one four-slot projector; all 80 first and 400 second source-coordinate entries, including mixed-partial symmetry, are checked: ALL PASS; \(t=1/10,H=1/5\) is an exact nonempty-cell witness |
| C33a–d, negative-height cusp projector cell (a,b) | on \(-1/4<H<-1/10\), \(|xy|<3/10\), \(y\ne0\), both cusp signs: the source deck-shifted frame has the same product \(xy\), its four-slot projector is unit rank one, and the source weight/deck exponential have exact two-jets; all 80 first and 400 second entries in the five independent variables \((\Re x,\Im x,\Re y,\Im y,H)\) obey the checked product rule: ALL PASS. A strict south round-meridian witness has \(t=1/10,H=-1/6\). |
| C34a–f, negative-real-log active-collar seam (a,b,c,d) | on \(3/10<|t_c|<1/2\), with south \(t=t_c\) and north \(t=1/t_c\): exact source Cramer solves reconstruct \(e^\rho\) in each cusp slot; deck changes preserve the surviving middle angular phase. The surviving angular entry and reconstructed cusp input have matching value/first/second jets in \((\Re t,\Im t,\rho,H,\mathrm{one},\mathrm{four})\); strict target witnesses and angular perturbations show nonempty seam sectors: ALL PASS. Manifest determinant bounds are \(>981/5000\) south and \(>40\) north; the separately labelled two numerical samples are 1.2640336891 and 125.5013159820. The source formulas imply projector agreement on this seam; they do not differentiate the earlier round-meridian map or cover other chart seams. |
| C35a–c, active upper-negative-axis collar projector (a,b,d) | On \(3/10<|t_c|<1/2\), \(3\pi/4<\Arg t<\pi\), and \(|H|<1/10\), the ordered source branches give a three-slot rank-one frame in the ten-slot embedding, supported on \((0,1,9)\) south and \((4,5,9)\) north. Its two positive normalization factors satisfy \(R_f>0\) and \(D=(1-s)^2+s^2\ge1/2\). Exact round-meridian witnesses and the five-variable value/first/second jet product rules pass: ALL PASS. Two separately labelled one-point numerical observations give \(D=0.9773063821\) south and \(0.6619583826\) north, margins \(0.4773063821\) and \(0.1619583826\) to the proved \(1/2\) bound. This is one projector cell, not the full meridian or free-slab tensor. |
| C36a–c, two active-collar free-slab determinant cells (a,b) | At the two strict C35 round witnesses, the source path from \((n_0,n_1,X,\theta,\phi_2,\phi_3)\) to the six real middle-\(W\) coordinates has a factored nonzero determinant: ALL PASS. The raw radial, ratio/profile, height/Cramer, finite-fold fibre, and torus-restoration factors are checked against the production branches. In the fixed target order \((\Re t,\Im t,\Re z,\Im z,\Re w,\Im w)\), the determinant is strictly negative in both cells; the adapted base determinant is positive. The manifest factor \(1-s^2\) is \(>15/16\) for \(9/100<s=|t_c|^2<1/4\). This uses the existing finite-projection nonvanishing bound and proves local invertibility only in these two neighborhoods; it does not supply their full Hessians or \(N_J\). |
| C37a–d, one all-mesh pure-cusp free-slab cell (a,b,c,d) | On the strict near-south chamber \(X<0,\ X>-1+10^{-9},\ 1+X>10^{-5},\ \Lambda(X)<10^{-4},\ |n_0|<1/100,\ n_1<-3/4,\ 1/2<n_2<2/3\), every source meridian sample, including the endpoints of the certified \(N=4\cdot10^{8000}+1\) fold, takes the same positive-\(x\) south-cusp branch. The fixed-count projection overlaps are uniformly \(>3/5\), the four-seed height is \(H_f=a/(1-q)\), and the adapted six-real Jacobian has the positive factorization below: ALL PASS. The exact round witness is \(Z=(1/100,i/200),\ X=-\sqrt{7999/8000}\). A separately labelled numerical initial-overlap sample is \(0.9999999999999999999512\), with measured excess \(0.3999999999999999999512\) over \(3/5\); an \(N=1\) endpoint sample is not a minimum over the certified mesh. This cell supplies neither the six Hessians nor a named \(J\). |
| N8, actual inner-pole charts (a) | both signs: the six real outputs of (11.17)–(11.18), including the four-seed \(H_f\), are differentiated to a 6×6 Jacobian and six Hessians; all 90 mixed-partial pairs per sign, the nonzero pole determinant, and 216 normal-frame \(N_J\) entries: ALL PASS |
| C38a–e, pure-cusp argument-selector overlap (a,b,c,d) | At both cusps, exact round witnesses and \(\pm1/1000\) angular perturbations lie on opposite sides of the source selector inside valid open cusp charts. In the lower wedge, the principal and half-centred heights differ by one deck unit; the two normalized ten-slot references, their projectors, and all first/second jets agree openly. The four-seed height shifts by one, and torus restoration gives a holomorphic cusp deck transition with complex determinant \(1\): ALL PASS. The exact bound is \(|t_c|<1/10\); separate numerical lower-wedge samples give \(|t_c|\approx0.0039838456\) south and \(0.0294368583\) north, measured gaps \(0.0960161544\) and \(0.0705631417\). Fold propagation assumes the already certified nonzero projection overlaps; other selectors and global \(J\) coverage remain open. |
| C39a–e, dominant-coordinate selector overlap (a,b,c,d) | On \(99/100<|u_1/u_2|^2<101/100,\ |\Arg(u_1/u_2)|<\pi/5,\ |\Arg u_2|<\pi/5\), both source positive-root formulas give the same effective \(t\) and all its real first and second jets. Explicit profile bounds place \(3/10<|t|<1/2\), \(1/100<q<1/8\), with Cramer determinant \(>981/5000\), torus determinant \(>15/16\), and raw-reference norm square \(\ge1/2\). Both branches select the same south active collar, middle data, ten-slot reference, projector and restored output on this open cell: ALL PASS. Exact round witnesses at \(\lambda=2-\sqrt2\) and \(\lambda\pm1/10000\) meet the selector and both sides. C39e is a separate one-point type-(d) source observation at the seam; measured gaps to the displayed bounds are \(0.0623622466\) and \(0.0526377534\) for \(q\), \(1.0927862054\) for Cramer, \(0.0441843611\) for torus, and \(0.1888800593\) for reference norm square. Other selector seams and global \(J\) coverage remain open. |
| C40a–f, north infinity/middle flat selector (a,b,c,d) | At the exact round witness \(|\mathrm{small}|=1/2,\ a=3/8,\ \Arg t=\pi/3\), and in a small open angular/height neighborhood, the source height and sparse projector have matching complete second jets across the selector. Conditional on matching finite-fold endpoint jets, a nonzero final phase, and nonzero fold denominators, the phase-corrected and restored middle formulas match through second order modulo local integer deck. C40a–e check source guards, exact Cramer and terminal-height identities, flat cutoff derivatives, all mixed projector jets, and phase cancellation: ALL PASS. Exact round witnesses at \(s=499/1000,1/2,501/1000\) lie on and strictly across the selector. Proven local bounds are \(q>1/2025\), Cramer \(>30\), torus \(>9/10\), reference norm square \(>1/4\); C40f's separate two-point type-(d) source sample has measured minimum gaps \(0.000709386291679574016\), \(19.9925817884944262\), \(0.036998497999\), and \(0.783279071736023484\), respectively, with a nonzero left interior height difference \(\Delta=-1.37624625388097242\cdot10^{-87}\) at 190 digits. Intermediate projector joins remain open; this is one seam chamber, not global \(J\) coverage. |
| C41a–f, cusp/middle \(|t_c|=3/10\) join (a,b,c,d) | Along the whole seam, at both cusp signs, for any argument chart of \(a\) and any fixed branch of \(\log t\): the core reference is \(g\) times the \(\lambda=0\) middle continuation \(M_0\), with \(g=t_c/|t_c|\) for \(a>0\) and \(g=1\) for \(a\le0\), so the ten-slot projectors are identical. The right-hand source differs from \(M_0\) only by terms with a factor flat at \(|t_c|=3/10\). For every \(\phi\) and a common final unit vector, the c0/ci restoration equals \(\Phi\)/\(\Theta\) of the middle restoration, because \(D-cc=H_fc\); the complex determinants are \(i\) and \(-i/t^4\). All parts: ALL PASS. Also checked: the source formulas, the paper's chart maps in `tfull.py`, the band \(2/25<|t_c|^2<1/10\) (where \(1/2500<q<1/8\), the Cramer determinant is \(>0\) and \(|a|<3/5\)), and exact witnesses at \(|t_c|=299/1000,\,3/10,\,301/1000\). C41f records 18 production observations (d) at 120 digits, over all cells and both argument charts: projector difference \(\le2.4\cdot10^{-121}\), restored-output difference \(\le3.9\cdot10^{-120}\), and right-side flat weight \(3.8\cdot10^{-87}\). Propagation through the fold assumes nonzero projection denominators and checks of the other seams. |
| C42a–f, all pure-middle seams (a,b,d) | On every open subset of \(1/2<|t|<2\), modulo the group generated by \(\mathrm{one}\mapsto\mathrm{one}+1\), \(\mathrm{four}\mapsto\mathrm{four}+1\) and \((\mathrm{one},L)\mapsto(\mathrm{one}-\mathrm{four},L+2\pi i)\), the source \(W\)-data equal data built from continuous branches of \(\arg u_1\), \(\arg u_2\), \(\arg t\) and \(\log t\). The ten-slot projector is that of the continuous data. For a final fold vector \(E\) that is smooth there (a premise resting on the other seams and on nonzero overlaps), \(\Psi(W)\) changes by lattice vectors \((m-kn)e_1+ne_4\), so the restored output changes by holomorphic deck translations (complex determinant \(1\)). This covers the \(|t|=1\) selector (the continuous formulas differ by an integer) and the cuts of \(u_1\), \(u_2\) and \(t\), with their corners: ALL PASS. Also checked: the pure middle is the inverse branch (\(|\mathrm{small}|\le e^{-1}\) in the direct branch), all cusp weights vanish there, and exact witnesses lie at \(|t|=999/1000,\,1,\,1001/1000\) and at the \(t\)-cut corner. C42f, numerical (d), sample size 32: across eight seams (two in the \(\lambda\) band) the jump of \(\Psi\) is an integer lattice vector plus a residual of order \(\varepsilon\) (\(10^{-12}\to10^{-24}\)). The collar joins at \(|t|=1/2,\,2\) are separate. |
| C43a–f, argument selector for \(|t_c|<1/2\) (a,b,c,d) | On the lower overlap \(-\pi<\Arg p<-5\pi/6\), where the half-centred chart has \(a_1=a_0+1\), at both cusp signs and for every \(|t_c|<1/2\): the two charts give the same ten-slot projector. For a smooth common final fold vector (a premise resting on the other seams), the restored outputs are related by the cusp deck \(T\) of C38d in the core and by \(+e_4\) (south) or \(-e_4\) (north) in the collar. Ingredients: \(K(u)+K(1/u)=1\) and linear seeds give \(H_f(\hat\psi p)=H_f(p)+1\) on \(1/4<Y_0<2\) for every \(q<1/8\); \(Y_0\in(1/4,2)\) in core and collar (collar \(q\in(1/1375,1/8)\) proved; measured \(Y_0\in[0.580,0.996]\)); \((x_1,y_1)=\hat\psi(x_0,y_0)\) for every moment \(f\); in the collar \(\rho_1=\rho_0+P\) and \((-iP,0)=\pm(e_4-c\cdot e)\mp(\theta/2\pi)e_1\), the ten-vector changes by \(e^{-i\theta_c}\), the corrected \(\mathrm{one}\) is unchanged mod 1, and \(D\) is invariant: ALL PASS. Exact witnesses at \(|t_c|=1/5,\,2/5\) on both sides of the selector, both signs. C43f, numerical (d), sample size 8: projector difference \(\le8.6\cdot10^{-122}\), deck residual \(\le1.7\cdot10^{-120}\), measured \(\Delta\mathrm{four}=+1\) south and \(-1\) north. |
| C44a–e, collar/pure-middle joins (a,b,c,d) | Along the whole seams \(|t|=1/2\) (south, selector `abs(t)<.5`) and \(|\mathrm{small}|=1/2\) (north, selector `inf`), for every argument chart of \(a\) and every fixed branch of \(\log t\). At \(\lambda=1\) every moment branch gives \(f=2\sqrt q\sinh(kb)\), so \(\rho=P(a+1/2)\), the source's \(\lambda=1\) return. With the C43 lattice identity, the \(\lambda=1\) collar data \(M_1\) have \(x_4=\xi_4+a+1/2\) (south) or \(-(a+1/2)\) (north); the actual collar data differ by flat terms, using §11.2's smoothness of the moment in \(\lambda\) up to \(1\). These equal the pure-middle \(x_4\) modulo an integer, while \(x_1\) differs by a real gauge. Hence the projectors agree, the corrected \(x_1\) agree modulo \(1\), and, for a smooth final fold vector, the restored outputs differ by a holomorphic lattice translation. \(1-\lambda\) and the cusp weight are flat at the seam: ALL PASS. Exact witnesses lie at \(|t|,|\mathrm{small}|\in\{499/1000,1/2,501/1000\}\). C44e, numerical (d), sample size 36: across both seams, for nine \((\Arg,a)\) choices including the half-centred chart and two with \(\Arg t\in(-3\pi/4,-\pi/4)\) (where the projector sees \(x_4\)), the \(\Psi\) jump is an integer lattice vector plus a residual of order \(\varepsilon\) (\(10^{-12}\to10^{-24}\)), as is the projector jump. |
| C45a–c, height moment and cells (a,b,c) | The executed bodies of `inverse_height_moment`, `_height_R` and `_height_terminal_forward` are bound to the source. Every alternative of a stable rewrite is the same function: \(x=e^2\operatorname{expm1}(2(a-\delta)/(c-a))=\varphi(a)\), the softplus forms, \(c\ell/(1+\ell)=\psi^{-1}(S)\), the positive root of \(x-1/x=L\), and \(\delta+c\ell/(2(2+\ell))=\varphi^{-1}(x)\). Both terminal branches equal \(z_c-\eta/\log(E_0+M)\), since \(1+F(e^U-1)=e^{E_z-E_0}(1+D)\). At \(|b|=2/5\) the interior composite equals the outer formula once \(T=\mathrm{id}\), and \(T-\mathrm{id}=\eta\ell/(L(L+\ell))\) is flat. At \(|b|=1/5\), \(R(a)-a=c\,l/(L(L+l))\) with \(l\) a smooth function of \(w=e^{-L}\) vanishing at \(0\). The \(|b|=1/10\) join is §11.2's written blow-up argument. The \(R^{-1}\) switch at \(|b|=b^*(q,\lambda)\in(0.15,0.2)\), where \(h(R(b^*))=c\), is flat because \(R-\mathrm{id}\) is; and \(h(2/5)>c\) since \(G_*(c)<0.446<0.4485<z_c\). Height-cell seams: core \(A=e^{2a-1}\) at \(a=1/2\) and \(Y=e^{-2a-1}\) at \(a=-1/2\) on both cells; collar \(\rho\) equal on both cells for every \(\kappa\): ALL PASS. |
| C46a–e, pole-cap internals (a,c,d) | Whole cap, both signs. Full bodies of `elementary_pole_patch`, `_elementary_pole_radial` (with its signature default) and `elliptic_pole_raw_inverse` are bound, and all 32 conditionals are classified. The identities are checked by executing the bound source statements under sympy, recording both alternatives of each symbolic if-expression. Both signs share the height \(\tan(2\pi h)=2\pi p_2/(1-2\pi p_0)\), and the executed height is \(\pm\operatorname{Arg}\) of the dominant coordinate; the \(AA\) and \(yy\) alternatives are identical (also for a symbolic \(q>0\)); \(|p_S|=s\); the executed chain gives \(|t_c|^2=K\,uu/\mathrm{den}\) with the positive root in den; `source_pc` is the Hopf vector of the slice pair \((\nu\bar x,y)\); and the returned pair has \(|\mathrm{out}|^2=s\). On the whole raw domain \(0<u\le10^{-4}\): \(|h|\le1.0008\cdot10^{-4}<1/10\), the height argument stays in the right half-plane, and \(|tc|^2<4.9\cdot10^{-7}\); the projection length (from (11.11)) satisfies \(|w|^2/u\ge0.9999\). The radial maps are \(S=\psi^{-1}T\varphi\) and \(S^{-1}=\varphi^{-1}T^{-1}\psi\), with identical rewrite alternatives; both outer seams have exact flat-factor forms \(t^2\ell/(c+t\ell)\), \(\ell=\log(1+QD(Q))\), \(Q=e^{-c/t}\); the inner-seam inputs of §11.3 are checked. Cap margins: \(r_L<1.6667\cdot10^{-6}\) (bound \(1.7\cdot10^{-6}\)); \(|z_S|^2<1.7003\cdot10^{-6}\) (bound \(2\cdot10^{-6}\)); four-seed \(A,Y<2\cdot10^{-6}\): ALL PASS. C46e, numerical (d), 32 production points: slice versus patch \(8.8\cdot10^{-105}\); flat defects at the seams scale as predicted. |
| C47a–f, slab phase band and its selectors (a,b,c,d) | The full bodies of `_elementary_corrected_section`, `elementary_phase_family`, `elementary_projection_section`, `_elementary_unit` and `_elementary_pole_inclusion` are bound. All 52 conditionals are classified, and all 8 raises sit in guards. Dataflow checks: `cap_band` feeds only the work estimate, and the diagnostics never reach an output. The band cutoff at \(1-|X|=\alpha\) is a flat join: the weight is \(s(1-u)\), with explicit derivative constants for \(k\le4\) and a written induction for all orders. The mismatch is unchanged under \(xx\mapsto\lambda xx\), since both sections are \(S^1\)-equivariant. The two stereographic chart values of \(\Theta\) agree exactly on the connected band \(|n_2|<3/5\): every identity of the argument is checked, including the radius bounds and the equator correction. The constant \(\operatorname{Arg}g(\mathrm{north},\mathrm{anchor})\) is one fixed point per sign: ALL PASS. C47f, numerical (d), reduced counts (a different map): the charts agree on the overlap, and the band selector behaves as stated. |
| C48, domain guards (a,b,c) | The functions executed by `elementary_global_coordinates` are listed (26) and the list is checked against the static call graph: every other reachable function of `cusp.py` (seven `elliptic_*` functions of the legacy path) is called only in an arm that the elementary path never enters. The 25 `raise` statements of the 26 functions are read from the source and each is matched with the bound that excludes it, so the list is complete (25/25). The bounds: the slab satisfies \(1-|X|\ge\alpha/2\); the fold overlaps exceed \(1/2\) (\(2L/N<1/2\) with the exact integers of (11.14)); the phase increments have \(\operatorname{Re}>7/8\) (\(4L_g/K<1/2\)); the hemisphere chart radius is \(\le1\); \(|a|<3/5\), and \(b<-21/50\) in core cell 1 (\(\cos(4\pi/25)<9/10\)); the core \(q<1/8\) (L13d) and the collar \(q>1/1375\) (C43b); the four-seed inputs satisfy \(\max(A,Y)<4\) in every cell; the cap radii are \(<2\cdot10^{-6}\); and on pure-middle data \(b_0=b_i=0\) exactly (C42a), so the height-free torus pair is never formed: ALL PASS. The dead-arm argument assumes `elementary` is True in every frame (C50h); implicit divisions, roots and `abs` are covered by C50f. |
| C49a–c, invertible differential everywhere (a,b) | On the slab \(G_W\) is \(S^1\)-equivariant (fold step, start vector, band correction and final phase), so \(\det DG_W=\det D\Lambda\). The orbit map \(\Lambda\) factors into: the raw base \((n,X)\mapsto R(X)O^{-1}n\) with \(R'/R=2/(\sqrt{2-X^2}(1-X^2))>0\); the conformal inverse stereographic map; \((u_1,u_2)\mapsto(u_1/u_2,\Arg u_2)\) and its north counterpart \((u_2/u_1,\Arg u_1)\); the profile \(t=Rg(|R|^2)\), with Jacobian \(g(g+2rg')>0\); the height maps with derivatives \(-f_a(1\mp f/R)/2\), \(-f_a/(2R)\) and, in core cell 1, \(2f_a/(R(R-f))\), nonzero since \(f_a>0\) by (11.13c) and the coordinate \(A\) or \(Y\) is chosen by the sign of \(f\) (C50e); and \(\partial x_4/\partial\rho=\pm1/P\) or \(\partial x_4/\partial\alpha=1/(2\pi)\). The restoration is invertible iff \(\det(e_2,e_3)\ne0\) (cusp) or the lattice determinant is \(\ne0\) (middle), independently of \(D\). On the cap, the radial map increases, the patch has \(\det=SS'/r>0\) (N10) and the raw \(g_j\) has \(\sigma_{\min}>0.99\) by (11.11): ALL PASS. The identities are type (a); the cited inputs (11.13c), (11.11) and §11.6 are types (b)/(c), N10 is a script check (types (a)/(b)), and C38–C47, C50 are type (a) with (d) observations. With §11.6's bijectivity, \(G_{\mathrm{cl}}\) is a diffeomorphism. |
| C50, C50a–h, the whole \(G\) path classified (a,b,c) | C50a pins the whole text of `cusp.py` (SHA-256 of its token stream; short digests localise changes in the 26 path functions) and checks that the runtime function objects are those the text defines; the certified object is that text read in exact arithmetic (§11.9): `elementary_global_coordinates(Z,X,phi,max_work=None)` with the shortcuts off and `force_slab=False`. For it: no `while`, `try`, `match`, `lambda`, comprehension filter, `%` or `//`; every comparison is a test or defines one of eight flags used only as tests or flag arguments; every call is a path function, local helper, listed builtin or listed mpmath function. C50: the 111 distinct tests (142 occurrences) are classified: 29 static, 24 guards (C48, C42a), 9 working-precision shortcuts (§11.7), 14 identical rewrites (C39a, C41a–b, C45a, C46b–c, C47c), 17 flat joins or weights (C41b–c, C44c, C45b, C46c, C47b, L13), 4 chart choices (C41, C44), 14 seams (N9, C38, C41–C47). C50b: the 29 `arg`/`log`/`atan` calls (implicit cuts) have positive or constant arguments, or a cut certified by C34 or C42. C50c: the 22 `floor`/`nint`/`sign`/`int`/`min`/`max` calls are precision settings, diagnostics, a guard, height cells (C45c), zero-weight torus-pair shifts (C41b), the principal-log normalization (C34, C42), a sign used only for \(|a|>1/10\), or unreachable. C50e: in the core (\(\lambda=0\)) the moment has the sign of \(b\) (\(f=b\), \(|f|\ge1/10\), \(|f|\ge e^{-1/5}-e^{1/5}/8\) on the three branches), the cell is 0 or 1, and the roots \(\sqrt{A/e}\), \(\sqrt Y\) are positive off the node, so the central fibres are not seams. C50f: all 280 divisions, square roots, `log1p`, `abs` and non-integer powers have their argument in the domain where they are smooth: 111 occurrences are recognised automatically (constant, inside a test, dead arm, \(|z|^2\)), and the other 169 occurrences (139 distinct) are listed with a reason. C50g: across the collar logarithm cut and the final-phase Arg cut the restored outputs differ by \(k\,e_4(L)+\tfrac{k(k-1)}2e_1\), \(0\) or \(e_1\). C50h: `elementary` is True in every executed frame. C50d: every check label cited in the four C50 tables exists: ALL PASS. |
| C51a–d, south-band phase shortcut (a,c,d) | C51a checks the Pauli projector norm factor, C51b the cancellation in each phase increment, C51c the monotone south radius and C51d the rational endpoint constants: ALL PASS. The type-(c) bound (11.15a) gives \(|\Theta_0|\le2300/2812499471<10^{-6}\) on the whole south band, independently of the fixed fold and phase counts. C47f's two reduced-count south samples are only type-(d) observations; a north-band shortcut remains open. |
| C52a–d, raw meridional derivative (a,b,d) | C52a checks the exact derivative of the stereographic radius, C52b the two manifest squares and C52c the conformal inverse stereographic identity, giving \(\|Du\|<6\) uniformly on the slab: ALL PASS. C52d is a separate type-(d) sample at 20001 values of \(X\), each with its worst angular direction; its maximum is below 3.814212, leaving a measured margin above 2.185788. The bound improves only the first factor of the deliberately coarse projector certificate (11.13). |
| C53a–d, inverse ratio/profile stage (a,b,c,d) | C53a checks the two exact radial derivative formulas. C53b checks the positive coefficients of both polynomial inequalities, yielding radial inverse norms below 12 and 40 on the source-selected domains; C53c checks the series lower bound, domain endpoints and angular comparison identities: ALL PASS. With \(\|Du\|<6\) from C52 and the dominant-coordinate ratio bound, \(\|D\mathrm{small}\|<480\), the reciprocal non-infinity branch has \(\|Dt\|<1920\), and the cusp coordinates satisfy \(\|Db\|<961\). C53d samples 20001 radii per branch (type (d)): minimum radial derivatives 0.2145580812 and 0.0318077198, measured margins above 0.1312 and 0.0068. This does not bound the remaining height and embedding derivatives. |
| m1–m3, formula-branch source audit and observations (d) | The checker identifies and gates all six working-precision identity returns, and checks nested/exception restoration: ALL PASS. Two strict interior type-(d) samples have exact-formula value gaps \(9.0129315017\cdot10^{-99}\) (pole radial) and \(3.23193414504\cdot10^{-92}\) (height radial) from the default identity returns, evaluated at 240 digits. Sample size two; m3 reproduces the observations. This does not evaluate the fixed-count global tensor. |
| N9a–i, cap/slab overlap (a,c) | for both signs, exact phase telescoping, source-linked raw-reference/projector identification, matching cusp inputs, torus restoration and first/second jets: ALL PASS, conditional on the written §11.3, §11.5 and §10.6 bounds; the post-raw radial margin is checked in N9b |
| N10a–d, radial shell (a,b) | source-linked finite elementary \(S,S^{-1}\) nodes, exact scalar inverse and derivative identities, 16 first and 64 second entries for each radial direction, and \(\det DD=S(r)S'(r)/r>0\) on \(10^{-8}<r<1.5\cdot10^{-8}\): ALL PASS; the raw pole and torus stages are not composed here |
| N11a–e, dominant raw-reference profile and selector seam (a,b) | the ordered source path selects \(T=u_2/u_1\) when \(|u_1|>|u_2|\), and \(T=u_1/u_2\) otherwise; its positive-root \(\mathrm{small}=T\sqrt{c/[A(1-z)+\sqrt{A^2(1-z)^2+4Bz}]}\), \(z=|T|^2\), \(c=2B\) or \(2\), has explicit value, 4 first and 8 second real-ratio entries on each branch: ALL PASS. The strict \(X=0\) witnesses have \(z=4/45\) and \(9/16\). \(A=(4\pi)^2>0\), \(B=e^4>0\), \(0\le z\le1\) make the root denominator positive. At \(|u_1|=|u_2|\), \(D_iD_d=4B\) proves the effective \(t\) and its first two jets agree; the source's inf test is false there and its phase input is \(u_2\) on both sides. The height and ten-slot projector are separate. |
| N12a–c, one active inverse-first height stage (a,b) | For the exact south collar witness \(a=b=3/20,\ |t|=2/5,\ \nu=1,\ \lambda=1/2\), the ordered source path through the shifted radial map and terminal-forward correction has a 30-node finite circuit. It computes the terminal output \(v\), its three first and six independent second partials in \((b,q,\lambda)\), with no unevaluated derivatives: ALL PASS. The rational bounds \(1/2500<q=(4/25)e^{-16/21}<1/8\), \(2<-\log q<4\), and \(0<\eta/(z_{\mathrm{cut}}-z)<9/40\) certify the selected open source branches and exclude the precision shortcuts. The separate one-point numerical check N12d gives value error \(<10^{-75}\), lower/upper \(q\) margins \(0.0742842,0.0503158\), and \(\eta/\mathrm{gap}\) upper margin \(0.2012902\). This is one height stage; the later inverse maps, \(Y,\rho\), and the six-real tensor are not included. |
| N13a–j, selected outer collar seams (a,b,c,d) | In a strict lower angular cell around \(\Arg t=-2\pi/3\), at \(|t|=1/2\) south and \(|t|=2\) north, the source positive-root profiles, Cramer correction and \(\lambda=1\) height give matching projector and phase-corrected middle \(W\)-slice two-jets modulo deck: ALL PASS. The raw reference tuples need not agree. The exact witnesses use \(a=1/20\) and \(9/20\); the Cramer determinants exceed \(5/8\) and \(40\), and both \(q\) values lie in \((1/2500,1/8)\). Separately labelled one-point numerical margins to those determinant bounds are \(0.6498254818\) and \(10.39720771\); lower/upper \(q\) margins are \(0.0654992845/0.0591007155\) south and \(0.0008069875/0.1237930125\) north. Other seam cells and the global tensor remain open. |
| N14a–d, completion of the N12 active height cell (a,c,d) | At \(b=3/20,\ |t|=2/5,\ \nu=1,\ \lambda=1/2\), the ordered source \(G^{-1}\), \(R^{-1}\), \(G\), positive quadratic \(Y\), and \(\rho=\log Y/2\) extend N12 to a 40-node exact circuit with three first and six independent second partials in \((b,q,\lambda)\): ALL PASS. Outward rational intervals prove \(v\in(19/100,21/100)\), \(G^{-1}(v)\in(1/10,19/100)\), \(R^{-1}G^{-1}(v)\in(1/10,1/5)\), \(f\in(19/100,21/100)\), \(Y\in(1/10,1/4)\), and \(\rho\in(-1,-1/2)\), selecting every source branch and positive denominator. The separately labelled one-point numerical comparison has \(|\rho_{\rm DAG}-\rho_{\rm source}|=3.1109\cdot10^{-91}\). Its measured lower/upper \(G^{-1}\) margins are \(0.06316123157/0.02683876843\), the \(R^{-1}<1/5\) margin is \(0.03670151219\), and the lower/upper \(Y\) margins are \(0.09136272386/0.05863727614\). This is one strict height cell, not the full free-slab tensor. |
| North raw-pole annulus, f1–f6 (a,b,c,d) | on \(3\cdot10^{-8}<|z_L|^2<10^{-7}\), \(1+2\beta>0\), \(0<\mathrm{height}<1/10\), \(t_c\ne0\): ordered north source guards, 28-node scalar 2-jet circuit, actual round and torus stages, 64+216 compositional Hessian skews, and named \(J=M^{-1}J_0M\): ALL PASS using the written (11.11) and §11.6 bounds. Proven raw singular-value margin \(>1999/2000\); one separately labelled numerical witness gives \(\|Dg-I\|_{\mathrm{op}}=3.9588551\cdot10^{-7}\) versus its bound \(3\cdot10^{-4}\), \(\sigma_{\min}=1.0000000086\), and source-value error \(3.53\cdot10^{-64}\). |
| South nonidentity radial transition, g1–g6 (a,b,c,d) | on \(1.3\cdot10^{-8}<|z_L|^2<1.4\cdot10^{-8}\), \(p_3(\operatorname{Hopf}z_L)>0\), \(\Im(2z_0\overline z_1)>0\): the source-bound \(h(Z,\phi)=\tau_f^{-1}(D^{-1}(g_0(D(z_L(Z)))),\phi)\) has a six-real value, Jacobian and six Hessian circuits; three 64-entry four-dimensional and one 216-entry six-dimensional compositional Hessian skews pass. The named \(J=M^{-1}J_0M\) satisfies the checked 36/216 square/Nijenhuis identity using (c8) and the stated invertibility bounds: ALL PASS. Both radial calls take nonidentity branches; the proved raw singular-value margin is \(>99993/100000\). At one rational witness, the separate numerical value error is \(4.3267\cdot10^{-84}\), measured raw defect \(1.1986845\cdot10^{-7}\) versus its pointwise bound \(6.6998742\cdot10^{-5}\), and sampled \(\sigma_{\min}=0.9999999867\). |
| North nonidentity radial transition, h1–h5 (a,b,c,d) | On the same strict radial and Hopf chamber, the source-bound north map \(h(Z,\phi)=\tau_f^{-1}(D^{-1}(g_\infty(D(z_L(Z)))),\phi)\) has separate scalar \(S,S^{-1}\) two-jet DAGs, a six-real value, Jacobian, and six Hessian circuits. Three 64-entry four-dimensional and one 216-entry six-dimensional compositional Hessian-skew checks pass; the named \(J=M^{-1}J_0M\) has the checked 36/216 square/Nijenhuis reduction using (c8), N10, (11.11), and §11.6: ALL PASS. Both radial calls take nonidentity branches; the proved raw singular-value margin is \(>99993/100000\). At one rational witness, the separate numerical patch-value error is \(3.6645\cdot10^{-84}\), measured raw defect \(8.7626736\cdot10^{-8}\) versus its pointwise bound \(6.6998742\cdot10^{-5}\), and sampled \(\sigma_{\min}=1.0000000019\). |
| `symbolic_identities.py`, pure middle (a,b) | actual six-coordinate \(\tau_f^{-1}\) on \(1/2\le |t|\le2\): explicit Jacobian, positive determinant \(2\pi(1+|t|^2)((|t|^2-1)\log|t|+2|t|^2)\), \(J^2=-I\), and all 216 \(N_J\) entries: ALL PASS; \((|t|^2-1)\log|t|\ge0\) is manifest |
| `symbolic_identities.py`, one active south collar chamber (a,b,c,d) | exact Cramer formulas, finite four-seed height and torus correction on \(3/10<|t|<1/2\), \(1/16<A,Y<1/2\); compositional 216-entry Hessian-skew check and positive determinant: ALL PASS; one 60-digit production comparison has error \(3.93\cdot10^{-61}\) and is numerical only |
| `symbolic_identities.py`, south raw-pole annulus (a,b,c,d) | on \(3\cdot10^{-8}<|z_L|^2<10^{-7}\), \(\Re u_2>0\), \(0<\mathrm{height}<1/10\), \(t_c\ne0\): source-linked raw inverse, round cap, and finite torus restoration form an explicit 43-node scalar second-jet circuit and a named six-real \(J=M^{-1}J_0M\); its 36 square and 216 Nijenhuis entries vanish by the checked compositional skew identity, using the written (11.11) and §11.6 bounds: ALL PASS. The proven raw singular-value margin is \(>1999/2000\); at one rational witness the separately labelled numerical values are \(\|Dg-I\|_{\mathrm{op}}=5.4871\cdot10^{-7}\) versus its bound \(3\cdot10^{-4}\), and \(\sigma_{\min}=0.9999999270\). |
| symbolic_identities.py (i1–i6), one source-linked free-slab fold (a,b,c) | On the strict C31 pure-middle projector cell, a six-variable Taylor ring carries a ten-slot rank-one source projector through one actual normalized projection and one right-half-plane phase increment, with a value, six first and 21 symmetric second entries per scalar. Generic product/unary coefficients, arbitrary local gauge cancellation, source AST call sites, exact certified \(N,K\), and an asymmetric-Hessian rejection pass: ALL PASS. The induction over fixed counts is conditional on valid jets and nonzero/right-half-plane bounds at every still-unchecked projector cell and seam; it is not the complete free-slab \(J\) circuit. |
| symbolic_identities.py (j1–j5), pure-middle negative-real log cut (a,b,d) | On the strict wedge \(1/2<|t|^2<3/4,\ \Re t<0,\ |\Im t|<-\Re t\), the source raw representative has an exact meridian witness and the canonical Cramer coordinate changes by \(\xi_4^+=\xi_4^-+1\) under \(L^+=L^-+2\pi i\). The sole active projector is \(E_{99}\) on both sides, with equal two-jets. Given a common nonzero endpoint vector, the corrected middle \(W\) data agree after the explicit \(e_1+e_4\) deck change: ALL PASS. The Cramer determinant is \(>3/2\). At \(t=-3/4\), its exact excess above \(225/128\) is \(175\log(4/3)/256>0\); the separately labelled one-point numerical value is \(0.196657666715085\). The final \(W\) comparison depends on the finite-fold denominator and endpoint premises; other paths and seams remain open. |
| symbolic_identities.py (k1–k4), exact finite scalar phase fold (a,b) | For unit source references \(r_j\), the production \(N+1\) rank-one updates give \(E=r_N(c_0/|c_0|)\prod_{j=1}^N(c_j/|c_j|)\), where \(c_0=\langle r_0,I_0\xi\rangle\) and \(c_j=\langle r_j,r_{j-1}\rangle\). Source AST, one-step algebra, induction at the certified exact \(N\), all 28 coefficients of a six-variable scalar-phase second-jet update, variable-reference gauge telescoping, and local argument/deck changes pass: ALL PASS. The actual C31 pure-middle same-cell reference gives one manifest positive adjacent overlap. This formula is conditional on nonzero source overlaps and valid reference jets along the whole path; it does not yet supply the full free-slab \(h\) or \(J\). |

| symbolic_identities.py (l1–l7), one all-mesh south free-slab tensor (a,b,c,d) | On the strict open chamber \(10^{-4}<|Z|^2<1/7000,\ |n_0|<1/100,\ n_1<-3/4,\ 1/2<n_2<2/3,\ X=-\sqrt{1-|Z|^2}\), the actual source meridian remains in one positive-\(x\) cusp cell for every \(0\le j\le N=4\cdot10^{8000}+1\). Forty scalar nodes carry six gradients and 21 symmetric Hessians each. The exact normalized-overlap product, source endpoint gauge, four-seed height, and torus restoration define the six-real \(h,M=Dh,D^2h\) and \(J=M^{-1}J_0M\); its 36 square and 216 Nijenhuis reductions pass: ALL PASS. C37d proves \(\det M\ne0\) on this subchamber. The independent local factor audit gives \(1-|t|^4>15/16\). The \(N=1\) six-output source comparison is a separately labelled one-point numerical transcription check, not a certified-count execution. Other free-slab cells and global gluing remain open. |
| N8, cap/source comparison (d) | three high-precision inner-cap values per sign: maximum coordinate errors \(1.1\cdot10^{-72}\) south and \(1.1\cdot10^{-75}\) north; these samples guard transcription only |
| New terminal-cut height, C23d | 132 value/deck/collar cases and 108 seam jets through order three: error 6.4·10⁻⁸⁶; sampled f_a∈[.466903,6.682791] versus proved type-(c) bound (11.13c) |
| New two-way height and forward evaluator, numerical | 276 cell/seam/axis/node compositions, error 1.0·10⁻⁵⁶; 38 reduced-count cap/free, exact half-cell seam, deck/lattice and independent-source compositions, error 1.2·10⁻⁵⁴ |
| C28e, numerical | 120 two-cap patch/inverse compositions, relative error <2·10⁻⁶⁸; sampled raw-collar residual zero at 85-digit working precision; phase-equivariance error <2·10⁻⁸⁶; sampled source squared radius <1.701·10⁻⁶ versus proved <2·10⁻⁶ |
| C28f, numerical | 10 two-cap inverse-Jacobian compositions, entrywise residual <2·10⁻²⁴; all sampled determinants positive |
| C28g, numerical | 120 production versus independent pole-patch/forward comparisons, relative error <2·10⁻⁶⁸ |
| Earlier elementary height, historical | 132 value/deck/collar cases and 108 seam jets through order three; error 6.4·10⁻⁸⁶; sampled f_a∈[1,4.536928] versus its earlier proved [8/9,14] |
| Open central-axis charts | 12 comparisons with the closed logarithmic inverse and 12 cell-seam jets, error 2.6·10⁻⁸⁶; 8 exact axis representatives also pass embedding and torus restoration checks |
| Inverse-first caps, earlier blend | 64 independent forward compositions, relative error 2.1·10⁻⁶³; sampled value/derivative defect constants 2.552/11.348 versus proved 200/5000; sampled earlier-blend defect 3.2·10⁻⁷ versus proved 13200ε=1.32·10⁻⁴<1/1000 |
| Concrete torch embedding | 4096 tests for each earlier height, error ≤4.7·10⁻¹⁴; per height/cusp 1024 chart comparisons, 256 exact-axis tangent comparisons (≤8.6·10⁻¹⁰) and 128 exact homogeneous cap models |
| Elementary scalar embedding/torus | 70 cusp/collar cases, including branch, deck and circle changes and an independent torch torus check; embedding error 7.4·10⁻⁵⁶ |
| Earlier raw global base/projector, historical | 96 round trips, error 1.4·10⁻⁴⁷; projector residual 2.0·10⁻⁵⁶; 16 derivative samples give maximum Frobenius norm 32.8697 versus the earlier proved operator bound 10¹¹⁴ |
| New global projector, numerical | 16 derivative samples give maximum Frobenius norm 32.8697 versus the proved operator bound 10⁸⁰⁰⁰; the analytic certificate is (11.13) |
| Finite real logarithm with interval parameter | 12 tests with multiple turns of phase and both real chart values, error 2.5·10⁻⁵⁴ |
| Full pole/slab overlap, **reduced counts N=512,K=32** | errors 3.1·10⁻⁶⁰ and 1.2·10⁻⁵⁴; sampled projection norm ≥.05504, phase increment real part ≥.91584 |
| Elementary pole shells | 8 production forward/inverse compositions, relative error <2·10⁻⁷⁶; 2 cap tensors inside the new shell, J²=−I and minimum sampled \|det DG\|=.00633258 |
| Elementary pole seam guards | 18 two-cap forward/inverse compositions including inner gaps 10⁻⁶, 10⁻¹⁰, 10⁻⁴⁰ and outer gaps 10⁻⁶, 10⁻³⁰; maximum relative round-trip error 3.0·10⁻⁵⁹; precision shortcuts certified in §11.7 |
| New local tensors | six pole/round-axis points and a free point with N=512; nonzero DG, J²=−I, independent sphere-chart change and integer torus changes pass; free-point \|det DG\|=2217.99 |
| Central-circle tensors | both a=1/2 circles, including transverse derivatives in the open axis charts; \|det DG\|=2.67101 and .0900846; J²=−I passes |
| Independent new-tensor finite-difference N | two caps and one free N=512 point, h=2·10⁻⁷ and 10⁻⁷; normalized residuals fall by a factor of four, maxima 4.73·10⁻¹⁰ and 1.18·10⁻¹⁰; nonintegrable controls .566/.0255 |

The N14 interval guard uses exact fractions rounded outwards to a grid of
\(10^{-30}\). For \(0\le x\le1\), the exponential series through degree 50
has remainder between zero and \(x^{51}/(51!(1-x/52))\). For
\(1\le m\le2\), put \(z=(m-1)/(m+1)\in[0,1/3]\); the logarithm series
through \(2z^{101}/101\) has remainder between zero and
\(2z^{103}/(103(1-z^2))\). Integer square-root floors enclose each
positive square root. Exact range reduction and these bounds propagate
through all 40 scalar nodes; the displayed rational brackets are proved
by interval inclusion, not by the one-point comparison. In particular
\(G^{-1}(v)<19/100<1/5\) and
\((1/5)/(1/5-G^{-1}(v))<20\), below the source precision-shortcut
threshold \(>(1+10)\log 10>22\); the certified inverse-\(R\)
intermediates satisfy \(y>1\) and \(L>0\).

For C36, write \(y=(u,v,0)\), \(S=1+u^2+v^2\),
\(Q=(1-u^2-v^2)^2+4u^2\), \(m=-n_2>0\), and
\(\Lambda=(1-X^2)/(1-X\sqrt{2-X^2})\). Let
\(s=|t_c|^2\), \(\Delta=\sqrt{f^2+4q}\), and \(f_a>0\) be
the selected source height derivative. The positive profile derivatives
are
\(z_s'=(A+2Bs+ABs^2)/(1+As)^2\) and
\(z_n'=(As^2+2Bs+AB)/(As+B)^2\), where
\(A=(4\pi)^2,\ B=e^4\). In the stated coordinate orders, the two
source-linked determinant factorizations are
\[
 \det Dh_s=
 -\frac{\Lambda^2\Lambda_X S(1-s^2)f_a}
 {4mv^4z_s'\Delta},\qquad
 \det Dh_n=
 -\frac{4\Lambda^2\Lambda_X S(1-s^2)f_a}
 {mQ^2s^4z_n'\Delta}.
\]
Here \(\Lambda_X/\Lambda=
2/(\sqrt{2-X^2}(1-X^2))>0\).
The source witnesses have \(u=1/100\) and \(v=12/5\) or \(4\);
their strict branch inequalities and all other denominator signs
are checked in C35–C36. The fibre phase has derivative \(+1\)
under the source finite projection recurrence, and the torus
restoration has identity \((\phi_2,\phi_3)\) block. These give the
remaining factors in the displayed full determinants.

For C37 let \(n_2=+\sqrt{1-n_0^2-n_1^2}\),
\(\Lambda=(1-X^2)/(1-X\sqrt{2-X^2})\), and
\(y=(\Lambda n_0,1-\Lambda n_2,-\Lambda n_1)\).
Put \(S=1+|y|^2\), \(v^2=y_1^2+y_2^2\),
\(s=|t|^2\), \(q=e^{-4s/(1-s)}s\),
\(\Delta=\sqrt{a^2+4q}\), and
\(z_s'=(A+2Bs+ABs^2)/(1+As)^2\), with
\(A=(4\pi)^2,\ B=e^4\).
In the period-one Hopf coordinate \(\theta\) and target order
\((\Re z_1,\Im z_1,\Re z_2,\Im z_2,\Re z_3,\Im z_3)\), the exact
source-linked factorization is
\[
 \det Dh=
 \frac{\Lambda^2\Lambda_X S(1-s^2)}
 {4n_2v^4z_s'\Delta}>0,\qquad
 \frac{\Lambda_X}{\Lambda}
 =\frac{2}{\sqrt{2-X^2}(1-X^2)}>0.
\]
The source inequalities give \(\Lambda/10<a<\Lambda/5\) and
\(|t|<\Lambda/15\), whence \(a>3\sqrt q/2\).
Every normalized reference has first-slot probability \(>4/5\);
the initial pole inclusion has first-slot probability \(>7/8\).
The reverse triangle inequality gives \(>3/5\) for the initial
overlap modulus and for
\(\Re\langle r_j,r_{j-1}\rangle\) between canonical references.
Consequently every source projection-overlap modulus is \(>3/5\).
These are uniform symbolic bounds at the certified count. The numerical
sample in the table checks transcription only. The production
evaluator's default work guard stops a literal run at this count;
the checker proves the source-defined finite expression by induction.

Section (l) of the symbolic checker composes this source path all the way to
the specific local tensor. Write \(r=|Z|^2\), \(X=-\sqrt{1-r}\),
\(\xi=Z/\sqrt r\), and \(n=\operatorname{Hopf}(\xi)\). On its stated open
subchamber, \(1+X>r/2>1/20000>10^{-5}\), and
\(\Lambda(X)<r/(1+|X|)<(10/19)/7000<10^{-4}\), so every C37 branch
and overlap bound applies. With
\(h_j=-1+10^{-9}+(X+1-10^{-9})j/N\), let \(r_j\) be the normalized
ten-slot source reference at \(h_j\). Its only nonzero entries are
\((x_j,t_j/x_j)/\sqrt{x_j^2+|t_j/x_j|^2}\) in slots 0 and 1;
the forty scalar assignments in (l2) give \(x_j,t_j\) by positive
roots from the displayed raw \(y\), ratio and cusp height formulas.
Define
\[
 \gamma_0=
 \frac{\langle r_0,I_0\xi\rangle}{|\langle r_0,I_0\xi\rangle|},
 \qquad
 \gamma_j=
 \frac{\langle r_j,r_{j-1}\rangle}
      {|\langle r_j,r_{j-1}\rangle|}\quad(1\le j\le N),
 \qquad \alpha=\prod_{j=0}^{N}\gamma_j .
\]
C37 proves the denominators are uniformly \(>3/5\), including at
the production integer \(N=4\cdot10^{8000}+1\). The checker assigns
all six first and 21 independent second derivatives to every source
stage and overlap. In particular, for real input indices \(i,k\),
\[
 \alpha_i=\alpha\sum_{j=0}^{N}\frac{\gamma_{j,i}}{\gamma_j},
 \qquad
 \alpha_{ik}=\alpha\left[
  \left(\sum_{j=0}^{N}\frac{\gamma_{j,i}}{\gamma_j}\right)
  \left(\sum_{j=0}^{N}\frac{\gamma_{j,k}}{\gamma_j}\right)
  +\sum_{j=0}^{N}\left(
   \frac{\gamma_{j,ik}}{\gamma_j}
   -\frac{\gamma_{j,i}\gamma_{j,k}}{\gamma_j^2}\right)\right].
\]
At the endpoint put \(x_W=\overline\alpha x_N\),
\(y_W=\alpha t_N/x_N\), \(t=x_Wy_W=t_N\),
\(q=e^{-4|t|^2/(1-|t|^2)}|t|^2\), and
\(H_f=(e^{-4|t|^2/(1-|t|^2)}|x_W|^2-|y_W|^2)/(1-q)\).
For \(C=\phi+2H_f(\Re t,\Im t)/(1-|t|^2)\) and
\(B=t(C_2-iC_3)\), the exact source output on this chamber is
\[
 h=(x_We^{-B},\,y_We^B,\,
       C_2(1-t^2)+iC_3(1+t^2)),\qquad
 J=(Dh)^{-1}J_0Dh .
\]
The script constructs \(Dh\) and all six Hessians in the input order
\((\Re Z_0,\Im Z_0,\Re Z_1,\Im Z_1,\phi_2,\phi_3)\).
The positive C37 determinant makes the inverse legitimate. Its
actual symmetric Hessians, the checked 36 square reductions and
216 transported Nijenhuis reductions give a direct executable
certificate on this open cell. The one-point \(N=1\) source
comparison checks transcription only; the exact production-count
claim is the finite product and its induction, not a literal run
of the astronomical loop. No global overlap claim follows from
this local computation.

C38 treats the raw reference's argument-selector join on an open
overlap. Let \(p=u_2\) at cusp \(0\) and \(p=-\overline{u_1}\) at cusp
\(\infty\). In the lower wedge
\(-\pi<\arg p<-5\pi/6\), the principal value
\(a_0=\arg p/(2\pi)\) and the half-centred value
\(a_1=1/2+\arg(-p)/(2\pi)=a_0+1\) both lie in valid source
charts. With \(b=a_0\), \(q=e\nu^2|t_c|^2\),
\(e=\exp(-4|t_c|^2/(1-|t_c|^2))\), and
\(E=\exp(-2b-1)\), the selected outer-core formula has
\(f=-E+q/E\) and \(\sqrt{f^2+4q}=E+q/E\). Thus the negative
and positive references are related exactly by
\[
 (x_0,y_0)=(t_c/\sqrt E,\sqrt E),\qquad
 (x_1,y_1)=\bigl(1/(m y_0),\,m x_0y_0^2\bigr),\qquad
 m=\nu\sqrt e .
\]
The first embedding weight is zero and the second is one on both
charts; the negative chart applies this same deck map before its
second frame. Hence the normalized ten-slot vectors and their
rank-one projectors agree as functions on the open overlap, with
all first and second jets. In the upper wedge
\(5\pi/6<\arg p<\pi\), the two argument values are already equal.
The exact source-round witnesses and strict perturbations in C38
meet both sides of the selector at both cusps.

On the lower overlap the source four-seed formula satisfies
\(H_f(x_1,y_1)=H_f(x_0,y_0)+1\). With
\(\sigma=1\) at cusp \(0\), \(\sigma=-1\) at cusp \(\infty\),
the restored holomorphic coordinates have the transition
\[
 (z_1,z_2,z_3)\longmapsto
 \left(\frac1{\nu z_2},\ \nu z_1z_2^2,\
 z_3+2\sigma z_1z_2\right),\qquad \det_{\mathbb C}DT=1 .
\]
Here \(z_2\ne0\) on the overlap. This equality of open chart
formulas transports all their jets; it does not assert that the
remaining projector path cells or source-chart seams have been
covered.

C39 treats the dominant-coordinate selector in the raw reference
on an open overlap. Set \(R=u_1/u_2\), \(r=|R|^2\),
\(A=(4\pi)^2\), \(B=e^4\), \(b=A(1-r)\), and
\(D=\sqrt{b^2+4Br}>|b|\). The direct source branch has
\(t=R\sqrt{2/(b+D)}\). On the inverse-ratio branch the source
first forms \(1/R\) and then takes its reciprocal; its effective
coordinate is
\[
 t=R\sqrt{\frac{D-b}{2Br}}.
\]
The identity \((D-b)(D+b)=4Br\), together with positivity of both
square roots, makes these the same smooth function for \(r>0\).
C39a checks the source assignments and the radial derivatives
through order two (identity type (a)).

For the strict band
\[
 \frac{99}{100}<r<\frac{101}{100},\qquad
 |\Arg R|<\frac{\pi}{5},\qquad |\Arg u_2|<\frac{\pi}{5},
\]
write \(s=|t|^2\). The inverse profile is
\(r=(Bs^2+As)/(1+As)\), with derivative
\((A+B(2s+As^2))/(1+As)^2>0\).
The explicit bounds \(A<160\), \(50<B<81\) show
\(r(9/100)<99/100\): after multiplying by its positive
denominator, the available gap is
\(99/100-(81\cdot81+160\cdot9)/10000=1899/10000\).
Likewise \(r(1/4)>101/100\), with gap
\(50/16-101/100-160/400=343/200\).
Thus \(3/10<|t|<1/2\). The source infinity flag is false
on both sides, so both use \(t_c=t\), the same argument \(a\)
and the same active south collar. Here
\(q=s\exp[-4s/(1-s)]\) satisfies
\[
 \frac1{100}<\frac{9}{100e^2}<q<
 \frac1{4e}<\frac18.
\]
The lower bound uses \(4s/(1-s)<4/3<2\); the upper
uses \((se^{-4s})'=e^{-4s}(1-4s)>0\) for \(s<1/4\).
The source Cramer determinant is \(>981/5000\), the
south torus determinant \(1-s^2>15/16\), and the
unnormalized reference norm square is
\(\beta^2+(1-\beta)^2=1/2+2(\beta-1/2)^2\ge1/2\).
These are explicit-constant inequalities and manifest squares
(types (b)–(c)); C39b and C39d check their algebra and source
guards.

At the exact round witness
\(\lambda=2-\sqrt2\),
\(X=(\lambda-1)/\sqrt{1+\lambda^2}\),
\(Z=(\sqrt{1-X^2},0)\), the source has
\(u_1=u_2=1/\sqrt2\) and \(t=e^{-1}\).
The values \(\lambda\pm1/10000\) give strict witnesses on
opposite selector sides inside the displayed band (C39c).
On this band the source normalized ten-vector has only slots
\(0,1,8\), with entries
\[
 \frac{\beta}{\sqrt{\beta^2+(1-\beta)^2}}
 \frac{(\overline{t/e^\rho},\,e^\rho)}
      {\sqrt{|t/e^\rho|^2+e^{2\rho}}},
 \qquad
 \frac{1-\beta}{\sqrt{\beta^2+(1-\beta)^2}}
 e^{2\pi i\,\mathrm{one}}
\]
in the first two and last slot, respectively. Since the
source data agree as functions on the open band, so do this
reference, its rank-one projector and all their second jets.
The same Cramer, height, phase and torus formulas then give
identical restored-output second jets (C39d). This equality
propagates through each finite meridian step wherever its
projection denominator is nonzero. C39e reports a separate
one-point numerical source comparison and the measured margins
listed in §11.8; it does not supply a global bound.

C40 treats the infinity/middle source selector
\(|\mathrm{small}|=1/2\) in one north angular/height
chamber. This is a **flat boundary join**, rather than
equality of the two raw formulas on an open overlap. Set
\(s=|\mathrm{small}|\), \(t=1/\mathrm{small}\),
\(r=|t|=1/s\), and \(\phi=\Arg t\). At the exact
seam witness choose \(\Arg u_1=\pi/4\), \(\phi=\pi/3\);
the infinity-side argument is \(a=3/8\).
C40b constructs actual round witnesses for
\(s=499/1000,1/2,501/1000\) from
\[
 p=s^2,\qquad
 \kappa^2=\frac{p(p+A)}{Ap+B},\qquad
 u_1=\frac{e^{i\pi/4}}{\sqrt{1+\kappa^2}},\quad
 u_2=\frac{\kappa e^{-i\pi/12}}{\sqrt{1+\kappa^2}},
 \qquad A=(4\pi)^2,\ B=e^4 .
\]
The inverse stereographic formulas recover this spinor from
the resulting round \((Z,X)\), and the source positive root
recovers \(\mathrm{small}=s e^{-i\pi/3}\).
The two strict side witnesses select different branches.
Explicit bounds give \(9/25<\kappa^2<49/100\),
\(-1/3<X<0\), \(n_2>59/100\), and
\(1/2025<q<13/2500<1/8\).
For the bound \(B=e^4<55\), the exponential series
through degree four is \(65/24\), while its tail from
degree five is
\(<(1/120)\sum_{j\ge0}6^{-j}=1/100\).
Hence \(e<65/24+1/100<68/25\) and
\((68/25)^4<55\).
The Cramer determinant exceeds \(30\), and the north torus
determinant exceeds \(9/10\) for
\(499/1000<s<501/1000\) at these fixed angular inputs.
All inequalities are strict, so this witness has a small
open angular/height neighborhood (explicit-constant
inequalities, type (c)).

On the infinity side let
\(P=\log(e^{-2}s)-2s^2/(1-s^2)<0\) and
\(\eta=\rho/P\). The source Cramer coordinates are
\[
 (\mathrm{one},\mathrm{two},\mathrm{three},\mathrm{four})
 =\left(\frac{\eta\phi}{2\pi},
 -\frac{2\eta\Re t}{r^2-1},
 -\frac{2\eta\Im t}{r^2-1},-\eta\right).
\]
At \(s=1/2\), the inverse-first height has
\(\lambda=1\), \(k=-\log q=-2P\), and
\(f(a)=2\sqrt q\sinh(ka)\).
Its positive root is \(Y=\sqrt q e^{-ka}\), giving
\(\rho=P(a+1/2)\) and \(\eta=a+1/2=7/8\).
The terminal-cut parameters satisfy
\(z_0=0,\ \Lambda=\mu=1,\ C_0=C_1=0\), so the
interior terminal formula is the identity at this parameter
value. Since \(1/5<a=3/8<2/5\), both radial maps
in that interior composition are on identity branches.
Moreover, \(z_{\rm cut}-\Delta z
=2\sqrt q[\sinh(2k/5)-\sinh(3k/8)]
>\sqrt q\,k/20>1/180\), using
\(q>1/2025\) and \(k>5\).
Consequently its source \(E_z\) satisfies
\(E_z<e^{9/5}<9<50\), selecting the low-\(E_z\)
terminal formula on a neighborhood of the seam.
C40a and C40c–c3 check these source signs and identities
(identity type (a)).

For \(s<1/2\), the north embedding weight
\(\beta=1-\operatorname{flatstep}((s-3/10)/(1/5))\)
and \(\Delta=\rho-P(a+1/2)\) have zero value,
first and second derivatives at \(s=1/2\);
on the other side \(\beta=0\).
C40c checks the flat-step derivative limits.
In the angular sector \(\pi/4<\phi<\pi/2\),
both surviving angular phase brackets reduce identically
to zero. The infinity-side reference is therefore,
through order two, \(e^{i\eta\phi}\) times the
right-side angular reference. Its cusp block and all
cusp/angular cross-projector entries have zero two-jet;
the angular projector block agrees. Here \(s>499/1000>2/5\)
gives \(\beta<1/2\) by monotonicity of the flat step,
so the unnormalized norm square is
\(\beta^2U+(1-\beta)^2(1+\omega^2)>1/4\), where
\(U\ge0\) is the cusp-block norm square and
\(0<\omega<1\) is the second angular weight.
C40d checks all
transverse, tangential and mixed second-jet entries using
a generic flat factor \(x^3b\) (identity type (a)).

If the finite-fold endpoint jets match, the final inner
product is nonzero, and one local argument branch is fixed,
the source phase correction cancels
\(e^{i\eta\phi}\). The corrected middle
\(\mathrm{one}\) coordinates then agree through order two
modulo an integer deck shift. At the seam,
\(q=e^{-16/3}/4\), and the four-seed torus inputs are
\(A_f=q^{1/8}\), \(Y_f=q^{7/8}\).
The strict \(1/2500<q<1/256\) implies
\(1/16<A_f<1/2,\ Y_f<1/16\), and
\(qA_f,qY_f<1/16\). Hence their source branches
and the torus correction are smooth locally; the restored
middle output has the same second jet (C40e).
Propagation through the full finite meridian fold is
conditional on all intermediate projectors having matching
two-jets and nonzero projection denominators.
C40 does not certify those premises along every global path.

C41 treats the cusp/middle source selector \(|t_c|=3/10\) at
both cusp signs. For \(|t_c|\le3/10\) the raw reference returns
the core reference of `elementary_cusp_reference`. For
\(3/10<|t_c|<1/2\) it returns middle data
\((t,\mathrm{one},\mathrm{four},\log t,a)\), built from
`elementary_cusp_logradius` and the real Cramer solve, with
\(t=t_c\) south and \(t=1/t_c\) north. At a given point both
sides use the same argument value \(a\), the same
\(q=e\nu^2|t_c|^2\) with \(e=\exp(-4|t_c|^2/(1-|t_c|^2))\), the
same cell \(n=\lfloor a+1/2\rfloor\) and the same moment
\(f=f(a-n,q,\lambda)\). Let \(M_0\) denote the right-hand
formulas with \(\lambda=s((|t_c|-3/10)/(1/5))\), the weight
deficit \(1-b_0\) (south) or \(1-b_\infty\) (north), and \(b_m\)
all replaced by zero. Unlike C39–C40, C41 is a join along the
whole seam, not in one chamber. It holds for each argument chart
of \(a\) and each fixed branch of \(\log t\); C34 relates the
branches of \(\log t\).

*Reference.* Put \(R=\sqrt{f^2+4q}\). The middle height is
\(e^{2\rho}=q^nY\) with \(Y=2q/(R+f)=(R-f)/2\), and C34e
reconstructs \(zz=e^\rho\), so the \(M_0\) reference is
\((t_c/e^\rho,e^\rho)\). In each of the cells \(n=-1,0,1\) the
core reference is
\[
 (x,y)=\bigl(\overline g\,t_c/e^\rho,\ g\,e^\rho\bigr),\qquad
 g=t_c/|t_c|\ (a>0),\quad g=1\ (a\le0).
\]
For \(a>0\) this reduces to \(|x|^2e^{2\rho}=|t_c|^2\), and for
\(a\le0\) to \(|y|^2=e^{2\rho}\). Both follow from
\((R-f)(R+f)=4q\). A unit scalar multiplies both source frames,
plain and deck-shifted, by \(g\). On \(|H|<3/5\) the core weights
\((\mathrm{hw}(H),\mathrm{hw}(H\mp1/2))\) and the shift integers
of `torus_pair` select the same frame wherever the weight is
nonzero. Hence the normalized ten-vectors differ by the factor
\(g\), and on the left of the seam the rank-one projectors
coincide identically (C41b, identity type (a)).

*Flatness on the right.* Near the seam put \(u=5|t_c|-3/2\).
South, \(\lambda=1-b_0=b_m=s(u)\); north, \(\lambda=1-b_\infty=b_m=s(u)\).
The other cusp weight is exactly zero, and the restoration
weight \(\beta\) is the same \(1-s(u)\). The flat step has value
and first two derivatives tending to zero as \(u\downarrow0\).
Every right-hand source quantity minus its \(M_0\) value either
carries such a factor, or is \(F(\lambda,\cdot)-F(0,\cdot)\) for
an \(F\) smooth in \(\lambda\ge0\). The latter holds because the
\(\lambda=0\) terminal data equal the general data at
\(\lambda=0\), the \(\gamma=0\) quotients are removable, and
§11.2 proves smoothness of \(f\) in \(\lambda\) on its second
parameter piece. These differences therefore have zero
complete two-jets on the seam (C41a, C41c; types (a) and (b)).

*Phase and restoration.* Let \(E\) be the common final unit
vector of the finite fold. On the left it lies in the image of
\(P_c=P_{M_0}\), so \(p=\langle\mathrm{ref},E\rangle\) has modulus
one. The source inner product is conjugate-linear in its first
slot, so \(p_m=g\,p_c\). (A linear first slot would give
\(\bar g\,p_c\), and the \(W\)-points would differ by \(g^2\);
C41d checks this control.) The restore loop rebuilds \(y\) from
the middle \(W\)-point by the same reconstruction as C34e. With
\(\mathrm{one}\) shifted by \(\Arg p_m/(2\pi)\), its exponent is
\(\rho+i\Arg p_m\) at both signs. So both final \(W\)-points have
\(y=g\,p_c\,e^\rho\) and \(x=t_c/y\), and both restorations call the
four-seed \(H_f\) at the same point. On the right, \(E\) lies in the
image of the actual projector, which differs from \(P_{M_0}\) by flat
terms, so there the identities below hold up to flat terms. The paper's holomorphic chart changes out of
the middle chart are
\[
 \Phi:\ y=e^{i(z-\pi s(s-1))},\ x=t/y,\ w_c=w-2ts,\ \
 s=\frac{\log t}{2\pi i};\qquad
 \Theta:\ y=e^{i(z-iw/t)},\ x=t^{-1}/y,\ w_c=w/t^2 .
\]
Their complex Jacobian determinants are \(i\) and \(-i/t^4\), and
`tfull.c0_to_mid` and `ci_to_mid` are their inverses. Let \(c\)
be the real coordinates of the cusp vector \(e_4\), and \(cc\)
those of \(-w_c\). The source middle correction then satisfies
\(D-cc=H_fc\) at both signs. South, this follows from the
lattice \(w\)-equation \(2ts=\xi_1e_2+\xi_2e_3+2t\xi_3\); north,
from \(cc=\mathrm{four}\cdot c\). Hence, for every \(\phi\), the
c0/ci restoration of the core \(W\)-point equals \(\Phi\)
(resp. \(\Theta\)) of the \(M_0\) middle restoration exactly
(C41d, identity type (a)).

Consequently, at every point of the seam the source projector
and the restored holomorphic output have matching complete
two-jets across \(|t_c|=3/10\), through a holomorphic chart
change. Higher jets would also match, given §11.2's smoothness
of \(f\) in \(\lambda\); C41 checks second order. Inside the finite fold the same
identity applies to every meridian sample that meets this
seam, so this selector creates no jet mismatch at any \(P_j\).
Propagation through the fold remains conditional on the nonzero
projection denominators and on the other source seams. The
argument selector of C38 is certified only for \(|t_c|<1/10\);
its continuation to \(1/10\le|t_c|<1/2\), which meets this seam
at a corner, is not covered by C41; C43 covers it.

*Domain and witnesses.* On the band \(2/25<|t_c|^2<1/10\) both
sources are in their certified domains (type (c)). There
\(1/25<q<1/10\) south and \(1/1375<q<1/10\) north, so
\(1/2500<q<1/8\). The Cramer determinant exceeds
\(2r^2(1+r^2)>0\). The torus determinants are \(1-|t_c|^4>0\) and
\(|t_c|^4-1<0\). Finally \(|a|<3/5\) in both argument charts,
because \(\cos(\pi/5)<9/10\). In fact production gives only
\(-0.43<a<0.58\), so cell \(-1\) never occurs there; C41b covers it
anyway. C41a binds the whole real Cramer solve
(`scalar_lattice_coords`, including its determinant helper). The exact round witnesses use
\(p=|t_c|^2\in\{(299/1000)^2,\,9/100,\,(301/1000)^2\}\),
\(A=(4\pi)^2\), \(B=e^4\) and \(Q=\sqrt{1+k^2}\). South:
\(k^2=(Bp^2+Ap)/(1+Ap)\) and \((u_1,u_2)=(-k/Q,i/Q)\), giving
\(t_c=i\sqrt p\). North: \(k^2=(p^2+Ap)/(Ap+B)\) and
\((u_1,u_2)=(i/Q,k/Q)\), giving \(t_c=-i\sqrt p\). Both have
\(a=1/4\) and \(\Arg t=\pi/2\), away from the log cut, and
\(|X|<3/5\) (C41e). C41f is a numerical observation (d) only. It compares the
production branches at 18 points, covering the cells with \(a\le0\)
(including the deck frame), \(a=0\), \(a>0\) in cells \(0\) and
\(1\), and both argument charts.

C42 treats every source seam of the pure middle
\(1/2\le|t|\le2\) at once. These are the selector \(|t|=1\) between
the two formulas for \(\mathrm{four}=x_4\), and the
principal-argument cuts of \(u_1\), \(u_2\) and \(t\), including
their corners. Here the raw reference returns
\((t,0,\mathrm{four},\log t)\) with the principal logarithm, and
\[
 \mathrm{four}=\frac{\Arg u_2}{2\pi}+\xi_4+\frac12\ \ (|t|\le1),\qquad
 \mathrm{four}=\frac{\Arg u_1}{2\pi}-\lambda\Bigl(\frac{\Arg t}{2\pi}-\xi_4\Bigr)
 -\frac{1-\lambda}{2}-\frac12\ \ (|t|>1),
\]
where \(\xi_4\) is the fourth Cramer coordinate of the canonical
point \((\pi s(s-1),2ts)\), \(s=\log t/(2\pi i)\), and
\(\lambda=1-s((|t|-6/5)/(2/5))\). Middle data
\(W=(t,\mathrm{one},\mathrm{four},L)\) give the point
\(\Psi(W)=(2\pi\,\mathrm{one}-iL\,\mathrm{four},\,2t\,\mathrm{four})\)
upstairs. Let \(G\) be generated by
(a) \(\mathrm{one}\mapsto\mathrm{one}+1\),
(b) \(\mathrm{four}\mapsto\mathrm{four}+1\), and
(c) \((\mathrm{one},L)\mapsto(\mathrm{one}-\mathrm{four},L+2\pi i)\).
These change \(\Psi\) by the lattice vectors
\(e_1=(2\pi,0)\), \(e_4(L)=(-iL,2t)\) and \(0\), respectively.

*The pure middle is the inverse branch.* In the direct ratio branch
\(|\mathrm{small}|^2=2rr/(b+\mathrm{root})\le\sqrt{rr/B}\le e^{-2}\),
so \(|t|\ge1/2\) forces \(|u_1|>|u_2|\) and
\(t=u_1/(c\,u_2)\) with \(c>0\). Hence
\(\arg t\equiv\arg u_1-\arg u_2\pmod{2\pi}\). On
\(1/2\le|t|\le2\) the cusp weights \(b_0,b_\infty\) and both
restoration weights vanish exactly (C42a).

*Invariance.* The embedding first normalizes \(L\) to the principal
logarithm, which is the move (c). Its slot 8 is
\(d_0e^{2\pi i\,\mathrm{one}}\), and its slot 9 is
\(d_1e^{i(2\pi\,\mathrm{one}-\beta_9\mathrm{four})}\) with
\(\beta_9=\Arg(-t)+\pi-\Arg t\in\{0,2\pi\}\). So the ten-vector
\(V\) is invariant under (a) and (b), and
\(V(\mathrm{one},\cdot)=e^{2\pi i\,\mathrm{one}}V(0,\cdot)\). With a
continuous logarithm \(L_c=P+2\pi ik\), where \(P\) is principal,
the vector \(V_c=V(k\,\mathrm{four}_c,\mathrm{four}_c)\) is
continuous across the \(t\)-cut. There \(k\) increases by one while
\(\beta_9\) jumps from \(0\) to \(2\pi\), and \(d_0=0\) (C42b).

*Continuous four.* Fix continuous branches of \(\arg u_1\),
\(\arg u_2\), \(\arg t\) and \(L_c\) on an open set. The principal
values differ from them by \(2\pi n_1\), \(2\pi n_2\), \(-2\pi k\)
and \(-2\pi ik\). By C34c, \(\xi_4(L+2\pi ik)=\xi_4(L)+k\). Hence
the left formula changes by \(n_2-k\) and the right formula by
\(n_1\), for every \(\lambda\), since \(\Arg t/(2\pi)-\xi_4\) is
unchanged. Near \(|t|=1\), where \(\lambda=1\), the two continuous
formulas differ by the integer \(-N-1\), where
\(\arg_ct=\arg_cu_1-\arg_cu_2+2\pi N\) (C42c).

*Conclusion.* Write \(\mathrm{four}=\mathrm{four}_c+n\). Then
\(V=e^{-2\pi ik\,\mathrm{four}_c}V_c\), so the source projector is
that of the continuous data. For a common final unit vector \(E\),
the source phase is \(e^{2\pi ik\,\mathrm{four}_c}\langle V_c,E\rangle\),
so \(\mathrm{one}=\mathrm{one}_c+k\,\mathrm{four}_c+m\), and
\[
 \Psi(W)-\Psi(W_c)=(m-kn)\,e_1+n\,e_4(L_c).
\]
On the pure middle \(C=\phi\), and the restoration commutes with
lattice translations. The translation
\((z,w,t)\mapsto(z+2\pi a-ib\log t,\,w+2bt,\,t)\) is a holomorphic
chart change with complex determinant \(1\) (C42d, identity type (a)).
Thus, on every open subset of \(1/2<|t|<2\), the projector is
smooth. If the final fold vector \(E\) is smooth there (this rests on
the projectors along the whole meridian, hence on the other seams
and on nonzero overlaps), the restored output equals a smooth one
modulo holomorphic deck translations. This holds across \(|t|=1\), across
all three argument cuts, and at their corners. Section (j) of
`symbolic_identities.py` is the special case of the \(t\)-cut on
\(1/2<|t|^2<3/4\). Exact round witnesses of \(|t|=1\) and of its
corner with the \(t\)-cut use
\(u_1=i/Q\), \(u_2=k\,e^{i(\pi/2-\theta)}/Q\), \(Q=\sqrt{1+k^2}\),
\(k^2=(p^2+Ap)/(Ap+B)\), \(p=|t|^{-2}\), giving
\(t=|t|e^{i\theta}\) and \(|X|<3/5\) (C42e). C42f is a numerical
observation (d) only. It crosses each of eight seams with production
points at distance \(\varepsilon\); the jump of \(\Psi\) is an integer
lattice vector up to a residual of order \(\varepsilon\), measured at
\(\varepsilon=10^{-12}\) and \(10^{-24}\). The joins of the pure middle
with the collars at \(|t|=1/2\) and \(|t|=2\) are separate (N13, C40).

C43 extends C38 to every \(|t_c|<1/2\), that is, to the rest of the
core and to the collar, at both cusp signs. On the lower overlap
\(-\pi<\Arg p<-5\pi/6\) the principal value \(a_0\in(-1/2,-5/12)\) and
the half-centred \(a_1=a_0+1\) are both valid source charts (C38a); on
the upper wedge they agree. Both charts use the same \(b=a_0\) and
\(\lambda\), with cells \(n=0\) and \(n=1\).

*Height.* The source seed \(K\) satisfies \(K(u)+K(1/u)=1\) for every
\(u>0\): on \([1/2,2]\) both terms are \(v/(1+v)\), and otherwise the
reciprocal branch returns \(1-K\) of the other. It is linear on
\(u\le1/16\). The deck \(\hat\psi(x,y)=(1/(my),mxy^2)\),
\(m=\nu\sqrt e\), sends \((A,Y)=(q/Y_0,Y_0)\) to \((1/Y_0,qY_0)\) with
the same \(t\) and \(q\). If \(1/4<Y_0<2\) and \(q<1/8\), both points lie
in the four-seed domain. The \(A\)-seed \(K(q/Y_0)\) of \(p\) cancels the
\(qA\)-seed of \(\hat\psi p\), the \(qY\)-seed \(K(qY_0)\) of \(p\) cancels
the \(Y\)-seed of \(\hat\psi p\), and the \(q^2\) seeds are linear. Hence
\[
 H_f(\hat\psi p)-H_f(p)=K(1/Y_0)+K(Y_0)=1
\]
exactly (C43a). This replaces C38c's seed-branch computation, which
needed \(q<1/100\). On both charts \(1/4<Y_0<2\): with \(|b|>2/5\) the
moment is on its outer branch
\(f=-(e^{\kappa(|b|-1/2)}-qe^{\kappa(1/2-|b|)})\), and
\(|f|\le Y_0\le|f|+\sqrt q\). The core (\(\kappa=2\)) gives
\(|f|>e^{-1/6}-e^{1/6}/8>1/4\). In the collar \(q\) is unimodal in \(s=|t_c|^2\in(9/100,1/4)\), so
\(q>1/25\) south and \(q>1/1375\) north, while \(q\le1/(4e)<1/8\). Hence
\(2\le\kappa\le k=-\log q<8\), which gives
\(|f|\ge e^{-k/12}(1-e^{-5k/6})>e^{-2/3}(1-e^{-5/3})>1/4\). In both,
\(Y_0<1+\sqrt{1/8}<2\) (C43b, type (c)). Measured margins (type (d)):
an independent scan found \(|f|\ge0.578\) and \(Y_0\in[0.580,0.996]\),
against the proven \(|f|>1/4\) and \(Y_0\in(1/4,2)\); C43f prints its own
values at its eight points.

*Core.* The \(a_0\) reference is \(y_0=\sqrt{(R-f)/2}\), \(x_0=t/y_0\),
and the cell-1 \(a_1\) reference is \(x_1=\sqrt{2/((R-f)e)}/\nu\). So
\((x_1,y_1)=\hat\psi(x_0,y_0)\) for every moment \(f\). Both heights
have zero first weight. The second weights agree, and the negative
chart applies \(\hat\psi\) before its frame, so the ten-vectors coincide
exactly. With a common final vector the phases coincide (they have
modulus one, since the fold vector lies in the image of the final
projector), the final
\(W\)-points are \(\hat\psi\)-related (C38e), and the restorations differ
by the holomorphic cusp deck \(T\) of C38d (C43c).

*Collar.* The logradius gives \(\rho_1=\rho_0+P\) with \(P=\log(m|t_c|)\).
With \(t=re^{i\theta}\) the middle variable, the raw point moves by
\((-iP,0)\), and
\[
 \text{south: } (-iP,0)=e_4-c\cdot e-\tfrac{\theta}{2\pi}e_1,\qquad
 \text{north: } (-iP,0)=-e_4+c\cdot e+\tfrac{\theta}{2\pi}e_1,
\]
where \(c\cdot e\) is a torus direction, dropped from the \(W\)-data. So
\(\mathrm{four}\) changes by \(\pm1\), and \(\mathrm{one}\) by
\(-\theta_c/(2\pi)\) with \(\theta_c=\Arg t_c\). The ten-vector is
multiplied by \(e^{-i\theta_c}\): in the cusp pair the height and the
shift integer change by \(+1\) and \(-1\), so
\(yy\mapsto yy\,|mt_c|/(mt_c)\); in the angular slots by
\(e^{2\pi i\Delta\mathrm{one}}\). Hence the projectors agree,
\(p\mapsto e^{i\theta_c}p\), and the corrected \(\mathrm{one}\) is
unchanged modulo \(1\). In the restore loop \(yy\) is multiplied by
\(mt_c\), so \(H_f\) rises by one. The correction \(D\) depends on
\(H-\mathrm{four}\) (south) or \(H+\mathrm{four}\) (north) and is
unchanged, so the restored outputs differ by \(+e_4\) (south) or
\(-e_4\) (north), a holomorphic lattice translation (C43d). These
conclusions about the restored output assume a smooth common final
fold vector, which rests on the other seams and on nonzero overlaps.

*Witnesses.* Exact round witnesses lie at \(|t_c|=1/5\) (core) and
\(2/5\) (collar), at both signs, with
\(p=\mathrm{pr}-i\sqrt{1-\mathrm{pr}^2}\),
\(\mathrm{pr}\in\{-9/10\mp1/1000,\,-9/10\}\), parametrized rationally
by \(\mathrm{pr}=-(1-\tau^2)/(1+\tau^2)\). The south collar witness lies
in the inverse ratio branch. All have \(1<R^2<4\), hence \(|X|<3/5\)
(C43e). C43f is a numerical observation (d) only. At eight production
points it evaluates the other chart with the production core reference
or logradius at the same point, and checks the projector and the deck
relation.

C44 treats the two joins of the collars with the pure middle along
their whole seams. These are \(|t|=1/2\) (south; the selector
`abs(t)<.5`) and \(|\mathrm{small}|=1/2\), i.e. \(|t|=2\) (north; the
selector `inf`). N13 and C40 had checked one chamber of each. Let
\(M_1\) be the collar formulas with \(\lambda=1\) and the cusp weight
(\(b_0\) south, \(b_\infty\) north) set to \(0\), so \(b_m=1\). From
the left, \(1-\lambda\) and the cusp weight equal \(1-s(u)\) with
\(u\uparrow1\). This is flat because \(s(1-u)=1-s(u)\) (C44c) and C41c, and §11.2
gives smoothness of the moment in \(\lambda\) up to \(1\). So the collar
source differs from \(M_1\) by flat terms.

At \(\lambda=1\) every moment branch gives \(f=2\sqrt q\sinh(kb)\), with
\(k=-\log q\): directly on the inner and \(\lambda=1\) branch; with
\(B=C=\sqrt q\), \(\kappa=k\) on the outer branch; and through the
identity terminal map on the interior branch (C40c). Then
\(\sqrt{f^2+4q}=2\sqrt q\cosh(kb)\), both source \(Y\) formulas give
\(Y=\sqrt q\,e^{-kb}\), and \(\rho=nP+\log(Y)/2=P(a+1/2)\). This is the
\(\lambda=1\) return of `elementary_cusp_logradius` (C44a). By the C43
lattice identity the \(M_1\) data are
\[
 \text{south: } x_4=\xi_4+a+\tfrac12,\ \ x_1=\xi_1-(a+\tfrac12)\tfrac{\theta}{2\pi};\qquad
 \text{north: } x_4=-(a+\tfrac12),\ \ x_1=(a+\tfrac12)\tfrac{\theta}{2\pi}.
\]
The pure middle has \(x_4=\Arg u_2/(2\pi)+\xi_4+1/2\) near \(|t|=1/2\),
and \(x_4=\Arg u_1/(2\pi)-1\) for \(|t|\ge8/5\), where \(\lambda=0\); in
both cases \(x_1=0\). Since \(a=\Arg u_2/(2\pi)+j\) south and
\(a=1/2-\Arg u_1/(2\pi)+j\) north, with \(j\) an integer, the two values
of \(x_4\) differ by \(j\) or \(-j\). The two values of \(x_1\) differ
by a smooth real gauge \(g\) (C44b). At the seam both ten-vectors
consist of the two angular slots only, and an integer change of \(x_4\)
is invisible there (C42b). So the \(M_1\) vector is \(e^{2\pi ig}\)
times the pure-middle vector, and the projectors agree. Moreover
\(p_{M_1}=e^{-2\pi ig}p_R\), so the corrected \(x_1\) agree modulo
\(1\). Both restorations have \(C=\phi\), so the restored outputs differ
by a lattice translation, which is holomorphic with complex
determinant \(1\) (C44c). This last conclusion again assumes a smooth
final fold vector. The witnesses use the C42e family at
\(|t|\in\{499/1000,1/2,501/1000\}\) and the C41e north family at
\(|\mathrm{small}|\in\{499/1000,1/2,501/1000\}\) (C44d). C44e is a
numerical observation (d) only.

C45 certifies that the height moment \(f(b,q,\lambda)\), and the cusp
height cells built from it, form one smooth function on the \(G\) path.
The executed bodies of `inverse_height_moment`, `_height_R` and
`_height_terminal_forward` are bound to the source. Each conditional is
of one of five kinds:
- a domain guard;
- a working-precision shortcut, bypassed in formula mode and certified
  in §11.7;
- a pair of identical rewrites;
- a flat join;
- a height-cell seam.

*Rewrites* (C45a, type (a)). With \(\delta=1/10\), \(c=1/5\),
\(E_h=e^{c/(c-a)}\), \(\varphi=E_h-e^2\), \(\psi=E_h-e\) and
\(S(x)=\log(1+e^{x-1/x})\): the source's
\(e^2\operatorname{expm1}(2(a-\delta)/(c-a))\) is \(\varphi(a)\); its two
softplus forms agree; \(c\ell/(1+\ell)=\psi^{-1}(S)\) for
\(\ell=\log(1+S/e)\); \(e\operatorname{expm1}(v/(c-v))=\psi(v)\); the two
forms of \(\log\operatorname{expm1}\) and of the positive root of
\(x-1/x=L\) agree; and \(\delta+c\ell/(2(2+\ell))=\varphi^{-1}(x)\). For the
terminal map, \(1+F(e^U-1)=e^{E_z-E_0}(1+D)\), with \(U=E_z-A_{\rm cut}\) and
\(D=C_0e^{E_0-E_z}+C_1e^{E_0+A_{\rm cut}-2E_z}\). Hence both of its branches
equal \(z_c-\eta/\log(E_0+M)\) of (11.3).

*Joins* (C45b; types (a)–(c)). For \(|b|\ge2/5\), with \(R=\mathrm{id}\) and
\(T=\mathrm{id}\), the interior composite equals the outer closed form.
The switch is flat because
\[
 T(z)-z=\frac{\eta\ell}{L(L+\ell)},\qquad
 \ell=\log\Bigl(1+\frac{\log(1+D)}{E_z}\Bigr),
\]
where \(D=O(e^{-E_z})\), \(E_z=e^L\) and \(L=\eta/(z_c-z)\). At \(|b|=c\),
\[
 R(a)-a=\frac{c\,l}{L(L+l)},\qquad L=\frac{c}{c-a},\qquad
 l=\log\bigl(w(e+S(1/w-e^2))\bigr),\quad w=e^{-L}.
\]
Here \(w(e+S(1/w-e^2))=1+(e-e^2)w-w^2/(1-e^2w)+(\text{a flat term})\), so
\(l\) is smooth in \(w\) at \(0^+\) and vanishes there. Every derivative of
\(R-a\) in \(a\) is therefore a polynomial in \(L\) times a flat function.
The join at \(|b|=1/10\) is §11.2's blow-up argument, cited here.
A further switch is internal to \(R^{-1}\), at \(|b|=b^*(q,\lambda)\), where
\(h(R(b^*))=c\). It is flat because \(R-\mathrm{id}\) is: if \(R=\mathrm{id}+\phi\)
with \(\phi\) flat at \(c\), then \(R^{-1}=\mathrm{id}-\phi\circ R^{-1}\).
Part (i) uses \(h(2/5)>c\), which follows from
\(G_*(c)\le(3/10)\sinh(8/5)/(8/5)<0.446<0.4485<z_c\) (§11.2).

*Height cells* (C45c, type (a)). In the core (\(\kappa=2\)), at \(a=1/2\)
both cells give \(A=e^{2a-1}\), and at \(a=-1/2\) both give
\(Y=e^{-2a-1}\). In the collar, for every \(\kappa\), the logradius
values \(\log Y(b)/2\) and \(P+\log Y(b-1)/2\) coincide. These are
equalities of analytic functions on \(|a\mp1/2|<1/10\), so all jets agree.

C46 certifies the pole-cap branch of the slice (\(1-|X|<\alpha/2\)) at
both signs as one smooth map. The bodies of `elementary_pole_patch`,
`_elementary_pole_radial` and `elliptic_pole_raw_inverse` are bound to the
source. Their 32 conditionals are each static, a guard, a precision
shortcut (gated by `elementary_formula_branches`), an identical
alternative, a flat join, a no-op, or unexecuted. The identities are
checked by executing the bound source statements under sympy: an
if-expression with a symbolic test records both alternatives, so the
identities hold for the production expressions themselves. This covers
the selected height and ratio, the \(t_c\) chain, `source_pc`, and the
returned pair; the default `inverse=False` of the radial map, on which
the call for \(D\) relies, is bound as well.

In \(p=\operatorname{Hopf}(z)\) both signs share the height
\(\tan(2\pi h)=2\pi p_2/(1-2\pi p_0)\). The real part of the dominant
coordinate has the sign of \(1-2\pi p_0\ge1-8u>0\), so the principal
argument never meets its cut. The two alternatives of \(AA\) and of \(yy\)
are the identical functions \((R\pm h)/2\), \(xx\cdot yy=\nu^2|t_c|^2\)
gives \(|p_S|=s\), and
\(|w|^2=r(1+n_S\cdot n_L)/2+(|p_S|^2-s^2)r/(4s^2)\). On the whole raw
domain \(0<u\le10^{-4}\) the explicit constants are:
\(|h|\le u/(1-8u)\le1.0008\cdot10^{-4}<1/10\); \(|\mathrm{ratio}|^2<1.283\cdot10^{-6}\);
\(|t_c|^2<4.9\cdot10^{-7}\); \(R,s>0\); and, citing (11.9)–(11.11),
\(|w|^2/u\ge0.9999\). So no guard fires (C46b–b4). The returned pair has
\(|\mathrm{out}|^2=s\) exactly (C46b4).

The production radial maps are exactly \(S=\psi^{-1}T\varphi\) and
\(S^{-1}=\varphi^{-1}T^{-1}\psi\), and their \(f<0\), \(y>1\), \(L<0\)
alternatives are identical. At \(r=c\) and \(v=c\), with \(Q=e^{-c/t}\),
\[
 S(c-t)-(c-t)=\frac{t^2\ell}{c+t\ell},\qquad \ell=\log\bigl(1+QD(Q)\bigr),
\]
with an explicit \(D_0\) for \(S\) and \(D_1=e^3+m+2Q/(K+1+Qm)\) for
\(S^{-1}\) (C46c–c2). The smoothness of \(D_0,D_1\) at \(Q=0\), hence the
flatness, is written calculus. The inner seam \(|z|^2=10^{-8}\) rests on
§11.3's blow-up argument; C46c3 checks its scalar inputs, namely the
simple zero of \(\varphi\) with \(\varphi'(a)=6e^3/a\), and \(\psi(0)=0\)
with \(\psi'(0)=e/c\). On the cap, \(|X|>1-\alpha/2\), the denominator
exceeds \(2-\alpha/2\), \(r_L<1.6667\cdot10^{-6}\) and
\(|z_S|^2<1.7003\cdot10^{-6}\), and every four-seed input is \(<2\cdot10^{-6}\),
on its linear branch (C46d–d2). C46e is a numerical observation (d).

C47 certifies the phase-correction band of the slab
(\(\alpha/2\le1-|X|<\alpha\)) at both signs. The bodies of
`_elementary_corrected_section` (including `mismatch`),
`elementary_phase_family`, `elementary_projection_section`,
`_elementary_unit` and `_elementary_pole_inclusion` are bound to the
source. Each of their 52 conditionals is static, a guard (with the bound
that keeps it inactive), diagnostic bookkeeping, or a smooth join. A
dataflow check shows that `cap_band` feeds only the work estimate.

*Band cutoff* (C47b). Outside the band the output is the uncorrected
section. The weight is \(1-s(u)=s(1-u)\) with \(u=1\) at \(1-|X|=\alpha\),
and \(e^{-iw\theta}-1=-iw\theta\int_0^1e^{-isw\theta}\,ds\), so the
correction has a flat factor. Flatness to all orders is a written
induction. Its pieces are checked, with explicit constants for
derivatives through order four and generic jet checks.

*Hemisphere reference* (C47c). The two source formulas for the Hopf
preimage `xx` differ on the overlap by a unit factor \(\lambda\). The
direct and transported sections are both \(S^1\)-equivariant, so the
mismatch \(\langle\lambda d,\lambda s\rangle=\langle d,s\rangle\) is
independent of the choice. A control shows that rotating only one of
them would change it.

*Two-chart logarithm* (C47d, types (a)–(c)). Each chart value of
\(\Theta\) is a sum of \(\arctan(\operatorname{Im}R/\operatorname{Re}R)\)
with \(\operatorname{Re}R\ge7/8\), so \(e^{i\Theta}=g\) in both charts by
telescoping. On the connected band \(|n_2|<3/5\), where both
stereographic radii are \(\le2\), the difference of the two values is a
continuous \(2\pi\mathbb Z\)-valued function, and it vanishes at
\(e=(1,0,0)\) because of the equator correction. The hemisphere charts
have radius \(\le1\), so the radius guard never fires.

*Start constant* (C47e). \(\operatorname{Arg}g\) is evaluated once, at
the fixed point \((\mathrm{north},\mathrm{anchor})\) of each band, so it
creates no seam. C47f is a numerical observation (d) with reduced counts,
which define a different finite map.

On the checked south raw-pole annulus, write \(r=|Z_0|^2+|Z_1|^2\) and
\(z_L=(Z_0-iZ_1,Z_0+iZ_1)/\sqrt{4\pi(1+\sqrt{1-r^2})}\). Let \(w=g_0(z_L)\)
be the explicit positive-root projector formula (11.10), and put
\(t=\overline{w_0}w_1\), \(E=\exp(-4|t|^2/(1-|t|^2))\),
\(H=(E|w_0|^2-|w_1|^2)/(1-E|w_0|^2|w_1|^2)\),
\(c=2(\Re t,\Im t)/(1-|t|^2)\), \(C=\phi+Hc\), and
\(A=t(C_2-iC_3)\). The specific six-real-coordinate map checked by the
script is the real and imaginary parts of
\((\overline{w_0}e^{-A},w_1e^A,C_2(1-t^2)+iC_3(1+t^2))\).
Its circuit explicitly defines \(M=Dh\) and \(J=M^{-1}J_0M\). The round
factor has positive determinant, (11.11) gives
\(\sigma_{\min}(Dg_0)>1999/2000\), and the restoration factor has determinant
\(-(1-|t|^4)\ne0\). These bounds make the displayed inverse legitimate
throughout this one annulus chamber.
The rational local witness is \(z_L=(1/5000,(1-i)/10000)\), with
\(|z_L|^2=6\cdot10^{-8}\); the positive round-cap radial factor gives it a
round input.

The corresponding north raw-pole annulus has a separate concrete circuit.
For \(z_L=(z_0,z_1)\), put
\(p=2z_0\overline z_1\), \(p_3=|z_0|^2-|z_1|^2\),
\((\alpha,\beta,\gamma)=\pi(p_3,-\Re p,\Im p)\), and
\(d=1+2\beta+2(\alpha^2+\beta^2+\gamma^2)\). The north stereographic
coordinates are
\[
u_1=(-1-2\beta+2i\alpha)/d,\qquad
u_2=2(\beta+\alpha^2+\beta^2+\gamma^2+i\gamma)/d.
\]
With \(\nu=e^{-2}\), \(q=|u_2/u_1|^2\), \(A=(4\pi)^2\),
\(B=e^4\), and \(D=A(1-q)+\sqrt{A^2(1-q)^2+4Bq}\), set
\(t_c=(u_2/u_1)\sqrt{2B/D}\) and
\(h=-\arg(-u_1)/(2\pi)\). The raw inverse is the explicit
four-real formula \(w=\sqrt{s/|V|^2}V\), where
\[
E=e^{-4|t_c|^2/(1-|t_c|^2)},\quad
R=\sqrt{h^2+4E\nu^2|t_c|^2},\quad
U=(R+h)/(2E),\quad Y=(R-h)/2,\quad s=U+Y,
\]
\[
V=(Uz_0+\nu\overline{t_c}z_1,\ \nu t_c z_0+Yz_1).
\]
The north round cap feeds this formula with
\(z_L=(Z_0+Z_1,Z_0-Z_1)/
\sqrt{2\pi(1+\sqrt{1-|Z|^4})}\).
For \(x=\overline w_0/\nu\), \(y=w_1\), \(t=xy\), the final
torus stage uses \(E_t=e^{-4|t|^2/(1-|t|^2)}\),
\(H_f=(E_t|w_0|^2-|w_1|^2)/
(1-E_t|w_0|^2|w_1|^2)\),
\(c=2(\Re t,-\Im t)/(1-|t|^2)\),
\(C=\phi+H_fc\), and
\[
h_{\infty}=(x e^{-t(C_2+iC_3)},\
y e^{t(C_2+iC_3)},\
C_2(t^2-1)+iC_3(1+t^2)).
\]
The code binds the first and second derivatives of this specific
composition, forms \(J=M^{-1}J_0M\), and checks its Hessian skews;
the 36 square and 216 Nijenhuis entries then vanish by the checked
compositional identity. The round factor has positive determinant,
(11.11) gives \(\sigma_{\min}(Dg_\infty)>1999/2000\), and the
restoration factor has determinant \((1-|t|^4)/\nu^2>0\).
The same rational witness \(z_L=(1/5000,(1-i)/10000)\) lies in this
north chamber. This is a local annulus check, not a check of the
radial transition or the free slab.

For the checked south radial-transition chamber, use the preceding
south six-real formula with
\[
 w=D^{-1}\!\bigl(g_0(D(z_L))\bigr),\qquad
 D(z)=\sqrt{S(|z|^2)/|z|^2}\,z,\qquad
 D^{-1}(w)=\sqrt{S^{-1}(|w|^2)/|w|^2}\,w.
\]
The exact finite scalar nodes for \(S,S^{-1}\), including their
positive-root and logarithm branches, are N10 and g2; g3–g5 bind
their first and second jets around the 43-node south raw inverse and
the same torus restoration. Both radial stages use their nonidentity
branches on the strict chamber in the table. The rational witness
\(z_L=(1/10000,(5-3i)/100000)\) has
\(|z_L|^2=1.34\cdot10^{-8}\). The square and Nijenhuis reductions
are compositional applications of (c8) to these actual bound jets.
This is one open part of the south shell.

The new second-jet checks are exact on their stated domains. The pole check
uses the identity piece \(H_j(z_L)=z_L\) for \(|z_L|^2\le10^{-8}\), and its
four-seed height simplifies there to \(H_f=(A-Y)/(1-AY)\). The pure-middle
check starts after the round-to-\(W\) map, in local coordinates
\((\Re t,\Im t,x_1,x_4,\phi_2,\phi_3)\), where both cusp collar weights vanish.
The checks distinguish a coordinate Hessian from an arbitrary frame: their
noncoordinate controls have \(J^2=-I\) but nonzero \(N_J\). C30 proves the
second-jet step identities for the exact finite projection and phase counts,
and checks the relevant source call sites; it does not symbolically execute all
scalar branches of (11.18). C31 and C32 cover two nonempty projector cells.
The active collar calculation begins after the round-to-\(W\) map and uses a compositional
second jet of its explicitly displayed torus correction; its numerical source
comparison is only a transcription check. The south raw-pole annulus check
starts with the actual round inputs, passes through the raw inverse and torus
restoration, and stores \(J=M^{-1}J_0M\) and all six Hessians as a finite
straight-line circuit; its 216-component conclusion uses the source-specific
jet symmetries and the C30/c8 algebraic identity, without expanding the full
radical expression. N9 identifies the raw cap projector with the source
reference and verifies the overlap algebra and jets, using the existing
analytic cap and phase bounds. Thus these scripts do **not** yet constitute a
direct executable verification of the specific \(J_{\mathrm{cl}}\) on all of
the round \(S^4\times T^2\). Section (l) now supplies the actual complete
free-slab \(J\) circuit on one strict chamber. The other free-slab
reference, height and selector cells, remaining collar and source-chart
seams, and their second-jet correspondence are still open. N10 supplies
the nonidentity radial shell's two explicit scalar inverse functions
and their radial first/second jets; (g) and (h) compose them with the
respective raw pole inverse and torus restoration on one strict
transition chamber per sign. The remaining shell cells are open.
The global \(J^2=-I\) and \(N_J=0\) conclusion above continues to use the
written diffeomorphism/gluing certificate and the generic C10 identity.

Run the new implementation checks with
`python3 extra/s4t2_explicit/check_inverse.py --closed` and
`python3 extra/s4t2_explicit/check_inverse.py --closed-nijenhuis`.
`elementary_tensor_numeric` is an independent high-precision derivative check of
(11.19), using `mp.diff` in interiors of source charts. It is not the definition
of the symbolic derivatives, and it does not imply a global guarantee for a
smaller projection count. In particular .05504 is **below** the full expression's
proved per-step bound 1/2: it belongs to a different, reduced-count numerical
instance. A preliminary N=64 instance actually had a zero projection. Neither
observation changes (11.14) or its proof; no smaller count is certified here.

The checks caught and fixed a genuine logarithm-branch metadata error in the
scalar embedding: the supplied height referred to a principal-log cusp, while
conversion used the unreduced logarithm. The explicit normalization in §11.1
fixes it, and the branch assertions remain unchanged. The central-fibre review
also replaced polar value formulas by the smooth open formulas (11.6a), including
neighbouring height cells. An initially ineffective N-control sheared inside a
single complex line; the checked control now mixes a sphere and a torus direction.
No legacy assertion or tolerance was weakened.

The complete legacy safety net was rerun: `check_tcore`, `check_tlift`,
`check_lift`, `check_full`, `check_inverse.py 3000`, the four numpy checks,
`check_degree.py 1121 128 128`, `nijenhuis_round.py 2000`, and
`check_jround.py 300` all pass. The degree is 1.000000 with positive determinant
at all 18,366,464 nodes; the legacy normalized N maximum is 2.59·10⁻¹⁰.
All four symbolic suites pass. Both manuscript Python checks pass; LaTeX is
**not run** because `latexmk` is unavailable, and the manuscript is unchanged.
These legacy tests protect the existing construction; they are not numerical
evaluation of the huge new expression.

The two-way finite elementary inverse-pair statement requires human review
under AGENTS.md. The proof is the blockwise identities and type-(c)
certificates above; the reduced-count forward and tensor checks are numerical
safety checks. The fixed products (11.14) have not been numerically executed,
and no smaller counts have been certified. Legacy I3–I5/T14 remain open for
their original fitted evaluator.

### 11.9 The direct global certificate of \(J_{\mathrm{cl}}\) (T16)

This subsection assembles the checks of §11.8 into a global statement. It
is a mathematical statement for human review. It takes a different route
from the one T16 specified, so it also needs the task owner's agreement:
instead of a branch-circuit checker that differentiates each actual
branch, it combines a generic integrability identity with seam-by-seam
smoothness and global invertibility certificates.

**The object.** \(G_{\mathrm{cl}}\) is (11.18): the function defined by
the text of `elementary_global_coordinates(Z,X,phi,max_work=None)` in
`cusp.py`, with `force_slab=False` and all working-precision shortcuts off
(`elementary_formula_branches()`, section (m) of `symbolic_identities.py`),
read in exact real arithmetic. Decimal literals are rationals
(`mp.mpf('.3')`\(=3/10\)), `mp.pi` and `mp.e` are \(\pi\) and \(e\), the `mp`
functions are the mathematical functions, and `workdps`, `max(dps)` and
unary `+` are identities. An mpmath run of the same text is a numerical
evaluation (d) of this function, not the function itself; C48 and C50
classify the exact reading (in it, for example, `lam==1` in the collar
and the `cell==-1` alternative are never taken). C50a pins the whole text
of `cusp.py` and checks that the runtime objects the scripts inspect are
the ones it defines. \(J_{\mathrm{cl}}=G_{\mathrm{cl}}^*J_Y\)
on the round \(S^4\times T^2\subset\mathbb R^5\times\mathbb R^2\times\mathbb R^2\).
In every local chart it is
\[
 J_{\mathrm{cl}}=(Dh)^{-1}J_0\,Dh,\qquad
 h=\phi_\beta\circ G_{\mathrm{cl}},
\]
where \(\phi_\beta\) is one of the paper's holomorphic charts (the cusp
charts c0, ci or the middle chart). On each open cell of its selectors,
each \(h\) is a finite composition of arithmetic, positive square roots,
\(\exp\), \(\operatorname{expm1}\), \(\log\), \(\operatorname{log1p}\),
\(\arctan\), \(\sinh\), \(\operatorname{asinh}\), the principal
\(\operatorname{Arg}\) off its cut, \(|\cdot|\) off zero, the flat step and
locally constant integers (`floor`, `nint`, `sign`), with the fixed integer
counts (11.14). (11.19) and its seven-tuple rules give every coefficient of
\(J_{\mathrm{cl}}\) without a differentiation operator.

**Hypotheses (written arguments, for human review).** The statement uses
the following written type-(b)/(c) certificates and calculus arguments,
which the scripts check only in their arithmetic or identities:
- the projector bound (11.13), the height bounds (11.13a)–(11.13c), and
  §11.2's terminal bounds and its blow-up argument at \(|b|=1/10\);
- the raw-lift bounds (11.8)–(11.11), §11.3's inner-seam blow-up (C46c3)
  and §11.5's cap estimate (G2's initial overlap, G14);
- §11.5's derivative propagation behind G3: \(\|D\xi\|\le10\),
  \(\|De_0\|\le4L\), \(|A_j-A_{j-1}|\le8L\), \(\|DE\|\le(5+8N)L\) and
  \(\|Dg_j\|\le(9+8N)L<L_g\);
- §10.6's two-chart logarithm (C47d) and §11.6's bijectivity;
- the all-order flatness of the flat joins: the written calculus of C46c2,
  the written induction of C47b, the smoothness of the moment in
  \(\lambda\) used by C41c and C44, and C45b;
- C43b's collar \(q\) bounds;
- the tensoriality step of item 3 below.

**Guard bounds.** C48 matches each `raise` on the path with one of:
G1 the slab margin \(1-|X|\ge\alpha/2>10^{-9}\); G2 fold overlaps
\(>1/2\) ((11.13)–(11.14)); G3 phase increments with
\(\operatorname{Re}>7/8\) ((11.14), (11.16)); G4 hemisphere chart radius
\(\le1\); G5 \(|X|<1\) on the slab; G6 only the node has \(r=0\), and it is
the pole; G7 \(|a|<3/5\) on the \(|t_c|\le3/10\) branch, \(b<-21/50\) in core
cell 1; G8 the core reference domain \(|t|\le3/10\), \(|a|<3/5\); G9 the
collar \(3/10<|t_c|<1/2\); G10 \(|b|\le1/2\), \(q<1/8\), and \(q>1/2500\) where
\(\lambda>0\); G11 four-seed inputs \(rr\le1/4\), \(\max(A,Y)<4\), \(q<1/8\);
G12 cap radius \(<2\cdot10^{-6}\) and \(|h|<1/10\) on the cap (C46); G13 the
radial maps see only radii \(>10^{-8}\) and nonzero raw outputs; G14
(11.11) and §11.5's cap estimate; G15 \(q>0\) in the collar (from G10);
G16 \(b_0=b_i=0\) exactly on pure-middle data (C42a).

**Theorem (direct certificate).** Under these hypotheses:

1. *Well defined.* \(G_{\mathrm{cl}}\) is a diffeomorphism
   \(S^4\times T^2\to Y\), so \(J_{\mathrm{cl}}\) is a smooth
   \((1,1)\)-tensor field, independent of the chart.
   - No guard on the path is ever triggered (C48). `elementary` is True
     in every executed frame (C50h), so the legacy arms are dead.
   - Every argument of a division, square root, `log1p`, `abs` or
     non-integer power lies where that operation is smooth (C50f), and
     every implicit branch cut and locally constant integer is classified
     (C50b, C50c).
   - Every one of the 111 distinct tests on the path (142 occurrences) is
     static, a guard, a precision shortcut, an identical rewrite, a flat
     join, a chart choice or a certified seam (C50), and no selector is
     hidden in another construct (C50a).
   - The genuine seams are:
     - the pole cap and slab (N9, C46);
     - the phase band and its hemisphere charts (C47);
     - the argument charts (C38, C43);
     - the cusp/middle switch (C41) and the core gauge seam \(a=0\)
       (C41b's identity, for all \(t_c\neq0\));
     - the pure middle, including its \(|t|=1\) selector and argument and
       logarithm cuts (C42);
     - the collar logarithm cut (C34, whose identities are generic in
       \((\operatorname{Re}t,\operatorname{Im}t,\rho,H,\mathrm{one},\mathrm{four})\)
       along the whole cut, with the restoration in C50g);
     - the final-phase \(\operatorname{Arg}\) cut (C42d in the middle,
       C50g in the collar);
     - the collar/middle joins (C44);
     - the height moment and height cells (C45).

   - The central fibres \(t_c=0\) are not seams: they lie in the core,
     where the reference has no division by a vanishing coordinate (C50e).
   - Across each seam the restored outputs in the two charts differ by one
     of the paper's holomorphic transitions: \(\Phi\), \(\Theta\), the cusp
     deck \(T\), or a lattice translation.
   - The differential is invertible at every point (C49).

   This also discharges the premise "smooth final fold vector" of
   C40–C44. Every seam of the projector \(P(n,h)\) is certified on the whole
   meridian domain \(S^2\times[X_-,X_+]\), \(X_\pm=\pm(1-10^{-9})\),
   including the band \(1-|h|<\alpha/2\) where the fold and the mismatch
   evaluate it (the certificates are generic in \(|t_c|<1/2\)). So \(P\) is
   smooth there, and each \(P_j=P(n,X_-+j(X-X_-)/N)\) is smooth in
   \((n,X)\). All denominators exceed \(1/2\) (G2), so \(E\) is smooth by
   induction over the fixed count (C30, section (k) of
   `symbolic_identities.py`). \(E_c\) is then smooth by G3, G4, G14 and
   C47b–e.
2. *\(J_{\mathrm{cl}}^2=-\mathrm{Id}\).* In each chart, \(J=(Dh)^{-1}J_0Dh\)
   with \(J_0^2=-\mathrm{Id}\).
3. *\(N_{J_{\mathrm{cl}}}=0\).* Let \(h\) be \(C^2\) near \(p\) with \(Dh(p)\)
   invertible. A constant linear change of the domain coordinates
   normalizes \(Dh(p)\) to the identity without changing the Nijenhuis
   equations (tensoriality, a written step). In that normalization C30i
   shows exactly that all 36 entries of \(J^2+\mathrm{Id}\) and its
   derivative, and all 216 Nijenhuis components of \((Dh)^{-1}J_0Dh\),
   vanish for an arbitrary symmetric second jet. C10a–b is the same
   identity in the forward convention \(A J_0A^{-1}\), and C10c is a
   non-coordinate control with \(N\neq0\), so the cancellation is not
   vacuous. These are exact polynomial identities in the first and second
   derivatives (type (a)). The direct circuits (l), (g), (h), (f), N8 and
   the pure-middle check apply them to the actual \(h\) on chosen cells.

T16 asked for a circuit computing the second jet of each actual branch.
Items 1–3 replace that per-branch computation. C30i holds for every
\(C^2\) map \(h\) with invertible \(Dh\), so no branch-by-branch jet
reduction is needed once each \(h\) is known to be smooth (C38–C50) with
invertible differential (C49).

**What is computed exactly.** Every identity named above is checked by
`sym/chart_identities.py` (C1–C50), `sym/node_remainder.py`,
`sym/lemmas_1d.py`, `sym/i2_margins.py` and `symbolic_identities.py`, each
printing `ALL PASS`. Inequalities are manifest (type (b)) or have explicit
constants (type (c)), written next to the check that verifies their
arithmetic. The reasons in C50f's table are types (b)/(c) and cite the
bounds they use. Numerical evaluations are labelled observations (d) and
are not used in the proof.

**Limits.** The counts (11.14) are fixed integers with 8001 and 16004
decimal digits. The formula is finite and explicit, and every claim above
is proved by identities and induction over these counts, not by
evaluating them. A practical evaluation of \(J_{\mathrm{cl}}\) at a point
of the slab is not feasible with these counts; smaller counts give
different maps, which are not certified. The checks of C48 and C50 are
censuses of the pinned text: any change to `cusp.py` invalidates them until
they are redone. The written hypotheses, the change of route
and the global statement require human review.

## 12. Direct matrix lift campaign (T17; practical radial gluing)

The collaborator permits a different explicit pullback of the same Y. The
current design retains §11's base, inverse, pole patches and torus restoration,
and constructs a direct circle lift. The results below do **not** replace the
runtime implementation of (11.15). Section 12.3 glues the local maps of §12.2 over
the pure middle; §12.4 supplies separate lifts on both punctured cusps.
Section 12.5 matches the existing pole formulas. Sections 12.7–12.9 give
a certified finite radial phase and global gluing for review. Section 12.10
assembles a value/inverse prototype and evaluates a specific tensor;
§12.11 supplies its finite jet circuit and source binding. The resulting
global certificate is a mathematical statement for human review; floating-point
checks and timing observations are kept separate from it.

### 12.1 A precise algebraic target

Let Π=ξξ* be the standard Hopf projector, |ξ|=1, and let P be (11.1),
evaluated at the corresponding base point. For a specified smooth elementary
10×2 matrix A(n,X), put

  B=P A Π,       q=tr(B*B),       L=B/√q.                  (12.1)

If q>0, the identities B*B=qΠ and BB*=qP give the circle map
ξ↦Lξ and its inverse e↦L*e. To derive the identities in any dimensions,
write B=c e f* in arbitrary nonzero frames: q=|c|²(e*e)(f*f).
C56a checks this factorization with independent symbolic frame components;
C56b checks the normalization and inverse (type (a)). Thus the remaining
design problem is a **specified small A with a proved positive margin**,
including exact agreement with I₀ and I∞ on the required pole collars.
The matrix dimensions alone impose no useful bound on the complexity of A.

For two normalized such maps L₀,L₁ over the same base, define
γ=tr(L₁*L₀). If Re γ≥0, their convex blend L=(1-w)L₀+wL₁, 0≤w≤1,
has norm square on the Hopf line

  (1-w)²+w²+2w(1-w)Re γ ≥ 1/2.                            (12.2)

Normalization then gives another circle isomorphism with adjoint inverse.
This is conditional: the transition inequality must be proved for the
actual local maps before using the blend. C56b checks the square identity
(1-w)²+w²=1/2+2(w-1/2)² (types (a)/(b)). Arbitrary local lifts cannot
be glued by an unchecked convex combination.

### 12.2 Two direct local matrices covering the pure middle

Here 1/2≤r=|t|≤2, so P has only its last two angular slots. Write the
raw S³ coordinates from §11.4 as

  u₁=a exp(i(θ+φ)), u₂=b exp(iφ), a,b>0, a²+b²=1,
  p=(a sin(θ+φ), -b sin φ, 1+a cos(θ+φ)-b cos φ),
  R=|p|, n=p/R, z=u₂/b,
  v₀=(1,0), v₁=(i(a+ib)(z-1),z+1)/2.                     (12.3)

Angles describe the proof; v₁ uses the single-valued z and has no argument
cut. Inverting the source stereographic formula gives
R²=2(1+Re u₁)(1-Re u₂) and β=R/(1+Re u₁),
X=(β-1)/√(1+β²). Its direction is exactly n in (12.3), with the source
sign convention n=(Im u₁,-Im u₂,1+Re u₁-Re u₂)/R.
C57a checks the squared-radius and frame identities:

  |v₁|=1,       Hopf(v₁)=(-a sin φ,-b sin φ,-cos φ).        (12.4)

The profile is |u₁/u₂|²=r²(A+Br²)/(1+Ar²), A=(4π)², B=e⁴.
Its derivative with respect to r² has numerator A+2Br²+ABr⁴>0.
Since A>144 and 16<B<55, its endpoint values give

  1<|u₁/u₂|²<1456/577,
  7/10<a<17/20,       1/2<b<3/4,       0<R<3.              (12.5)

For the upper endpoint, 4(A+4B)/(1+4A) increases with B and decreases
with A when B>16. The a,b bounds follow by dividing by 1+|u₁/u₂|²;
R²<2(37/20)(7/4)<9. C57c–d check the derivative, endpoints and rational
comparisons. The elementary bounds π>3, 16<e⁴<55 are also used in C53.

Let η₀,η₁ be the ninth and tenth standard coordinate vectors of C¹⁰
(indices 8,9 in Python). Define **two local candidates**

  A₀=η₀v₀*, A₁=η₁v₁*, B_j=P A_j Π,
  U₀={|Arg t|<5π/8}, U₁={|Arg(-t)|<5π/8}.               (12.6)

The two open angular sets cover the annulus, including either argument
cut. On U_j the corresponding source angular weight d_j is >1/2,
since its flat-step argument is <1/2. The same weak weight bound holds
on the closed angular domains. Therefore

  η_j*Pη_j=d_j²/(d₀²+d₁²)≥1/8,
  q_j=tr(B_j*B_j)=(η_j*Pη_j)(v_j*Πv_j).

The following explicit estimates prove

  v₀*Πv₀>1/60000 on cos θ≥-3/4,
  v₁*Πv₁>1/10000 on cos θ≤1/2.                            (12.7)

These domains contain the respective closed domains of (12.6):
cos(5π/8)=-√(2-√2)/2>-3/4 and
cos(3π/8)=√(2-√2)/2<1/2. Consequently

  q₀>1/480000 on U₀,       q₁>1/80000 on U₁.              (12.8)

Here is the complete type-(c) proof of (12.7). Put s=sin φ, c=cos φ,
D=p₃, ε=1/20, c_*=399/400. For a unit Hopf vector m,
v*Πv=(1+n·m)/2=|n+m|²/4. For v₀, if D≥0 the bound is 1/2.
If D<0, then c>0 because D≥1-a-bc. When |s|≥ε,

  v₀*Πv₀≥n₂²/4=b²s²/(4R²)>1/57600.

Otherwise c>c_*. The alternative cos θ≥0 would imply
D≥1-(17/20)ε-3/4>0. Thus -3/4≤cos θ<0 and |sin θ|>3/5,
so |p₁|>(7/10)((3/5)c_*-ε)=:h₀ and
v₀*Πv₀≥n₁²/4>h₀²/36>1/60000.

For v₁ use m from (12.4). If |s|≥ε, the second components of n and m
have the same sign and |n₂+m₂|=b|s|(1+1/R)>b|s|, giving
v₁*Πv₁>1/6400. Otherwise |c|>c_*. If c<0, then D≥1-a>0;
n₃ and m₃ are positive and |n₃+m₃|>|c|. If c>0 and cos θ≤-9/10,

  D≤1-((9/10)(7/10)+1/2)c_*+(17/20)ε<0;

now n₃ and m₃ are both negative and the same estimate applies.
Both cases give norm square >c_*²/4. Finally, if
-9/10<cos θ≤1/2, then |sin θ|>2/5, and

  |p₁|>(7/10)((2/5)c_*-ε)=:h₁,
  |n₁+m₁|>h₁/3-(17/20)ε>1/30.

The norm square exceeds 1/3600. These exhaustive cases prove (12.7).
C57b checks the projection formula, and C57d checks every rational
comparison in this argument. C56 supplies each normalized local inverse.

**Measured margins, type (d).** C57e uses 17 radial, 64 θ and 64 φ values
(69632 base points); each closed angular domain contains 43520 samples.
The minimum q₀ and q₁ are 0.06960246379 and 0.02865159099, with measured
excesses over (12.8) greater than 0.0696003 and 0.0286390. C57f separately
evaluates the actual `elliptic_raw_base_reference` and
`elliptic_bundle_embedding` at 192 base points, then applies B_j and its
adjoint inverse on each applicable patch (120 checks each). At 60 digits,
the maximum profile/image/norm/inverse residual is <2.2·10⁻⁶⁰.
These samples check the transcription; they are not used in (12.7).

### 12.3 A direct lift and inverse on the whole pure middle

The transition is essential. Put h=v₀*Πv₁ and k=i(a-ib), so |k|=1.
For L_j=B_j/√q_j from (12.6), the source's two angular entries give

  tr(L₁*L₀)=h/|h|                         when θ>0,
  tr(L₁*L₀)=exp(-2πi x₄) h/|h|           when θ<0.          (12.9)

Indeed tr(B₁*B₀)=P_(9,8) h, using zero-based indices, and
P_(9,8)/|P_(9,8)| is respectively 1 or exp(-2πi x₄).
This follows directly from (11.1); h is nonzero on the overlap by
the rank-one identity |h|²=(v₀*Πv₀)(v₁*Πv₁).

Use the real θ₁=Arg(-t)+π on U₁ and the logarithm ℓ₁=log r+iθ₁.
In this paragraph ξ₄(t,ℓ₁) denotes the fourth real lattice coordinate
of (πs(s-1),2ts), s=ℓ₁/(2πi), using that same logarithm in the
lattice. Define

  λ=1-smooth_step((r-6/5)/(2/5)),
  C(r,θ₁)=(1-λ)θ₁/(2π)+λ ξ₄(t,ℓ₁)+λ/2,
  w=smooth_step((θ₁-5π/8)/(3π/4)),
  d=conj(k) exp(i w(π-2π C(r,θ₁))).                        (12.10)

Here `smooth_step` is the existing flat step. Its name distinguishes it
from the lattice variable s. C42's lattice identity gives
ξ₄(θ+2π)=ξ₄(θ)+1, hence C(θ+2π)=C(θ)+1. The raw source formulas give
x₄=φ/(2π)+C(r,θ) modulo integers on the principal θ chart. For r≤1
this is ξ₄+1/2; for r>1 substitute Arg u₁=θ+φ modulo 2π into
the source formula. C58c checks the substitution and the shift, including
the harmless integer difference at r=1. Since λ=1 up to r=6/5, C is
one smooth formula across r=1.

Put χ=smooth_step((|Arg t|-3π/8)/(π/4)). Define

  M=(1-χ)L₀+χ d L₁,       L_mid=M/√tr(M*M).                (12.11)

Evaluate L₀ only when χ<1 and d L₁ only when χ>0. The corresponding
closed domains satisfy (12.8). This is a fixed finite formula: no
projection chain, quadrature, phase walk or root solve is present.

**Nonvanishing certificate.** On the upper overlap 3π/8≤θ≤5π/8,
w=0, and on the lower overlap -5π/8≤θ≤-3π/8, w=1. Thus the real
part of the corrected transition is respectively
Re(kh)/|h| and Re(-k exp(-iφ)h)/|h|. Both exceed 1/48, as proved
below. Equation (12.2) therefore gives tr(M*M)≥1/2. It follows that
L_mid*L_mid=Π and L_mid L_mid*=P; its circle inverse is L_mid*.
These are type-(a)/(b) consequences of C56 and the following certificate.

Write c=cos φ, s=sin φ, Cθ=cos θ, Sθ=sin θ, D=p₃. C58a checks

  4R Re(kh) = F(R),
  F(R)=(R+D)(1-c)+s²
       +a(Sθ c+(Cθ+1)s)(b(1+c)-a s).                     (12.12)

On the upper band |Cθ|<2/5 and Sθ>9/10. Also
R²>2(3/20)(1/4)=3/40>(27/100)² by (12.5), so
F(R)≥F(27/100). The following finite polynomial identity certificate
proves F(27/100)>1/4 on the larger parameter box

  a∈[7/10,17/20], b∈[1/2,3/4], Cθ∈[-2/5,2/5], Sθ∈[9/10,1].

For σ=±1 substitute c=σ(1-y²)/(1+y²), s=σ 2y/(1+y²), |y|≤1.
These two rational charts cover the entire φ circle. The expression

  Q_σ=(1+y²)²(F(27/100)-1/4)

is a polynomial. On each y-cell [-1+j/4,-1+(j+1)/4], j=0,…,7,
rescale all five variables to [0,1]. Expand Q_σ in the tensor Bernstein
basis of degrees (2,1,1,1,4): the basis functions are products of
binomial(d,i) x^i(1-x)^(d-i), hence are nonnegative and sum to one.
For power coefficients p_α the exact Bernstein coefficients are

  b_β=Σ_(α≤β) p_α Π_l binomial(β_l,α_l)/binomial(d_l,α_l).

C58b computes all 120 rational coefficients in each cell and converts
them back to check equality of the entire polynomial. Their minima are:

| σ | minima for j=0,…,7 |
|---|---|
| +1 | 29971/32000, 2939/16000, 18459/128000, 32491/128000, 3019/8000, 67337/128000, 3341/4000, 219087/128000 |
| -1 | 140259/64000, 5081/2560, 15943/8000, 140141/64000, 261/100, 791509/256000, 14457/4000, 881/250 |

Every coefficient exceeds 1/8. This is a manifest positive combination
(type (b)), verified by polynomial identities (type (a)); no numerical
sampling or interval arithmetic supplies the sign. Thus
Re(kh)>1/(16R)>1/48. For the lower band, C58a checks
h(-θ,-φ)=-exp(-iφ)h(θ,φ). It reduces its corrected real part to
the same positive upper-band expression. Finally |h|≤1 because Π is
an orthogonal projector and v₀,v₁ are unit. This proves the normalized
transition bound used above.

**Smoothness and scope.** The source P is smooth through the middle
argument and lattice selectors by C42. Formula (12.3) for v₁ is periodic
in φ. Near the θ=0 cut of θ₁, χ=0 on an open set; near the θ=±π cut
of Arg t, χ=1 on an open set and only the θ₁ formula is used. Flatness
of the two weights gives smooth joins at their plateau endpoints.
All denominators are positive by (12.8) and (12.11), so (12.11)
defines a smooth circle-bundle isomorphism throughout 1/2<|t|<2,
with the same formulas and bounds at the two boundary tori.
This proves a regional lift, not an extension to the cusps.

**Implementation and observations (d).** C58's `direct_middle_matrix58`
evaluates (12.11) using the actual source projector. At 324 base points
(three radii, nine θ values, twelve φ values), the minimum blend norm
square is 0.6904810986 and the minimum evaluated corrected overlap is
0.3035001081. Their measured margins over 1/2 and 1/48 are >0.19048
and >0.28266. Both inverse matrix identities have residual <3.2·10⁻⁶¹
at 60 digits. The 84 comparisons across the θ and φ argument cuts
have maximum residual <6.0·10⁻⁶⁰. The finite coefficient table is used
only for proof, not during evaluation of the lift.
C58e times 100 lift-value evaluations at r=1, θ=π/2, φ=j/100,
j=0,…,99, at 60 digits, including the raw base and source projector.
A separate run of this same loop took 0.1772 seconds on the working
machine (type (d)); this is not a benchmark of J or its derivatives.

**Next gate.** Extend the direct map across the two cusp collars,
including the central fibres, and match the existing pole maps exactly.
Then certify stable derivative evaluation and build the final
source-linked J circuit. No new global tensor or global diffeomorphism
is claimed by this middle-region result.

### 12.4 Direct lifts on both cusps, including the nonnode central fibres

This section gives separate regional maps, not their radial gluing to
L_mid. Use the existing raw base, 0≤r=|t_c|<1/2, and exclude the node.
Put ν=1 or e⁻², m=ν exp(-2r²/(1-r²)), φ=2πH, with
-1/2<H<1/2 on the first height chart. In the south write
u₁=a exp(i(θ+φ)), u₂=b exp(iφ); in the north write
u₁=-a exp(-iφ), u₂=-b exp(i(θ-φ)). Here θ=arg t_c when r>0.
The formulas below use u₁,u₂ themselves and extend to r=0.
In both cases p=(Im u₁,-Im u₂,1+Re u₁-Re u₂), R=|p|, n=p/R.

Let U_j be the existing unitary pole matrix, and put F=U_j(m conj x,y).
C59a checks the exact Hopf vectors

  Hopf(F)=(2m Im t_c,-f,2m Re t_c)              (south),
  Hopf(F)=(f,2m Im t_c,2m Re t_c)               (north),
  f=m²|x|²-|y|².                                        (12.13)

On the first height chart f has the sign of H by the strictly increasing
odd moment (11.4). Thus F's normalized second component (south), or first
component (north), has the same strict sign as the corresponding component
of n whenever H≠0. These vectors cannot be antipodal. At H=0 and r>0,
their remaining sine components have the same sign. If sin θ=0, their
third components have the same strict sign: 1+a-b>0 at θ=0 in the south,
1-a+b>0 in the north, and 1-a-b<0 at θ=π in either case. The last
inequality uses a,b>0 and a²+b²=1. At r=0, H≠0 the signed-height
argument still applies. Therefore

  (F/|F|)*Π(F/|F|)=|n+Hopf(F)/|F|²|²/4>0.              (12.14)

This is a sum-of-squares nonvanishing proof (type (b)), using the checked
identities (12.13); it gives strict positivity, not a uniform constant.

Let E_{j,0}:C²→C¹⁰ insert into the first pair of cusp j, and let η_{j,1}
be the first coordinate of its second pair (indices 2 and 6). Define

  A_c,0=E_{j,0} diag(m/ν,1) U_j*,
  A_c,1=η_{j,1} v(H₁)*,
  H₁=1/2+Arg(-z)/(2π),
  τ=(π/2)smooth_step((H₁-7/20)/(3/10)),
  v(H₁)=U_j(cos τ,sin τ),                               (12.15)

where z=u₂ in the south and z=-conj(u₁) in the north. H₁ lies in (0,1)
on the second height chart. The maps B_c,j=P A_c,j Π, normalized as in
(12.1), are denoted L_c,0 and L_c,1. Use the first on |H|<7/20 and
the second on 3/10<H₁<7/10. These two domains cover the height circle.
Their source weights are nonzero; the second height weight is exactly 1.
The first map is nonzero by (12.14).

Here are quantitative bounds for the second map and their derivation
(type (c), rational identities/comparisons in C59c). The profile is
increasing. At r=1/2 its south ratio square is (A+B/4)/(A+4),
A=(4π)²>144, 16<B=e⁴<55. This decreases in A and increases in B;
its value at A=144,B=55 is <5/4. Hence a<3/4 and b>2/3.
In the north the ratio square is (A+1/4)/(A+4B)<1, giving b<3/4
and a>2/3. Also R<3 by the triangle inequality. On the second height
domain cos φ<0, so p₃>1/4 and n₃>1/12. C59b checks

  Hopf(v)=(0,-cos(2τ),sin(2τ))             (south),
  Hopf(v)=(cos(2τ),0,sin(2τ))              (north).

Since 0≤τ≤π/2, Hopf(v)₃≥0. Thus v*Πv>1/576>1/600.
The chosen source second-pair representative has positive height
H₁∈(3/10,7/10), so m²|x₁|²>|y₁|², also across H₁=1/2 by the
explicit endpoint chart in §11.2. As ν≥m, its first squared normalized
component exceeds 1/2. Using V*V≤2 from (11.2),

  tr(B_c,1*B_c,1)>β_j²/2400>0.                           (12.16)

In particular this exceeds 1/2400 in the pure cusp r≤3/10.
No positive constant independent of r is asserted as r approaches 1/2,
where the source radial weight β_j vanishes.

**Height overlap, without a phase walk.** On 3/10≤H≤7/20, τ=0
and v=U_j(1,0). Choose the source circle gauge x>0, and write
(m conj x,y)=(A₁,B₁ exp(iθ)), A₁>B₁≥0. In the U_j frame put
Q=U_j*ΠU_j. Its first diagonal entry is (1+ε)/2, where
ε=-n₂>0 in the south and ε=n₁>0 in the north. Its off-diagonal
entry has magnitude √(1-ε²)/2. The source slot of A_c,1 contributes
a positive factor νx. Consequently the phase of tr(L_c,1*L_c,0)
is the phase of (A₁,B₁ exp(iθ))*Q(1,0), whose real part is at least

  [A₁(1+ε)-B₁√(1-ε²)]/2 ≥ A₁ε/2.

Its modulus is at most √(A₁²+B₁²)≤√2 A₁, so the normalized real
part is ≥ε/(2√2). The inequality 1-√(1-ε²)≥0 is also the manifest
identity ε²/(1+√(1-ε²))≥0, checked by C59c.
On -7/20≤H≤-3/10, τ=π/2 and v=U_j(0,1). Choose y>0.
Now B₁>A₁, and the second source chart has x₁=1/(my)>0 by the
explicit deck map. Interchanging the two components gives the same bound.
On either overlap sin(2π|H|)>4/5. For example its minimum is sin(3π/10),
whose square is (3+√5)/8>16/25. With the amplitude bounds above and R<3,
ε>8/45; since 2√2<3, this proves

  Re tr(L_c,1*L_c,0)>1/20.                               (12.17)

Set χ_c=smooth_step((|H|-3/10)/(1/20)) and

  M_c=(1-χ_c)L_c,0+χ_c L_c,1,
  L_c=M_c/√tr(M_c*M_c).                                  (12.18)

Evaluate only the terms with nonzero weights. Equations (12.2) and
(12.17) give tr(M_c*M_c)≥1/2. Thus L_c*L_c=Π, L_c L_c*=P,
and the inverse is L_c*. These are consequences of C56 (type (a)/(b)).
The first argument cut is surrounded by χ_c=1; the second by χ_c=0.
The source projector is smooth across the central axes and deck seams by
C41/C45. The matrices (12.15), positive denominators and flat weights
therefore give a smooth map on each cusp with its node removed.

**Executable observations (d).** C59d uses the actual source projector
at 1232 points: both cusps, six radii (including zero and two collar
radii), eight θ values and thirteen heights, excluding the two nodes.
The minimum local norm square is 0.0001323368056; the minimum blended
norm square is 0.9753386158 and the minimum overlap is 0.9506772316.
The latter two exceed the proved bounds by >0.47533 and >0.90067.
The two inverse matrix identities have residual <3.7·10⁻⁶¹; the 96
height-seam comparisons have residual <4.6·10⁻⁵⁹, at 60 digits.
The local minimum includes collar points, so it is not compared to
the pure-cusp constant in (12.16).

The remaining geometric tasks are radial gluing to §12.3 and exact
agreement with the existing pole lifts before applying (11.17).
Regional circle inverses do not establish a global diffeomorphism or
the requested executable J and its derivatives.

### 12.5 Exact agreement with the existing pole lifts

The first cusp matrix admits a small correction that closes the regional
pole-matching task. Put α=10⁻⁵, σ₀=-1, σ_∞=1 and d_X=1-σ_j X.
Replace diag(m/ν,1) in (12.15) by

  D_j=diag(1-λ_p(1-m/ν),1),
  λ_p=smooth_step((d_X-α)/α).                             (12.19)

It equals the previous diagonal for d_X≥2α and the identity for d_X≤α.
Here is the explicit perturbation certificate (type (c), all rational
comparisons checked by C60a). Throughout the support d_X<2α, the existing
linear pole coordinate satisfies

  r_L=|z_L|²≤2d_X/3<4α/3<10⁻⁴,
  r_S=|p_S|≤r_L+80r_L²<1/70000,
  |t_c|≤r_S/(2ν)<4r_S<10⁻⁴.

These are (11.9) and (11.17); ν≥e⁻²>1/8 follows from e⁴<55<64.
Also |H|≤r_S+38r_S³<1/10. Thus χ_c=0 here and the source embedding
is entirely in the first cusp pair. The existing raw pole projection
has squared norm at least 1-10000r_L²>(99/100)². Meanwhile

  ||D_j-I||≤1-exp(-2|t_c|²/(1-|t_c|²))
           ≤2|t_c|²/(1-|t_c|²)<10⁻⁷.

The inequality 1-e⁻ˣ≤x for x≥0 follows by integrating
1-e⁻ˣ≥0, or differentiating x-1+e⁻ˣ. Multiplying by orthogonal
projectors and unitary matrices cannot increase this operator norm;
the input Π has rank one, so it is also the needed Frobenius bound.
Consequently the corrected matrix satisfies

  tr(B_c,0*B_c,0)>(99/100-10⁻⁷)²>0.98.                  (12.20)

Outside this band (12.14) supplies positivity. Flatness gives smooth
joins between the formulas.

On d_X≤α, (12.19) gives B_c,0=P I_j Π exactly, where
I_j=E_{j,0}U_j*. Its normalization is the same minimal raw pole lift
g_j used in (11.10). Therefore use the existing explicit (11.17) for
d_X<α/2, and the corrected cusp formula elsewhere. On
α/3<d_X<α/2 one has r_L>α/21>2ε, ε=10⁻⁸, hence H_j=g_j.
The formulas agree on this open overlap, with all derivatives. The
elementary inverse and positive differential of the pole replacement
are already certified in §11.3; neither is changed here.

C60b evaluates forty actual-source matrices: both pole signs, four
distances and five Hopf directions, including axis directions. Their
minimum norm square is 0.99999999850476, a measured excess >0.01999999
over (12.20). The inverse identities and equality to the existing inner
raw lift have maximum residual <6.5·10⁻⁶² at 65 digits (type (d)).

This supplies capped regional lifts with elementary inverses. It does
not join them to the middle lift. A global radial transition, followed
by the actual tensor and derivative circuit, is still required.

### 12.6 Radial overlaps and a specific remaining phase inequality

The middle lift extends to 49/100≤r=|t|≤100/49 by correcting its
phase offset in the two short cusp collars. Write φ=Arg u₂, use the
same θ₁ as in (12.10), and let ρ(H,t_c,ν) be (11.5), including its
deck extension. In the south collar r<1/2 replace C by

  H=φ/(2π), P=log r-2r²/(1-r²),
  C_s=ξ₄(t,log r+iθ₁)+ρ(H,t,1)/P-H.                     (12.21)

In the north collar r>2 put r_c=1/r and replace C by

  H=(π-θ₁-φ)/(2π), P=-2+log r_c-2r_c²/(1-r_c²),
  C_n=1-ρ(H,1/t,e⁻²)/P-φ/(2π).                         (12.22)

C15's exact slice identities give x₄=ξ₄+ρ/P in the south and
x₄=-ρ/P in the north. Thus the lower angular transition still uses
exp(-2πix₄)=exp[-iφ-2πiC]. C61b checks that C is unchanged under
φ→φ+2π and increases by 1 under θ₁→θ₁+2π, using ρ(H+1)=ρ(H)+P.
At either collar endpoint ρ=P(H+1/2). Equation (12.21) then equals
ξ₄+1/2, and (12.22) equals θ₁/(2π), exactly as in (12.10).
The **+1 in (12.22) is necessary**: dropping it would change the
interpolated factor exp[-2πiwC] even though exp[-2πiC] is unchanged.
C45's flat radial matching makes all endpoint derivatives agree.

The bounds 7/10<a<17/20, 1/2<b<3/4 and R<3 still hold on this
larger annulus. C61a checks the explicit endpoint inequalities

  profile(49/100;144,16)>49/51,
  profile(100/49;144,55)<289/111,
  profile(r;A,B)=r²(A+Br²)/(1+Ar²).

For the lower bound first decrease B to 16; the expression then
increases in A at r=49/100, since 16r⁴<1. For the upper bound increase
B to 55; it decreases in A at r=100/49. Combine with radial
monotonicity. These are type-(c) bounds, with the derivative identities
already checked in C57/C53. The source middle weight is ≥1/2 throughout
the enlarged annulus. Therefore the local norm-square bounds in
(12.8) may be divided by 4, and the same angular phase certificate
still gives corrected real overlap >1/48 and blended norm square ≥1/2.
This proves actual radial overlap domains for the middle and cusp maps.

At a common base point define δ=tr(L_mid*L_c). Both maps are unit
isomorphisms between the same complex lines, so |δ|=1 and
L_c=L_mid δ (type (a), the rank-one algebra of C56). At fixed south
radius r=49/100 call this δ_s(θ,φ); at fixed north radius
r=100/49 call it δ_n(θ,φ), using the same global θ=Arg t and
φ=Arg u₂ in both. These are single-valued smooth functions on T²
by the seams already checked above. C61's `frozen_transition61`
evaluates their explicit finite formulas.

**Remaining proof gate.** A particularly short global gluing would
follow from the specific inequality

  Re g(θ,φ)>0 for every (θ,φ)∈T²,
  g=δ_n/δ_s.                                             (12.23)

This inequality is **not proved**. Conditional on it, the elementary
real phase Θ=atan(Im g/Re g) is smooth and periodic, with exp(iΘ)=g.
Use a radial flat step w_r=smooth_step((r-3/4)/(1/2)). In the region
0<w_r<1 multiply L_mid by δ_s exp(iw_rΘ). This equals δ_s on
the lower plateau and δ_n on the upper plateau.

To match the actual radial dependence, rather than only the frozen
boundary values, on the south end use δ_s evaluated at

  R_s(r)=49/100+(r-49/100)
                 [1-smooth_step((r-49/100)/(1/200))].

It agrees with r to all orders at 49/100 and is exactly 49/100
for r≥99/200, staying inside the common collar. At the north end
use the same retreat in r_c=1/r. Both retreats are already constant
on a neighbourhood of 0<w_r<1. Thus the plateau formulas join
smoothly, with no phase walk. This is a conditional construction,
not an assertion of (12.23) or a completed global map. Its algebraic
inverse would again use the adjoint circle map and the existing
two-way base and cap formulas.

**Measurements, type (d).** C61c evaluates both frozen transitions
on a 32×32 grid (2048 source points). The middle blend norm square
is at least 0.6859041071 and its corrected angular overlap at least
0.3364921174, with margins >0.18590 and >0.31565 over their proven
bounds. Inverse and scalar-transition residuals are <3.6·10⁻⁶¹;
32 seam pairs have residual <2.4·10⁻⁶⁰ at 60 digits.
C61d records the separate grid values of g: its minimum real part is
about 0.227442 and maximum principal-argument magnitude about 1.34135.
All sampled sums of adjacent principal phase increments around either
grid circle are zero to floating-point precision. These observations
do **not** prove (12.23), a continuous winding number, or a global
argument-cut margin. A proof of (12.23), or another certified radial
phase construction if it fails, is the next geometric task. The
source-linked executable J and its derivatives remain to be built.

### 12.7 Cancelling the source embedding and freezing at the exact boundary

There is an alternative to proving (12.23). It uses the exact radii
`r_c=1/2` at both ends, a two-component scalar comparison, and a finite
phase sum with **65536 terms**. The proof and value evaluator below do
not require a sign for Re g on the whole torus. The earlier, stronger
inequality (12.23), at the different radii 49/100 and 100/49, remains open.

First cancel the positive source weights in the cusp matrices. For
49/100≤r_c≤51/100 put t_c=r_c exp(iθ_c), ν=1 or e⁻²,
m=ν exp[-2r_c²/(1-r_c²)] and P=log(mr_c). Use the existing height
ρ for r_c<1/2 and set ρ=P(H+1/2) for r_c≥1/2. C45 proves flat
matching of these formulas. Let y=exp ρ, x=t_c/y, and let e be the
unit source embedding of this physical reference, expressed in the
middle chart. With U_j as in (12.15), define unit rows

  c₀ = F*Π / ||F*Π||,             F=U_j(m conjugate(x),y),
  c₁ = exp(iη) v*Π / ||v*Π||,     η=θ_c if H≥0, otherwise 0,
  v=U_j(cos τ,sin τ),             τ as in (12.15),
  c = [(1-χ_c)c₀+χ_c c₁] / ||(1-χ_c)c₀+χ_c c₁||.         (12.24)

Then L_c=e c. Only rows with nonzero weights are evaluated. Below
1/2 this is exactly (12.18): in the first source pair the factor
`diag(m/ν,1)` turns its row into F* times a positive scalar. In the
second pair the first source coordinate has phase conjugate(x)/|x|
for H≥0, and phase 1 for H<0, since the shifted x is 1/(my)>0.
The projector column conjugates that phase. Normalization cancels
all the remaining positive factors. Thus (12.24) extends the map
even where those source slots have zero radial weight; it does not
divide by a zero slot. This is the outer-product identity C62b (a).

The positivity proof of §12.4 extends to this small collar. The
profile bounds in §12.6 still hold there. In cusp coordinates the
dominant amplitude is >2/3 and the other amplitude is <3/4:
in the south, b²=1/(1+profile(r_c;A,B))>4/9, using
profile(51/100;144,55)<5/4; in the north,
|u₂/u₁|²=r_c²(r_c²+A)/(Ar_c²+B)<1, since r_c⁴<B,
so a>1/√2>2/3 and b<1/√2<3/4. The south bound also gives
a<√5/3<3/4. The first height has the sign of H; the second-frame
hemisphere argument and overlap estimate (12.17) therefore apply
unchanged (b,c). At the height cut e changes by exp(-iθ_c), while
c₁ changes by exp(iθ_c), so their product agrees exactly. The
existing logarithm transitions give the angular seam. Consequently
the extended L_c is smooth throughout this collar.

Here are explicit lattice formulas that remove a matrix solve from
the comparison. Write t=r exp(iθ), ℓ=log r and P₀=ℓ-2r²/(1-r²).
For (z,w)=(-iρ,0),

  x₄=ρ/P₀,                    x₁=-θρ/(2πP₀).

For Ξ=(πs(s-1),2ts), s=(ℓ+iθ)/(2πi),

  ξ₄=θ/(2π)-ℓ/(2P₀),
  ξ₁=-(θ²+ℓ²)/(8π²)-θ/(4π)
        -ℓr²/[2π²(1+r²)]+θℓ/(4πP₀).                    (12.25)

C62a verifies these by substituting into the full complex lattice
equations, using the rotated real coordinates
u=cosθ x₂+sinθ x₃, v=-sinθ x₂+cosθ x₃:

  w/t=(1/r-r)u+i(1/r+r)v+2x₄,
  z=2πx₁-rv+θx₄-i(ru+ℓx₄).

At the exact boundary let h=log 2, c_s=1/(3h+2), k_s=h+2/3,
k_n=h+8/3. Set θ=Arg t, φ=Arg u₂, θ_c=θ in the south and
θ_c=-θ in the north. In the south H=Arg u₂/(2π); in the north
H=Arg(-conjugate(u₁))/(2π). The unit version of F in (12.24) is

  F=U_j (exp(k_j H-iθ_c), exp(-k_j H)) / √[2cosh(2k_j H)]. (12.26)

The relevant physical source coordinates are

  south: x₄=H+θ/(2π)+c_s,  x₁=ξ₁-θ(H+1/2)/(2π),
  north: x₄=-H-1/2,        x₁=θ(H+1/2)/(2π),             (12.27)

where ξ₁ in this display is (12.25) at r=1/2. The north sign uses
P₀(1/r_c)=-[log r_c-2r_c²/(1-r_c²)-2], also checked by C62a.

For the middle rows put m₀=e₁*Π/||e₁*Π||, and, with v₁ from
(12.3),

  m₁=(b-ia) exp{i[w(π-2πC)+(θ₁-θ)x₄]}
                 v₁*Π/||v₁*Π||,
  C=θ₁/(2π)+c_s (south), or θ₁/(2π) (north),
  m=[(1-χ)m₀+χm₁]/||(1-χ)m₀+χm₁||.                     (12.28)

The actual middle matrix is e exp(-2πix₁)m. Thus the boundary
transition is the two-row expression

  δ_j=exp(2πix₁) Σ_l conjugate(m_l)c_l,     g=δ_n/δ_s.   (12.29)

This g is at r=1/2 and 2, unlike (12.23). It is a smooth unit
function on T², by the regional identities and seams. C62b checks
the exact cancellation of e. `boundary_transition62` implements
(12.26)–(12.29) with two-component rows, avoiding a Hopf spinor
chart and the ten-component embedding. Its float64 outputs are
numerical evaluations of these elementary formulas, not proof data.

**Measurements (d).** C62c checks 720 actual-source points at five
cusp radii from .49 to .51. Below .5 the cancelled and original
matrices differ by <3.6·10⁻⁶¹; inverse identities have residual
<3.2·10⁻⁶¹. Sixty angular/height seam pairs have residual <8.4·10⁻⁶¹.
At the exact boundary, 144 comparisons of the float64 scalar formula
with the independent 60-digit matrix calculation differ by <9·10⁻¹⁶.

### 12.8 A practical derivative bound and an exact winding certificate

All estimates in this subsection are type (b)/(c); their algebra and
rational comparisons are checked in C63. They concern only the two
fixed boundary tori, not the derivatives of the full tensor.

**Projection denominators.** The amplitude box in §12.6 gives
R²=2(1+Re u₁)(1-Re u₂)>3/40>1/16. Hence R>1/4,
||∂θ n||<4 and ||∂φ n||<6: the corresponding derivatives of p
have norms a<1 and at most a+b≤√2<3/2. Consequently the operator
norms of ∂θΠ and ∂φΠ are <2 and <3. C63b checks
(∂Π)²=||∂n||²I/4.

The first cusp row has the additional uniform bound

  ||F*Π||²>1/576,                   ||F||=1.             (12.30)

To see this, set E=tanh(2k_jH), T=sech(2k_jH), and ψ=2πH.
C63a gives the unnormalized Hopf pairings

  south: T[a cosψ+(1-b cosψ)cosθ_c]+bE sinψ,
  north: T[b cosψ+(1-a cosψ)cosθ_c]+aE sinψ.

If |H|≤1/12, E sinψ≥0 and each bracket is at least
(a+b)cosψ-1>(6/5)(√3/2)-1>0. The projection norm square is
therefore >1/2. If 1/12≤|H|≤7/20, the distinguished component
of the target Hopf vector has the same sign as E and magnitude
ε>(1/2)(1/2)/3=1/12. Two unit vectors with these same-sign
components have scalar product at least -√(1-ε²). Thus their
spinor projection norm square is at least
(1-√(1-ε²))/2≥ε²/4>1/576. These cases cover the active first row.
The second cusp row retains ||v*Π||²>1/600 from §12.4. The middle
target row bounds are >1/60000 and >1/10000 from §12.2; source
slot weights have already been cancelled.

**Derivative budget.** For any nonzero complex row B,

  ||D(B/||B||)||≤||DB||/||B||.

C63b checks the derivative and its squared norm in four real
coordinates, with the removed term a square. For (12.26), the
H derivative has norm ≤k_j<4, and the θ_c derivative has norm ≤1.
For example the first bound reduces to
4k²z²/(1+z²)²≤k², with difference k²(z²-1)²/(1+z²)².
The elementary bounds 2/3<h<7/10 follow from
log 2=2Σ_{n≥0}(1/3)^(2n+1)/(2n+1), whose tail after 2/3 is
<1/36. This also proves 0<c_s<1/4.

On each active branch, |H_θ|,|H_φ|≤1/(2π)<1/6.
The first cusp row derivatives are both bounded by
24(5/3+2)=24(2/3+3)=88. The rotating v has angular derivative
norm ≤5/3, since s′≤2 (L1); its phase has θ derivative ≤1 and
φ derivative 0. The second row derivatives are therefore bounded
by 25(5/3+2)+1<118 and 25(5/3+3)<118. The height cutoff has
angular derivatives ≤20/π<7. Its blend has norm ≥1/√2, so

  ||c_θ||, ||c_φ|| < (118+14)(3/2) < 200.

For m₀ the bounds are 490 and 735, using √60000<245. For m₁
the normalized projection contributes <200 and <400. Its additional
phase in (12.28) has derivatives <5 and ≤1: on the support of w′,
|π-θ₁-2πc_s|≤7π/8, |w′|≤8/(3π), |2πx₄,θ|=1, and
|2πx₄,φ|=1. For the north put c_s=0 in this calculation. Thus
the phase θ derivative is ≤7/3+1+1<5. The angular cutoff χ has
|χ_θ|≤8/π<3, χ_φ=0. After normalization,

  ||m_θ||<(490+6)(3/2)<745,     ||m_φ||<735(3/2)<1103.

Finally (12.25)–(12.27) give |2πx₁,θ|≤2 and
|2πx₁,φ|≤1/2. Applying the product rule to (12.29) yields

  |∂θ g|<1900,                     |∂φ g|<2700.          (12.31)

The estimates hold up to all seams by the already proved smooth
regional matching. Inactive rows are omitted, so a zero inactive
denominator is not differentiated.

**The two winding numbers are zero.** Here is an elementary proof,
independent of the numerical phase grids. Write A_j=atan(a_j/b_j),
where j refers to the corresponding boundary radius.

On φ=0, choose the continuous target spinor
(cos(β_j/2),sin(β_j/2)), where
β_j=arg[(1-b_j)+a_j exp(iθ)] is continuous on [-π,π] with
β_j(±π)=±π. It exists with |β_j-θ|<π/2, because
Re[a_j+(1-b_j)exp(-iθ)]>a_j+b_j-1>0. The middle row's phase
μ_j is 0 for |θ|≤3π/8; on the two angular overlaps it is the
continuous argument of (1-χ)+χ exp(-iA_j). In the remaining
positive and negative portions respectively it is

  μ_j=-A_j+w(π-θ-2πc_j),
  μ_j=-A_j+(1-w)(π+θ+2πc_j),     c_j=c_s or 0.

The formulas agree at 5π/8 and -5π/8; their endpoint difference
is μ_j(π)-μ_j(-π)=-2πc_j.

In the south H=0, χ_c=0. C63e gives exactly
F=exp(-iθ/2)(cos(θ/2),sin(θ/2)), so F* times this target spinor
has phase θ/2 and positive modulus cos((θ-β_s)/2). A continuous
real phase of δ_s is therefore 2πξ₁-μ_s. Equation (12.25) gives
2π[ξ₁(π)-ξ₁(-π)]=-2πc_s (C63f); this phase has equal endpoints.

In the north H=1/2-θ/(2π) for θ>0 and H=-1/2-θ/(2π) for θ<0.
Put c=cos(β_n/2)>0, s=sin(β_n/2), A=exp(k_nH)/√(2cosh(2k_nH)),
B=exp(-k_nH)/√(2cosh(2k_nH)). The first-row scalar has

  Re = [A cosθ(c+s)+B(c-s)]/√2,
  Im = -A sinθ(c+s)/√2.                                 (12.32)

For 0<θ<π its imaginary part is negative. The second-row scalar
is a positive real multiple of exp(-iθ), so the blended phase has
a continuous value in (-π,0). For -π<θ<0, A≤B, s<0, c>0,
and c-s>|c+s|, so (12.32) has positive real part. The second-row
scalar is positive real as well. Here the blended phase lies in
(-π/2,π/2). The two values join at 0 with phase 0; at θ=π the
phase is -π and at θ=-π it is 0. The source phase
2πx₁=θ(H+1/2) has endpoint difference π, and μ_n has endpoint
difference zero. Thus δ_n also has a continuous real phase with
equal endpoints. This proves zero θ winding of g. C63d–f check
all the half-plane and endpoint identities used here.

For the other generator set θ=0. Then n₃>0 and a target spinor
is (C,(a+ib)S) with C>0 and S having the sign of sinφ. In each
active cusp row its coefficients before U_j are real; writing the
result as (A,iD) in the south and (A,D) in the north gives A>0
and DS≥0. This holds for both the exponential row and the rotating
row: the sign of D is the sign of H, which agrees with sinφ on
this circle. Their scalar projections are respectively

  AC+bDS-iaDS,                     AC+aDS+ibDS.          (12.33)

Consequently the south blended row has phase in [-atan(a_s/b_s),0]
and the north one in [0,atan(b_n/a_n)]. Here m=m₀ has phase 0.
The source phases are α₀=h(1/5-h/4)/π in the south and 0 in
the north, with 0<α₀<1/100 (C63f). The profile is strictly
increasing, so a_s/b_s<a_n/b_n and
atan(a_s/b_s)+atan(b_n/a_n)<π/2. Hence g on this circle has
phase in (-1/100,π/2), a single continuous periodic argument.
Its φ winding is zero. C63g checks (12.33). This proves both
winding statements using inequalities and endpoint identities;
no sampled winding number is an input.

### 12.9 A certified 65536-term phase and the radial gluing it supplies

For θ,φ∈[-π,π], N=65536, define g_j=g(jθ/N,jφ/N), with g from
(12.29), and

  Θ(θ,φ)=-α₀+Σ_{j=1}^N atan[Im(g_j/g_{j-1})/
                                      Re(g_j/g_{j-1})]. (12.34)

By (12.31), every adjacent chord has length less than
(22/7)(1900+2700)/65536<1/4. Since |g_j|=1, its ratio has
real part >31/32. Thus each atan is on its smooth right-half-plane
branch. C64a checks that its exponential is that ratio. The
products telescope, and g(0,0)=exp(-iα₀), so exp(iΘ)=g exactly.
This is a fixed finite elementary expression, with no implicit
inverse, quadrature, infinite series or adaptive stopping rule.

It is smooth and periodic. Indeed the sum is a continuous lift
of g on the angle rectangle. Its values on opposite edges differ
by a continuous integer multiple of 2π and therefore by a constant
multiple. At φ=0 and θ=0 those multiples are zero by §12.8. Hence
the values match. Locally any two lifts of the smooth function g
differ by a constant multiple of 2π, so their derivatives of every
order also agree. This proves smooth matching across both edges
and their corners, rather than just matching sampled values.

Here is the resulting radial correction. For 1/2≤r≤51/100 use

  R_s(r)=1/2+(r-1/2)[1-smooth_step((r-1/2)/(1/100))],

and set R_s=1/2 beyond 51/100. It stays in [1/2,51/100], agrees
with r to all orders at 1/2, and is flatly constant after 51/100.
Let Δ_s(r,θ,φ) be tr(L_mid* L_c) at radius R_s(r), using the
extended (12.24). Define Δ_n by the same retreat in r_c=1/r,
using the north matrices. Let w_r=smooth_step((r-3/4)/(1/2)).
For 1/2≤r≤2 the global middle map is L_mid multiplied by

  Δ_s(r,θ,φ),                         r≤3/4,
  δ_s(θ,φ) exp[i w_r Θ(θ,φ)],          3/4<r<5/4,
  Δ_n(r,θ,φ),                         r≥5/4.             (12.35)

Both retreats are constant on neighbourhoods of the two interior
joins, and (12.34) gives exp(iΘ)=δ_n/δ_s. Thus (12.35) matches
to every order. At r=1/2 and 2 its jets equal the corresponding
cusp lift because L_mid tr(L_mid*L_c)=L_c and the retreat matches
the actual radius to all orders. Use the capped cusp maps of
§§12.4–12.5 beyond these boundaries. This supplies a global smooth
unit circle isomorphism covering the unchanged explicit orbit map;
its circle inverse is its adjoint. The cap inverse and the orbit
inverse remain the explicit ones in §11. The formula is a new
construction for review; it makes no assertion that (12.23) holds.

**Executed full-count phase check (d).** `boundary_phase64` evaluates
(12.34) in vectorized float64 arithmetic. Twenty-four calls, each
with the certified 65536 terms, took about 0.83 seconds in this
environment. Eight endpoint exponential checks have residual
<9.1·10⁻¹⁵; eight opposite-edge pairs differ by <4.3·10⁻¹⁵.
The smallest measured adjacent real part is 0.9999999944859101,
an excess >0.03124999 over the proved 31/32 bound. These are
observations about phase values, not a tensor or derivative benchmark.
Other term counts are labelled uncertified experiments by the helper.

The regional algebra and the phase proof close the radial gluing design.
The next subsection records the assembled value prototype and its first
tensor evaluation, separately from the remaining symbolic implementation
certificate.

### 12.10 Global values and endpoint phase derivatives

`direct_global_coordinates65` and `direct_global_forward65` assemble
(12.35), the existing orbit and torus maps, and the existing cap pair.
The numerical prototype evaluates the finite phase sum with NumPy;
all other blocks use mpmath. In the collars it uses the cancelled
(12.24) directly: a floating-point evaluation of the old expression
`1-smooth_step` can round the source slot weight to zero immediately
below r_c=1/2, before its algebraic cancellation. This changes no exact
formula. It prevents that artificial zero denominator in the prototype.

C65 checks 65 complete G-then-F values, with two nodes, eight nonnode
central-fibre points, ten further pole-band points, and forty-five
points at fifteen middle/cusp radii. All three holomorphic chart kinds
occur. The maximum round-coordinate/torus residual is <5.7·10⁻⁶⁰
at 60 digits (d). Because forward and inverse use the same numerical
phase evaluator, this is an inverse-consistency observation, not a
60-digit accuracy claim for the float64 phase.

The first and second derivatives of the finite phase sum need no loop.
Writing subscripts for derivatives in θ,φ, they are exactly

  Θ_a = Im(g_a/g),
  Θ_ab = Im(g_ab/g-g_a g_b/g²).                           (12.36)

C66a differentiates atan(v/u) symbolically and verifies these formulas
in a right-half-plane chart. C66b checks that the first and second
logarithmic derivatives of each quotient g_j/g_{j-1} are the
differences of those for g_j and g_{j-1}. Summation telescopes;
g₀ is constant. The formulas consequently apply to the specific
65536-term (12.34), including across its periodic seams (a).

This gives a practical evaluation strategy for the actual tensor.
In a round coordinate chart q, let h(q) be the holomorphic Y-chart
coordinates returned by the assembled G. Put A=Dh, A_k=∂_k A and
let J₀ be the three standard real 2×2 complex blocks. Then

  J=A⁻¹J₀A,       ∂_kJ=A⁻¹(J₀A_k-A_kJ).                (12.37)

The established generic identities C10/C30 apply to these expressions
when A is invertible and A_k are the symmetric second jets of h.
The global circle inverse, unchanged orbit inverse, torus inverse
and cap inverse give invertibility for the mathematical construction.
Binding all of the new implementation's derivatives to those formulas
is still necessary for the requested direct executable certificate.

`direct_tensor66` is a **numerical local evaluator**, not that completed
certificate. It evaluates (12.36) at the endpoint with high-precision
numerical derivatives, then substitutes a quadratic phase germ having
those same two jets into the surrounding formula. Other derivatives
use `mp.diff` in one fixed holomorphic chart. This avoids differentiating
the float64 phase routine, which would incorrectly erase infinitesimal
variations. The germ is only a device for evaluating the jets at its
centre; it is not a replacement global map.

At the strict interpolation point r=1.07, θ=.81, φ=-.63, Hopf
phase .37 and torus coordinates (.17,-.23), C66c evaluates the
specific J, its first derivative, and all 216 Nijenhuis components.
The round chart omits coordinate 0 with positive sign. The measured
square residual is <5.6·10⁻⁶⁰ and the maximum N component <5.3·10⁻⁵⁸.
Computing the value, first and second coordinate jets took about
0.58 seconds here. Two independent centred differences of the actual
J evaluator, at step 10⁻⁵, agree with ∂J to relative error <3.5·10⁻⁷
(d). A separate half-step comparison gives <8.8·10⁻⁸, consistent
with the expected second-order finite-difference error. This is one
interior tensor point, not global numerical coverage or a symbolic
check of every new branch.

**Next acceptance gate.** Implement the finite first/second jet circuit
without numerical differentiation, bind it to the new global source
formulas and every selector/overlap, and check the resulting specific
J in the cusps, caps, central fibres and radial/argument seams. The
large phase count is now a practical value computation; its derivatives
are the fixed-size endpoint expressions (12.36). The legacy runtime
and its certificates remain unchanged.

### 12.11 Finite tensor circuit and the resulting certificate

The next gate of §12.10 is addressed by C67–C68. `direct_tensor67`
uses a finite second-jet algebra instead of numerical differentiation.
For each scalar it stores its value, six first derivatives and a 6×6
Hessian. Products use the two cross terms, and a unary function uses

  D f(v)=f′(v)Dv,
  D²f(v)=f′(v)D²v+f″(v)Dv⊗Dv.                          (12.38)

C35 already verifies these generic composition rules. C67a checks the
first and second derivatives of every analytic primitive in the new
backend: powers, exp, expm1, log, log1p, sin, cos, sinh, cosh and asinh.
C66 supplies Arg and the endpoint phase rule. C67b/C68b match the
actual product and unary implementations to those checked expressions.
Conjugation and real/imaginary parts act componentwise on real jets.
Small matrix operations are the corresponding finite sums and products.
The final 6×6 solves in (12.37) are rational operations, equivalently
`adj(A)/det(A)`; no nonlinear equation or iterative inverse is used.

The value functions are executed in an isolated namespace with this
arithmetic. The source adapter changes the mpmath binding and rewrites
`abs(z)**(2k)` as `(Re(z)²+Im(z)²)**k` for positive integers k.
C68c checks that transformation and preservation of the other operations.
This rewrite is essential at a central axis: differentiating abs and
then squaring would introduce an artificial singularity. A genuinely
singular intermediate jet is marked invalid; it may be discarded by an
inactive flat branch, but a nonfinite output jet raises an error.
No live module is monkey-patched and the frozen runtime source is unchanged.

The inherited smooth base, source embedding, torus and cap blocks retain
the C38–C50 certificates. C68 pins the twenty new value/jet definitions
by their complete AST digests. It also enumerates **all 56 if/if-expression
tests** in the new geometric value path and **all seven `where` selectors**
in the vectorized boundary evaluator. Their classification is as follows;
these written arguments are inputs to the source-identity checks, not
consequences of a digest alone.

| selector family | justification |
|---|---|
| `infinity`, `inf`, cusp chart kind | A fixed cusp label on each region. The cusp/middle representation changes are the existing holomorphic changes, C38–C45 and C50. |
| Middle χ and cusp χ_c active-row tests | The weights are flat at 0 and 1. All active denominators are positive by §§12.2–12.4 and (12.30). The extra overlap computation only supplies diagnostics. |
| Pole correction present/absent and its flat weight | In the actual cusp path the correction is supplied; outside its support it is identically the uncorrected matrix. Equation (12.20) and §12.5 give positivity and exact cap agreement. |
| `phase_offset` and the two tests in `collar_offset61` | It is supplied precisely for r_c<1/2 in the extended collar. The south route has r<1/2; the north route has r>2. C61 proves exact phase and flat endpoint matching, so its exception is unreachable. |
| `rc>=1/2` | The linear height extends the original height with all jets, C45 and (12.24). |
| `H>=0` in the second cusp row | The row is inactive near H=0. At H=±1/2 its phase cancels the exact source deck phase, §12.7. |
| Boundary `th>=0` | The logarithm representative switches inside the first angular plateau; the phase-corrected middle formula has the exact angular transitions of C58/C61. |
| Vectorized step interval tests | On an inactive interval the result is exactly 0 or 1; its surrogate interior argument has no effect on the output. The endpoint joins are the flat-step joins of L1/C40/C45. |
| Vectorized active-row normalization | Its zero-norm guard cannot fire on an active row by the four projection bounds in §12.8. Inactive rows have zero weight. The jet circuit uses the equivalent matrix formula and never differentiates an inactive zero norm. |
| Phase count, angle and adjacent-ratio guards | The actual calls use exactly 65536, principal input angles, and Re ratio >31/32 from (12.31). C68e binds both call sites to that fixed default. |
| Global r=1/2,2,3/4,5/4 and retreat selectors | Equation (12.35) and flatness prove complete jet matching. The additional r_c=49/100 implementation selector switches two exactly equal forms of L_c, by source-weight cancellation (12.24). |
| `phase_local` callback | The value evaluator uses (12.34); the finite jet evaluator supplies its same value and its exact jets (12.36). It is not an arbitrary callback in the certified object. |
| Cap distance and forward source-radius tests | These are the existing §11.6 cap domains and inverse selectors. C60 makes the new regional map agree exactly with their raw lift on the required open overlap. The forward cap returns the existing explicit cap inverse directly. |
| Source circle Arg and middle/cusp return kind | The unit phase is nonzero by C56. Its argument changes x₁ by an integer or applies the existing cusp deck action. C30k and C38–C50 supply the holomorphic changes and torus restoration. The numerical unit-phase guard is identically false in exact arithmetic. |

This classification includes tests inside local helpers. The census finds
no while loop, try block or match statement in that geometric value path.
The inherited numerical precision shortcuts have the same status as in
C46/C48/C50: **they are off in the exact mathematical object**. Decimal
literals denote rational numbers; pi/e and the elementary functions have
their mathematical meanings; array storage, arithmetic precision casts
and `workdps` have no mathematical effect. Float64/multiprecision runs
are approximate evaluations of that object. This is not a claim of exact
floating-point arithmetic or a certified rounding-error bound.

**Certificate statement, for human review.** Assume the written elementary
arguments for the unchanged base, two-way heights, torus maps and caps
listed in §11.9, and the new inequalities and seam arguments in §§12.2–12.9.
Then (12.35) with the capped cusp maps defines a smooth explicit
diffeomorphism G_dir from round S⁴×T² to the same Y, with explicit inverse
F_dir assembled in C65. In each of the ten round sphere charts, together
with local torus angles, (12.37) defines a finite elementary tensor
J_dir=G_dir*J_Y. The domains cover because some round coordinate is nonzero;
the torus formulas are periodic by the inherited torus equivariance.
The maps agree through the holomorphic changes of Y, so this tensor is
globally well defined. Its coordinate first/second jets are obtained by
the finite circuit (12.36)–(12.38), bound to the specified source by
C67–C68 and the unchanged source certificates. Substitution in the exact
six-dimensional identity C30i gives all 36 components of J_dir²+I and
all 216 components of N_{J_dir} equal to zero. These conclusions use the
written inequalities and source correspondence, not numerical residuals.

There are no astronomical projection counts in this construction.
Only phase values use the fixed 65536-term sum; first and second phase
derivatives use the fixed-size endpoint expressions. The forward and
inverse maps use the same circle isomorphism and its adjoint. The original
§11 construction and its source remain a separate, unchanged reference.

**Executable observations (d).** The finite jet circuit agrees with the
independent C66 high-precision differentiation at its interpolation point
to relative errors <6·10⁻¹⁶ in J and <1.2·10⁻¹⁴ in dJ. A call takes about
1.4 seconds here, including source adaptation. C67c runs 31 more tensor
points: nine radial points, eight nonnode central-fibre points, and fourteen
cap/pole-band points including both nodes and four points in the actual
cap compression annuli. All three target chart kinds
occur. The largest observed scaled N residual is <1.2·10⁻¹⁵ and scaled
square residual <7.7·10⁻¹⁶; the computed Hessians are symmetric. The N
scale is `1+max|J| max|dJ|`; the square scale is `1+||J||²`. The maximum
sampled coordinate-Jacobian condition number is <223.

An initial binary64-only jet implementation failed near the height
transition H=±1/10: underflow in the conjugated flat functions erased
the recovered height derivative. The retained circuit uses arbitrary
precision/exponent arithmetic for its jets; all eight central-fibre
tests, including those four transition points, pass. The phase values
remain vectorized float64 evaluations, and the final small matrices are
converted to float64. The timing and residual figures are measurements,
not uniform runtime, conditioning or floating-point error guarantees.

Run the certificate with the existing command:

    .venv/bin/python3 -B extra/s4t2_explicit/sym/chart_identities.py

It must finish with `ALL PASS`, including C62–C68. For an evaluation
session, `runpy.run_path` on that same file returns a namespace containing
`direct_global_coordinates65`, `direct_global_forward65` and
`direct_tensor67`. The latter takes four retained round coordinates and
two torus coordinates, with `missing` in 0,…,4 and `sign` in {-1,1} for
the reconstructed sphere coordinate. It returns J, dJ, Dh, D²h, the
Nijenhuis components, and the selected holomorphic chart. The script's
fixture `input66, missing66, sign66` is a reproducible interior example.

The next work is review of this new global certificate, broader independent
implementation checks (especially near the nested flat transitions), and
simplification. The stronger single-argument inequality (12.23) is still
open; proving a suitable inequality at the exact boundary could remove
the 65536-term value sum altogether. No such simplification is assumed here.

### 12.12 A standalone function and its exact symbolic definition

The collaborator requested a callable whose subfunctions carry explicit
mathematical contracts, followed by executable square and Nijenhuis checks,
and authorized the two new files `extra/s4t2_explicit/direct.py` and
`extra/s4t2_explicit/sym/direct_certificate.py`. The implementation is
independent of the historical check script at import and evaluation time.
It implements the same construction as §12.11; it does not change the
global statement's written hypotheses or their human-review status.

**Public interface.** From the repository root:

```python
from math import sqrt
from extra.s4t2_explicit.direct import complex_structure

point = [1 / sqrt(5)] * 5 + [0.17, -0.23]
J = complex_structure(point)                 # a 6 by 6 numpy array
symbolic = complex_structure(point, symbolic=True)
J_exact = symbolic.matrix                   # a factored SymPy expression
```

The input is `(x0,x1,x2,x3,x4,t1,t2)`, with sum x_i²=1 and torus
coordinates of period one. The numerical interface checks the sphere
constraint to tolerance 10⁻¹² and reconstructs the largest sphere
coordinate in its hemisphere chart. Exact mathematics uses the sphere
constraint itself. The symbolic call keeps six chart coordinates free:
the supplied point is a witness selecting formula branches, not a
substitution turning the formula into six fixed constants. Branch choices
are available in `symbolic.graph.choices` and must be interpreted with
the domain and seam contracts below.

**One global matrix convention (a; D1).** For x in S⁴, use the smooth
identification

  Φ_x(v,a,b)=(v+a x,b),       x·v=0,
  Φ_x⁻¹(w,b)=(w-(x·w)x,x·w,b).                          (12.39)

Here a,b are tangent components in the period-one torus coordinates.
The public matrix is Φ_x J Φ_x⁻¹. If T is the matrix taking hemisphere
coordinate components to (12.39), the returned expression is
T A⁻¹J₀A T⁻¹. D1 checks both inverse identities for every omitted
sphere coordinate, using the relation sum x_i²=1. It also differentiates
the actual `frame_matrices` formula and checks `frame_derivatives`.
The hemisphere domains cover: otherwise all five sphere coordinates
would vanish, contradicting their squared norm one (b).

This frame is not a coordinate frame. Its vector fields are
E_i=(e_i-x_i x,x_i,0), i=0,…,4, and E_5=(0,0,1). Direct polynomial
differentiation in D1 gives

  [E_i,E_j]=x_i E_j-x_j E_i  (i,j<5),   [E_i,E_5]=0.     (12.40)

For coordinate calculations, `coordinate_tensor(q,missing=…,sign=…)`
returns J, its six derivatives, A, D²h and all Nijenhuis components.
`global_frame_data(x,coordinate_tensor(...))` also returns the public
matrix, its frame directional derivatives, and N with the bracket terms
included. Thus the public matrix is not accidentally tested using a
coordinate formula in a noncommuting frame. D9 compares that calculation
with the tensor transformation of the coordinate N (numerical, d).

**One definition, two arithmetic backends.** `Formulas` contains 39
geometric methods. The finite `Jet` arithmetic propagates first and
second derivatives by (12.38); its numerical shortcut subclass is described
below. The `SymbolicArithmetic` backend executes
the same methods using pairs of exact real SymPy expressions. Constants
are rational numbers, π and e; even powers of a complex modulus are
polynomials in its real and imaginary parts. `ExpressionGraph` retains
named assignments rather than expanding the same large subexpression
many times. To inspect the definition and its derivatives:

```python
values = symbolic.graph.active_assignments(symbolic.outputs)
A, Hessians, derivative_assignments = symbolic.differentiate()
J_coordinate, dJ_coordinate, _ = symbolic.tensor_expressions()
phase_definition = symbolic.graph.finite_phase_definition()
```

Every assignment's first and second partials are computed by SymPy
differentiation; the chain rule then assembles the six input derivatives.
The six Hessians are exact named expressions, not numerical differences.
D4/D10 check the actual primitive rules, complex arithmetic and mixed
chain rule (a). D10 additionally compares the actual named derivative
compiler with direct differentiation of an expanded six-variable
expression, and deliberately corrupts a second derivative in an isolated
copy: the checker rejects it. D5 builds actual symbolic expressions in
both poles, a cusp interior, the phase interpolation, and both collars.
Those six traces contain no floating constants or unevaluated derivatives.
They are implementation checks on those branches; they do not constitute
an exhaustive proof of global coverage.

The phase value is the specified 65536-term elementary sum (12.34),
available as an unexpanded finite `Sum`; the function called
`boundary_mismatch` in that sum is the concrete `Formulas` method.
Its local differential is that of an argument of its nonzero endpoint,
as in (12.36). D3 differentiates this endpoint expression twice and
checks the quotient cancellation underlying telescoping (a).
The replacement is a rule for derivatives, not a replacement of the
global phase value by principal Arg. Its validity uses the written
nonvanishing, increment and winding proofs of §§12.8–12.9 (b,c).
The scalar assignments and derivative assignments, together with this
finite-sum definition, specify the factored matrix completely. Matrix
inverse in this expression means the finite rational operation adj(A)/det(A).

**Exact tensor check without choosing A=I (a; D2).** Write
H_{aij}=∂_i∂_j h_a. Substitution of (12.37) in the coordinate N formula,
followed by multiplication by A, gives

  (AN)_{aij}=H_{aji}-H_{aij}
    + Σ_{b,l}(J₀)_{ab}(J_{li}H_{bjl}-J_{lj}H_{bil}
                         +H_{bli}J_{lj}-H_{blj}J_{li})
    + Σ_{l,m}(-J_{li}H_{aml}J_{mj}+J_{lj}H_{aml}J_{mi}).  (12.41)

To obtain the last two terms of the first sum, use AJ=J₀A on the
outer J in the last two terms of N; use J₀²=-I on the resulting
H_{aji}-H_{aij} term. Symmetry H_{aij}=H_{aji} cancels the first
line and the first sum. Interchanging l,m cancels the last sum.
D2 expands and reduces all 216 polynomials to zero with independent
symbolic entries of J and symmetric H. The same check reduces
(A⁻¹J₀A)² to -I, and executes the actual symbolic tensor constructor
on formal A,H to bind its returned J and dJ to these formulas.
It also executes the actual N routine on a symbolic
symmetric-H example and a nonintegrable negative control. The hypotheses
are that H is the Hessian of the implemented h and det A≠0; the named
derivative program supplies the former, and the explicit inverse
contracts supply the latter. Numerical residuals are not premises here.

**What a reviewer must still read.** `proof_contracts()` associates every
one of the 39 geometric methods with exactly one of the following ten
groups. D6 verifies that coverage and compares every complete method body
with its earlier source, using a displayed, limited AST rewrite: specialize
the elementary construction, remove numerical shortcuts and unreachable
legacy alternatives, rename dependencies, and replace even modulus powers.
This establishes correspondence; neither the comparison nor the contract
list proves inequalities or smoothness.

| contract group | written facts to check | earlier checks and arguments |
|---|---|---|
| Flat step | Range, constant arms, flat endpoint joins. | L1; §§11.2, 11.9; C40/C45/C47. |
| Two-way height | Positive arguments, monotonicity, explicit inverse, smooth central axes and all-order joins. | §11.2; C29/C43–C45/C49–C50. |
| Pole inverse pair | Domains of the raw map, identity extension of the radial conjugation, both inverse identities and cap overlap. | §§11.3, 11.6; C25–C28/C46; N10. |
| Orbit coordinates | Covering by the dominant-coordinate choices, their agreement and invertibility. | §§11.2, 11.6, 11.9; C38–C45/C49–C50; N11. |
| Source line | Equivariance, chart changes, nonzero reference and normalization. | §11.1; C21–C22/C38–C45/C50. |
| Torus restoration | Deck shifts, periodicity and the holomorphic changes of target chart. | §§10.2, 11.6, 11.9; C30k/C38/C41–C44/C50g. |
| Regional lifts | Active rows nonzero, positive overlap bounds, normalized circle isomorphisms and adjoint inverses. | §§12.1–12.5; C56–C60. |
| Boundary extension | Cancellation of vanishing source weights, angular transitions and all-order height joins. | §§12.6–12.7; C45/C61–C62. |
| Periodic finite phase | Derivative bounds, zero winding, positive increments and periodic smooth lift. | §§12.8–12.9; C63–C64/C66. |
| Global assembly | Flat radial retreats, cap agreement and reversal of orbit, circle and torus maps on the covering. | §§12.5, 12.9, 12.11; (12.20), (12.35). |

These are the same written identity/inequality arguments used by §12.11,
not new assumptions introduced by packaging. With them, each local h is
smooth with invertible derivative, and target-chart changes are holomorphic.
Composition and the chain rule then identify the local matrices on overlaps;
(12.39) expresses their common tensor in one smooth global frame. Individual
chart selectors need not be continuous as labels. Their mathematical
outputs agree through the stated changes of representation. Once those
contracts are checked, (12.41) is the executable exact algebraic end of
the proof. The focused script alone is not a standalone proof of all ten
contracts, and a proposed global certificate still requires human review.

**Numerical evaluation and observed costs (d).** `NumericalFormulas` adds
the six inherited flat-tail shortcuts explicitly, separately from exact
`Formulas`. Literal nested exponentials can otherwise become impractical:
an exploratory evaluation of the height radial map at 1/5-10⁻⁸ exceeded
a three-second guard, and a terminal-tail evaluation encountered an
mpmath exponent-conversion overflow. The numerical shortcut path passes
the affected regional checks. D11 compares a less extreme height-tail
value with a separate 180-digit evaluation. The value cutoffs described
in §11.7 and `cusp.py` do not provide a uniform error bound for J or dJ;
no such bound is claimed. Exact symbolic tracing uses the full formulas
and can likewise be costly at extreme witness points.

The focused check ran in about 49 seconds here. Its six actual symbolic
graphs have 123–1521 active value assignments and 1722–20773 derivative
assignments. Independently evaluating the named derivative expressions
agrees with the jet backend after rounding to float64 at all six witnesses;
this is not a 60-digit accuracy claim. D8 checks 32 regional tensor/inverse
pairs, including both nodes, eight nonnode central-fibre points and the
pole compression annuli. The largest scaled N residual is <1.2·10⁻¹⁵,
scaled square residual <6.5·10⁻¹⁶ and inverse residual <4·10⁻⁵⁴.
The inverse comparison shares the same float64 phase in both directions.
D9 compares all five sphere charts at three points and integer torus
changes: relative differences are <3.9·10⁻¹⁵ and <5.2·10⁻¹⁵.
The public example above takes about 0.36 seconds for one call here.
These are measurements, not uniform conditioning, runtime or rounding
guarantees. Phase values and the final matrices are float64; intermediate
jets use mpmath with 60 decimal digits by default.

Run the focused verifier from the repository root:

```sh
.venv/bin/python3 -B extra/s4t2_explicit/sym/direct_certificate.py
```

It must print `ALL PASS`, with exact checks and numerical observations
identified in its output. The earlier C1–C68 script remains the broader
record of the construction's identities. The next review task is the ten
written contracts and their source correspondence; a later numerical
analysis task is a certified error bound for flat-tail derivative shortcuts.
Neither task is replaced by observing small square or Nijenhuis residuals.
