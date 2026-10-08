# The closed-form construction of a complex structure on S⁴ × T² (full write-up)

> **How to read this file.**
> * This is a Markdown version of the construction and proofs of the manuscript
>   *An Explicit Complex Structure on S⁴ × T²* (E. Hirst, T. Schettini Gherardini, T. S. R. Silva,
>   J. Stein). It follows `manuscript_draft/main.tex` at commit `74362a6`, the last commit that changed
>   `main.tex`. **The manuscript is authoritative.** If the two ever differ, the manuscript holds, and
>   this file should be corrected through a ledger entry.
> * Numbers follow the compiled manuscript. For example, Lemma 2.4 is `lem:toric` and Proposition 4.1 is
>   `prop:filling`. Sections 2–12 match the manuscript's sections; §0 and §13 are additions, and §1 has a
>   shortened introduction. Citation keys such as [AC26] or [Matsumoto85] are those of
>   `manuscript_draft/references.bib`.
> * The identities listed in the manuscript's Appendix A are checked in exact arithmetic by
>   `manuscript_draft/code/verify_construction.py`. **[check k]** refers to check k of that script, and the
>   full map is in Appendix A below. The other computations (for example the inverse of Θ, the ℂ²-action at
>   ∞, or the Pontryagin–Chern relations of Lemma 11.11) are checked by hand in the manuscript, not by the
>   script.
> * The two scripts live in `manuscript_draft/code/`. `extra/s4t2_analytic/` holds their environment, a
>   runner and the reference outputs.
> * The proof that Y ≅ S⁴ × T² is not constructive: it uses a classification argument (W is a homotopy
>   sphere, then the h- and s-cobordism theorems, or Iwase's theorem). An explicit diffeomorphism is the
>   goal of the symbolic campaign (`SYMBOLIC_CAMPAIGN.md`, `docs/s4t2_explicit_diffeo.md`).

---

## 0. Results at a glance

**The construction in one line.** For t ∈ ℂ*,

    X_t = (ℂ*_ζ × ℂ_w) / ⟨ g₂, g₃, g₄ ⟩,
    g₂(ζ, w) = (e^t ζ,     w + 1 − t²),
    g₃(ζ, w) = (e^{−it} ζ, w + i(1 + t²)),
    g₄(ζ, w) = (t ζ,       w + 2t).

This is a holomorphic family of complex 2-tori. Its two cusps, at t = 0 and t = ∞, have opposite
chirality. Fill both cusps with Mumford–Tate models, then perform one Kodaira logarithmic transform of
multiplicity one at t = 0. The result is a compact complex threefold Y.

**Main Theorem.** There is a compact complex manifold Y of complex dimension 3 that is diffeomorphic to
S⁴ × S¹ × S¹. In particular, the smooth manifold S⁴ × T² admits an integrable almost complex structure.

This answers Calabi's question [Calabi58] for V = T². Calabi asked whether V² × S⁴ is complex for every
closed orientable surface V. The case V = S² was settled in 2026 [GFL26, LLZ26], building on the complex
S⁶ of [AC26, Engel26]. **The Main Theorem uses none of the 2026 works.** It uses only the classical facts
(F1)–(F4) of §2.3 and explicit computations.

**Theorem B.** Y has the following properties.
1. There is a proper surjective holomorphic map f: Y → ℙ¹.
   - Its fibres over ℂ* are complex 2-tori X_t.
   - Its fibres over 0 and ∞ are I₁ × E₀ and I₁ × E_∞. Here I₁ is a nodal rational curve and
     E₀ ≅ E_∞ ≅ ℂ/(ℤ + iℤ).
   - The two cusps differ in the orientation of the marked lattice bases: (1, i) at 0 and (−1, i) at ∞.
2. Torus actions.
   - T² acts freely and holomorphically on Y, preserving f.
   - W = Y/T² is a smooth homotopy 4-sphere, and Y ≅ W × T² as principal T²-bundles.
   - ℂ² acts holomorphically by fibrewise translations, transitively on the smooth fibres.
3. c₁(Y) = 0 in H²(Y; ℤ), and K_Y ≅ f*𝒪_{ℙ¹}(−4). So K_Y is topologically trivial but not torsion in Pic(Y).
4. The smooth fibres are null-homologous in H₄(Y; ℤ). Hence Y is neither Kähler nor balanced.
5. ℳ(Y) = f*ℂ(ℙ¹), so the algebraic dimension of Y is 1.

In addition, W is diffeomorphic to the standard S⁴ (Proposition 8.1, via Iwase's classification
[Iwase84]). This gives a second proof of the Main Theorem.

**Theorem C (conditional).** Assume input (I3) of §11.1: the complex S² × S⁴ of [LLZ26], with its line
bundle L₊. Then:
- for every q ≥ 0, S⁴ × S³ × S^{2q+1} carries a complex structure, as a holomorphic elliptic-curve bundle
  over X₊ × ℙ^q;
- S⁴ × S³ × S⁴ × S³ carries one as a holomorphic elliptic-curve bundle over X₊ × X₊.

**Corollary D (conditional).** Assume the inputs (I1)–(I3) of §11.1. Then:
- a product of spheres of dimension ≤ 8 is complex if and only if it is almost complex;
- the only such products with no almost complex structure are S⁴, S⁸ and S⁴ × S⁴.

Corollary D uses the Main Theorem for exactly two products: S⁴ × T², and S⁴ × T⁴ ≅ Y × E. It uses Theorem C
for S⁴ × S³ × S¹.

---

## 1. Idea of the construction and conventions

### 1.1 Idea of the construction (manuscript §1.2–1.3)

**Matsumoto's fibration of S⁴.** Matsumoto showed that S⁴ is the total space of a genus-one *achiral*
Lefschetz fibration over S² [Matsumoto82].
- It has two singular fibres, of types I₁⁺ and I₁⁻ (opposite signs), with the same vanishing cycle a
  [Matsumoto85, Thm. 4.2]; see also [GompfStipsicz99, Ex. 8.4.7].
- It has no section: a twisted gluing kills the dual cycle b.
- It cannot be holomorphic, because a holomorphic family of elliptic curves only has positive nodes.

**For complex 2-tori the situation is different.** Write the fibre near a cusp as
X_t = ℂ²/⟨e₁, e₂, e₃, e₄⟩, with vanishing cycle e₁ = (2π, 0).
- The intermediate quotient ℂ²/⟨e₁, e₂, e₃⟩ is an extension of the elliptic curve E_t = ℂ/⟨ē₂, ē₃⟩ by
  ℂ* = ℂ/2πℤ. Here ē₂, ē₃ are the second coordinates of e₂, e₃.
- X_t is the quotient of this extension by the infinite cyclic group generated by e₄.
- The central fibre is I₁ × E, with smooth locus ℂ* × E.
- The orientation of E relative to the complex orientation of the fibre is an extra sign, which
  holomorphy does not fix.

**The family of §3 has three properties:**
- the elliptic curve is ℂ/⟨1, i⟩ at t = 0 and ℂ/⟨−1, i⟩ at t = ∞, so the two cusps have opposite
  chirality;
- the full rank-4 lattice stays nondegenerate for every t ∈ ℂ*, thanks to a null conic in the period map;
- the vanishing cycle is the same at both ends.

**The steps of the construction:**
1. fill both cusps (§4);
2. perform one logarithmic transform of multiplicity one, with class the coinvariant [e₄] (§5);
3. use translations along e₂, e₃ to get a free holomorphic T²-action (§6).

The quotient W = Y/T² carries a genus-one achiral fibration f_W: W → ℙ¹. It has one node of each sign and
a common vanishing cycle, which are the singular-fibre data of Matsumoto's fibration. We do **not** identify
f_W with Matsumoto's map. Instead we prove directly that π₁(W) = 1 and χ(W) = 2 (§7).

**Why two cusps with monodromies T^{±1}.** Let Y → ℙ¹ be a 2-torus fibration with Y ≃ S⁴ × T². Assume that
global sections of R¹f_*ℚ are detected on a smooth fibre. This holds, for instance, when every singular fibre
is of the type I₁ × E used here; the manuscript does not check it for other degenerations. Then b₁ = 2 forces every local monodromy to be a
transvection with a common vanishing direction (Remark 10.4). So the order-3 and order-4 monodromies of the
(3,4,∞) family of [AC26, Engel26, GFL26] are excluded under this assumption.

### 1.2 Conventions (end of manuscript §1)

- Manifolds are smooth, Hausdorff and second countable.
- ℂ² ≅ ℝ⁴ via (z, w) ↦ (Re z, Im z, Re w, Im w). This is the complex orientation, and det_ℝ is the real
  4 × 4 determinant.
- ω(u, v) = Im(ū v) = Re u · Im v − Im u · Re v for u, v ∈ ℂ.
- An action of a discrete group Γ on a locally compact Hausdorff space X is *properly discontinuous* if
  {γ : γK ∩ K ≠ ∅} is finite for every compact K. A free, properly discontinuous action by biholomorphisms
  has a complex manifold as quotient, and the quotient map is a covering.
- D_ε = {|t| < ε} and D*_ε = D_ε ∖ {0}.

---

## 2. Preliminaries (manuscript §2)

### 2.1 Lattices in ℂ²

**Lemma 2.1 (determinant formula).** Let e₁ = (2π, 0) and e_j = (x_j, ē_j) ∈ ℂ² for j = 2, 3, 4. Then

    det_ℝ(e₁, e₂, e₃, e₄) = 2π · det [ Im x₂   Im x₃   Im x₄
                                      Re ē₂   Re ē₃   Re ē₄
                                      Im ē₂   Im ē₃   Im ē₄ ]
                          = 2π · Im X · (Re ē × Im ē),

where X = (x₂, x₃, x₄), ē = (ē₂, ē₃, ē₄), and × is the cross product in ℝ³.

*Proof.* The real coordinate vector of e₁ is (2π, 0, 0, 0). Expand along the first column. ∎

**Lemma 2.2 (local determinant identity)** [check 10]. Let α₂, α₃, ē₂, ē₃, ē₄, ℓ ∈ ℂ with ω(ē₂, ē₃) ≠ 0.
Let φ: ℂ → ℝ be the ℝ-linear functional with φ(ē_j) = Re α_j for j = 2, 3, and put P = Re ℓ − φ(ē₄). Then

    det_ℝ((2π, 0), (−iα₂, ē₂), (−iα₃, ē₃), (−iℓ, ē₄)) = −2π P ω(ē₂, ē₃).

*Proof.* By Lemma 2.1 the left side is 2π L(Im X), where L(v) = det[v; Re ē; Im ē] and
Im X = −(Re α₂, Re α₃, Re ℓ).
- Write φ(x + iy) = px + qy. Then u := (φ(ē₂), φ(ē₃), φ(ē₄)) = p Re ē + q Im ē, so L(u) = 0.
- Hence L(Im X) = −L(u + (0, 0, P)) = −P L(0, 0, 1).
- Finally L(0, 0, 1) = Re ē₂ Im ē₃ − Im ē₂ Re ē₃ = ω(ē₂, ē₃). ∎

**Lemma 2.3 (families of lattices).** Let B be a complex manifold, and let v₁, …, v₄: B → ℂ² be holomorphic
with det_ℝ(v₁(b), …, v₄(b)) ≠ 0 for all b. Then ℤ⁴, acting on ℂ² × B by (v, b) ↦ (v + Σ n_j v_j(b), b),
acts freely and properly discontinuously. The quotient is a complex manifold, and the projection to B is a
proper holomorphic submersion with fibres the complex tori ℂ²/⊕ ℤ v_j(b).

*Proof.* Let V(b) ∈ GL₄(ℝ) have columns v_j(b).
- On a compact L ⊂ B, ‖V(b)⁻¹‖ ≤ C_L, so |Σ n_j v_j(b)| ≥ |n|/C_L. So only finitely many n move a compact
  subset of ℂ² × L into itself. Freeness is clear.
- The preimage of a compact L is the image of the compact set
  {(Σ s_j v_j(b), b) : s ∈ [0,1]⁴, b ∈ L}, so the projection is proper. ∎

### 2.2 The toric surface 𝒳

**The fan.** Let N = ℤ², v_n = (n, 1) and σ_n = cone(v_n, v_{n+1}) for n ∈ ℤ.
- The cones σ_n and their faces form an infinite fan Σ. Its support is the open upper half-plane
  together with the origin.
- Any two cones meet in a common face. Every compact subset of the open upper half-plane meets only
  finitely many cones; this is the finiteness used in Mumford's construction [Mumford72].
- The fan is **not** locally finite at the origin, which lies in every cone.
- Each σ_n is unimodular, since det(v_n, v_{n+1}) = −1.

𝒳 = 𝒳_Σ is the smooth toric surface obtained by gluing the affine charts U_{σ_n} ≅ ℂ² along the open
subsets corresponding to common faces [Fulton93]. (ζ, t) are the coordinates on its dense torus (ℂ*)²,
given by the characters (1, 0) and (0, 1).

**Lemma 2.4.**
0. *Separation.* 𝒳 is a Hausdorff, second countable complex manifold of dimension 2.
1. *Charts* [check 8]. 𝒳 = ⋃_n U_n with U_n ≅ ℂ², with coordinates x_n = ζ⁻¹ t^{n+1} and y_n = ζ t^{−n}.
   - On U_n: t = x_n y_n and ζ = x_n^n y_n^{n+1}.
   - On U_n ∩ U_{n+1} = {x_n ≠ 0}: x_{n+1} = t x_n and y_{n+1} = 1/x_n.
2. *The map t.* The character t extends to a holomorphic function t: 𝒳 → ℂ. Every point with t ≠ 0 lies
   in the dense torus, so t⁻¹(t₀) ≅ ℂ*_ζ for t₀ ≠ 0.
3. *Central fibre.* t⁻¹(0) = ⋃_n D_n, where D_n ≅ ℙ¹ is the divisor of v_n.
   - D_n ∩ U_n = {x_n = 0} and D_{n+1} ∩ U_n = {y_n = 0}.
   - So t⁻¹(0) is an infinite chain of ℙ¹'s: D_n ∩ D_{n+1} = {P_n}, the origin of U_n, and
     D_n ∩ D_m = ∅ for |n − m| ≥ 2.
4. *Torus action.* λ ∈ ℂ* acts by λ·(ζ, t) = (λζ, t), i.e. λ·(x_n, y_n) = (λ⁻¹x_n, λy_n). It preserves t,
   each D_n and each P_n.
5. *Shift* [checks 8, 15]. The lattice automorphism A(a, b) = (a + b, b) maps σ_n to σ_{n+1}.
   - The induced automorphism ψ of 𝒳 is (ζ, t) ↦ (tζ, t) on the torus.
   - It maps U_n to U_{n+1}, with x_{n+1} ∘ ψ = x_n and y_{n+1} ∘ ψ = y_n.
   - It commutes with the ℂ*-action, preserves t, and maps D_n to D_{n+1} and P_n to P_{n+1}.
6. *Log form* [check 8]. On U_n, (dζ/ζ) ∧ dt = −dx_n ∧ dy_n.

*Proof.* Items (1)–(6) are standard toric computations [Fulton93].
- (1): the dual basis of (v_n, v_{n+1}) is (m₁^{(n)}, m₂^{(n)}) = ((−1, n+1), (1, −n)), and x_n = χ^{m₁^{(n)}},
  y_n = χ^{m₂^{(n)}}.
- (2): ⟨(0, 1), v_n⟩ = 1 ≥ 0 for all n.
- (5): with A = [[1, 1], [0, 1]] one has Aᵀ m_i^{(n+1)} = m_i^{(n)}. Since χ^m ∘ ψ = χ^{Aᵀm}, this gives
  x_{n+1} ∘ ψ = x_n and y_{n+1} ∘ ψ = y_n.

(0) Let m = n + k with k ≥ 1, and write t_j = x_j y_j on U_j. On the dense torus, (1) gives x_m = t^k x_n and
y_m = t^{−k} y_n. So the pairs (p, q) ∈ U_n × U_m identified in 𝒳 form

    R_{n,m} = { (p, q) : t_n(p) = t_m(q),  x_m(q) = t_n(p)^k x_n(p),  t_n(p)^{k−1} x_n(p) y_m(q) = 1 }.

- For k = 1 this is the graph of the gluing map (x_n, y_n) ↦ (t x_n, 1/x_n) on {x_n ≠ 0}.
- For k ≥ 2 the last equation forces t ≠ 0 and x_n ≠ 0. So R_{n,m} is the graph of the identification of the
  two dense tori, which is the overlap U_n ∩ U_m (the cones meet only in {0}).

Each R_{n,m} is defined by polynomial equations, hence closed. By symmetry the same holds for m < n, and
R_{n,n} is the diagonal. So the equivalence relation on ⊔ U_n is closed. The quotient map is open, because
the identifications are homeomorphisms between open subsets. A quotient of a Hausdorff space by an open map
with closed equivalence relation is Hausdorff. The atlas is countable, so 𝒳 is second countable. ∎

**Slabs.** For 0 < ε < 1 put 𝒳^ε = t⁻¹(D_ε). For ρ ≥ 1 put

    W_n(ρ) = { p ∈ U_n : |x_n(p)| < ρ, |y_n(p)| < ρ, |t(p)| < ε }.

**Lemma 2.5 (slabs).** Let 0 < ε < 1/16.
1. The slabs W_n(2) cover 𝒳^ε.
2. If |n − m| ≥ 2, then W_n(2) ∩ W_m(2) = ∅.
3. Each W_n(2) is contractible. Each W_n(2) ∩ W_{n+1}(2) is connected and homotopy equivalent to S¹.
4. ψ(W_n(ρ)) = W_{n+1}(ρ), and λ·W_n(ρ) ⊂ W_n(ρ max(|λ|, |λ|⁻¹)).
5. If p ∈ W_n(2) and t(p) ≠ 0, then log|ζ(p)| / log|t(p)| ∈ (n − ½, n + 3/2).

*Proof.*
- (5) On the torus, |x_n| < 2 means |ζ| > |t|^{n+1}/2, and |y_n| < 2 means |ζ| < 2|t|^n. Take logarithms,
  divide by log|t| < 0, and use log 2/|log|t|| < ½.
- (1) Let p be a torus point and n = ⌊log|ζ|/log|t|⌋. Then |t|^{n+1} < |ζ| ≤ |t|^n, so |x_n| < 1 and
  |y_n| ≤ 1. Now let p ∈ D_n, so p ∈ U_{n−1} ∪ U_n.
  - If p ∈ U_n and |y_n(p)| < 2, then p ∈ W_n(2).
  - Otherwise either p = P_{n−1}, or p ∈ U_n with |y_n(p)| ≥ 2. In the second case p ∈ U_{n−1} ∩ U_n and
    |x_{n−1}(p)| = 1/|y_n(p)| ≤ ½. In both cases p ∈ U_{n−1} with y_{n−1}(p) = 0 and |x_{n−1}(p)| ≤ ½, so
    p ∈ W_{n−1}(2).
- (2) By Lemma 2.4(3), W_n(2) ∩ t⁻¹(0) ⊂ D_n ∪ D_{n+1}, and P_m ∈ W_n only for m = n. So the central parts
  are disjoint for |n − m| ≥ 2. For torus points with m ≥ n + 2, a common point would need
  |t|^{n+1}/2 < |ζ| < 2|t|^m, i.e. |t|^{m−n−1} > 1/4. This is impossible since |t| < 1/16.
- (3) In the coordinates (x, y) = (x_n, y_n), W_n(2) = {|x|, |y| < 2, |xy| < ε} is star-shaped about 0. By
  Lemma 2.4(1), W_n(2) ∩ W_{n+1}(2) = {½ < |x| < 2, |y| < min(2, ε/|x|)}, which retracts onto
  {|x| = 1, y = 0}.
- (4) This follows from Lemma 2.4(4),(5). ∎

**Lemma 2.6 (topology of 𝒳^ε).** For 0 < ε < 1/16, 𝒳^ε is connected and simply connected, and

    H_q(𝒳^ε; ℤ) ≅ ℤ, 0, ℤ[τ, τ⁻¹], 0   for q = 0, 1, 2, ≥ 3,

where τ acts on H₂ through ψ_*. Moreover, any homeomorphism p ↦ κ(t(p))·ψ(p), with κ: D_ε → ℂ* admitting a
continuous logarithm, acts on homology as ψ_*.

*Proof.* Put Z_N = ⋃_{|n|≤N} W_n(2). By Lemma 2.5, adding W_{N+1}(2) (or W_{−N−1}(2)) to Z_N attaches a
contractible set along a connected set ≃ S¹.
- Van Kampen shows that each Z_N is simply connected.
- Mayer–Vietoris shows that each such step adds one free generator to H₂ and changes nothing else. When
  W_{n+1} is attached to Z_n along W_n ∩ W_{n+1}, the boundary map H₂(Z_n ∪ W_{n+1}) → H₁(W_n ∩ W_{n+1}) ≅ ℤ
  is onto, and H₂(Z_n) is its kernel. (The space is Z_n ∪ W_{n+1}, not Z_{n+1}; ledger entry
  2026-10-01-TSG-3.)
- D_{n+1} lies in W_n ∪ W_{n+1}: on D_{n+1}, y_{n+1} = 1/x_n, so every point has |x_n| < 2 or |y_{n+1}| < 2.
  The circle {|x_n| = 1} ⊂ D_{n+1} cuts it into two discs, one in each slab. So the boundary map sends
  [D_{n+1}] to a generator of H₁(W_n ∩ W_{n+1}).

So H₂(Z_N) is free on the classes [D_m] of the compact components in Z_N. The direct limit H₂(𝒳^ε) is free on
{[D_n]}_{n∈ℤ}.
- ψ is holomorphic and maps D_n to D_{n+1}, so it preserves complex orientations and ψ_*[D_n] = [D_{n+1}].
- Hence H₂ ≅ ℤ[τ^{±1}] as a module, with [D_n] ↔ τⁿ.
- p ↦ κ(t(p))^u · p, with κ^u = exp(u log κ) and u ∈ [0, 1], is an isotopy from the identity to
  p ↦ κ(t(p)) · p. This gives the last statement. ∎

### 2.3 Standard facts used

- **(F1)** Principal T^k-bundles over a CW complex B are classified by H²(B; ℤ^k). Topologically trivial
  smooth principal bundles are smoothly trivial.
- **(F2)** (Milnor [Milnor68]) For an infinite cyclic covering X̃ → X with deck generator τ there is an exact
  sequence
  ⋯ → H_q(X̃) →(τ_* − 1)→ H_q(X̃) → H_q(X) → H_{q−1}(X̃) → ⋯.
- **(F3)** (Wall [Wall64]) Two smooth closed simply connected 4-manifolds with isomorphic intersection forms
  are h-cobordant. For homotopy spheres this is Θ₄ = 0 [KervaireMilnor63], where Θ₄ is the group of
  h-cobordism classes of homotopy 4-spheres. **It does not assert the smooth 4-dimensional Poincaré
  conjecture.**
- **(F4)** Wh(ℤ²) = 0 [BHS64]. The s-cobordism theorem [Kervaire65, Milnor66]: an h-cobordism of dimension
  ≥ 6 with vanishing Whitehead torsion is a product.

---

## 3. The family of 2-tori over ℂ* (manuscript §3)

**Construction 3.1.** On ℂ*_ζ × ℂ_w × ℂ*_t consider the commuting automorphisms

    g₂(ζ, w, t) = (e^t ζ,     w + 1 − t²,   t),
    g₃(ζ, w, t) = (e^{−it} ζ, w + i(1 + t²), t),
    g₄(ζ, w, t) = (t ζ,       w + 2t,       t).

Put Γ° = ⟨g₂, g₃, g₄⟩ and Y° = (ℂ* × ℂ × ℂ*)/Γ°.

**The lattice.** In the additive coordinate z, with ζ = e^{iz}, and for a choice of log t, the lifts of
g₂, g₃, g₄ together with z ↦ z + 2π are the translations by

    e₁ = (2π, 0),  e₂(t) = (−it, 1 − t²),  e₃(t) = (−t, i(1 + t²)),  e₄(t) = (−i log t, 2t).      (3.1)

Changing the branch of log t by 2πi replaces e₄ by e₄ + e₁. So Λ_t = ⊕_j ℤ e_j(t) does not depend on the
branch.

**Lemma 3.2** [check 1]. For all t ∈ ℂ*, with r = |t|,

    det_ℝ(e₁, e₂, e₃, e₄) = 2π (1 + r²) (2r² − (1 − r²) ln r) > 0.

*Proof.* Write t = a + ib. Then Im x₂ = −a, Im x₃ = −b, Im x₄ = −ln r, and

    Re ē = (1 − a² + b², −2ab, 2a),   Im ē = (−2ab, 1 + a² − b², 2b).

Then Re ē × Im ē = (1 + r²)(−2a, −2b, 1 − r²), so by Lemma 2.1 det_ℝ = 2π(1 + r²)(2a² + 2b² − (1 − r²) ln r).
The factors 1 − r² and ln r have opposite signs for r ≠ 1 and both vanish at r = 1, so (1 − r²) ln r ≤ 0.
Hence det_ℝ ≥ 4π r²(1 + r²) > 0. ∎

The determinant does not depend on arg t; check 1 verifies this too.

**Lemma 3.3** [check 2]. Γ° ≅ ℤ³ acts freely and properly discontinuously, and Y° is a complex manifold.
- f: Y° → ℂ*, [ζ, w, t] ↦ t, is a proper holomorphic submersion with fibres X_t = ℂ²/Λ_t, which are
  compact complex 2-tori.
- The monodromy of {H₁(X_t; ℤ) = Λ_t} around t = 0 (counterclockwise) is
  T: e₁ ↦ e₁, e₂ ↦ e₂, e₃ ↦ e₃, e₄ ↦ e₄ + e₁.
- The zero section t ↦ [1, 0, t] is a holomorphic section of f.

*Proof.* Let B ⊂ ℂ* be a simply connected open disc, with log t fixed on B. The map
(z, w, t) ↦ (e^{iz}, w, t) from ℂ² × B is the quotient by translation by e₁. Under it, Γ° corresponds to the
translations by e₂, e₃, e₄ modulo e₁.
- By Lemmas 2.3 and 3.2, the ℤ⁴-action on ℂ² × B is free and properly discontinuous. The quotient
  f⁻¹(B) → B is a proper holomorphic submersion with fibres ℂ²/Λ_t.
- Γ° preserves t. So freeness and proper discontinuity on all of ℂ* × ℂ × ℂ* follow by covering the image
  in ℂ* of a compact set by finitely many such discs.
- Γ° ≅ ℤ³ because e₁, …, e₄ are ℝ-independent.
- The monodromy statement is the branch change above. ∎

**Remark 3.4 (opposite chirality)** [check 9]. The *projected period vectors* ē₂(t) = 1 − t² and
ē₃(t) = i(1 + t²) satisfy ω(ē₂, ē₃) = 1 − |t|⁴.
- For |t| ≠ 1 they span a lattice in ℂ, with modulus τ_E(t) = ē₃/ē₂ = i(1 + t²)/(1 − t²). This is i at
  t = 0 and −i at t = ∞. Near the two cusps these are the lattices of E₀ and E_∞, which therefore have
  opposite orientations (Lemma 4.2).
- On |t| = 1 they are ℝ-linearly dependent and span no lattice.
  - For t = e^{iθ} with θ ∉ πℤ, τ_E = −cot θ is real. If cot θ is irrational, they span a dense subgroup of
    a real line.
  - At t = ±1, ē₂ = 0.
- By Lemma 3.2, the full rank-4 lattice Λ_t stays nondegenerate across |t| = 1.

The vector ē(t) = (1 − t², i(1 + t²), 2t) parametrises the null conic {x² + y² + z² = 0} ⊂ ℙ². This conic has
no real points, which is why Re ē × Im ē is nowhere zero.

### 3.1 Periods and Néron–Severi groups

Let B ⊂ ℂ* be an open disc with a branch of log t, so that e₁, …, e₄ is a basis of Λ_t = H₁(X_t; ℤ) depending
holomorphically on t.
- For j < k let Σ_jk ⊂ X_t be the 2-torus spanned by e_j, e_k, oriented by (e_j, e_k). The classes [Σ_jk]
  form a basis of H₂(X_t; ℤ) = Λ²Λ_t.
- The periods of σ = dz ∧ dw are p_jk(t) = ∫_{Σ_jk} σ = z_j w_k − z_k w_j, where e_j = (z_j, w_j).

**Lemma 3.5 (periods)** [check 12].
1. The periods are

       p₁₂ = 2π(1 − t²),   p₁₃ = 2πi(1 + t²),            p₁₄ = 4πt,
       p₂₃ = 2t,           p₂₄ = i(1 − t²) log t − 2it²,  p₃₄ = −(1 + t²) log t − 2t².

2. If c = (c_jk) ∈ ℤ⁶ and Σ_{j<k} c_jk p_jk vanishes identically on B, then c = 0.

*Proof.* (1) Parametrise Σ_jk by (a, b) ↦ a e_j + b e_k on [0, 1]². Then dz ∧ dw pulls back to
(z_j w_k − z_k w_j) da ∧ db. Insert (3.1).

(2) Write F = Σ c_jk p_jk = A(t) + C(t) log t with the polynomials

    C(t) = i c₂₄ (1 − t²) − c₃₄ (1 + t²),
    A(t) = 2π c₁₂ (1 − t²) + 2πi c₁₃ (1 + t²) + (4π c₁₄ + 2 c₂₃) t − 2(i c₂₄ + c₃₄) t².

- If C ≢ 0, then log t = −A/C on B. Continuing analytically once around 0 changes the left side by 2πi and
  leaves the right side unchanged, which is impossible. So C ≡ 0, and then A ≡ 0.
- The constant term of C is i c₂₄ − c₃₄. Since c₂₄ and c₃₄ are real, c₂₄ = c₃₄ = 0.
- The constant term of A is 2π(c₁₂ + i c₁₃), so c₁₂ = c₁₃ = 0.
- The linear term of A is 4π c₁₄ + 2 c₂₃. Since π is irrational, c₁₄ = c₂₃ = 0. ∎

**Proposition 3.6 (Néron–Severi groups of the fibres)** [check 12 for part 2; part 1 rests on Lemma 3.5(2)].
1. There is a countable subset ℰ of ℂ* such that, for every t ∈ ℂ* ∖ ℰ:
   - NS(X_t) = 0;
   - X_t contains no curves;
   - X_t has no non-constant meromorphic functions, i.e. a(X_t) = 0.
2. The fibre over t = 1 is a product of elliptic curves, X₁ ≅ ℂ/(2πℤ + iℤ) × ℂ/(2ℤ + 2iℤ), hence a projective
   abelian surface. In particular 1 ∈ ℰ, so not every smooth fibre is non-algebraic.

*Proof.* (1) By the Lefschetz (1,1)-theorem, NS(X_t) is the group of integral classes of type (1,1).
- A real class γ ∈ H²(X_t; ℝ) is of type (1,1) if and only if ∫ γ ∧ σ = 0. Only the (0,2)-part of γ
  contributes to this integral. That part is a multiple of σ̄, ∫ σ̄ ∧ σ ≠ 0, and γ^{2,0} is the conjugate of
  γ^{0,2}.
- Write the Poincaré dual of an integral class γ as Σ c_jk [Σ_jk] with c ∈ ℤ⁶, using the convention
  ∫_{X_t} α ∧ γ = ∫_{PD(γ)} α. Then ∫ γ ∧ σ = Σ c_jk p_jk(t).
- For c ≠ 0 this holomorphic function is not identically zero (Lemma 3.5), so its zeros in B are discrete,
  hence countable.
- ℂ* is covered by countably many discs, and ℤ⁶ is countable. So ℰ = {t : NS(X_t) ≠ 0} is countable.

Now let t ∉ ℰ. A curve C ⊂ X_t would have [C] ≠ 0 in NS(X_t), since ∫_C ω₀ > 0 for a flat Kähler form ω₀. A
non-constant meromorphic function has a non-empty polar divisor. So neither exists.

(2) At t = 1 the generators are e₁ = (2π, 0), e₂ = (−i, 0), e₃ = (−1, 2i) and e₄ = (0, 2).
- In the coordinates (z′, w) = (z − (i/2)w, w) they become (2π, 0), (−i, 0), (0, 2i) and (−i, 2).
- Replacing the last vector by e₄ − e₂ = (0, 2) exhibits Λ₁ as the product lattice
  (2πℤ + iℤ) × (2ℤ + 2iℤ). ∎

---

## 4. Filling the cusps (manuscript §4)

**Proposition 4.1 (Mumford–Tate filling).** Let ε₀ ∈ (0, 1/16). Let α₂, α₃, ē₂, ē₃, ē₄ be holomorphic on
D_{ε₀}, and let ν ∈ ℂ*, such that
- **(H1)** α₂(0) = α₃(0) = ē₄(0) = 0;
- **(H2)** ω(ē₂(t), ē₃(t)) ≠ 0 for t ∈ D_{ε₀}.

On 𝒳^{ε₀} × ℂ_w, with t = t(p), define

    G_j(p, w) = (e^{α_j(t)} · p, w + ē_j(t))   (j = 2, 3),     G₄(p, w) = (ν · ψ(p), w + ē₄(t)).

Then there is ε ∈ (0, ε₀) with the following properties, where Ũ = 𝒳^ε × ℂ and Γ = ⟨G₂, G₃, G₄⟩.
- (a) Γ ≅ ℤ³ acts on Ũ freely and properly discontinuously by biholomorphisms preserving t.
- (b) U = Ũ/Γ is a complex threefold, and f = t: U → D_ε is proper.
- (c) f⁻¹(0) ≅ I₁ × E, where E = ℂ/(ℤ ē₂(0) + ℤ ē₃(0)) and I₁ = t⁻¹(0)/⟨νψ⟩ is a nodal rational curve. Near
  every point of f⁻¹(0) there are local coordinates in which f = xy.
- (d) f⁻¹(D*_ε) is the open subset (ℂ*_ζ × ℂ_w × D*_ε)/Γ. On it G₂, G₃ act by
  (ζ, w) ↦ (e^{α_j} ζ, w + ē_j), and G₄ by (ζ, w) ↦ (ν t ζ, w + ē₄). Its fibres are the tori with lattice
  ⟨(2π, 0), (−iα₂, ē₂), (−iα₃, ē₃), (−i log(νt), ē₄)⟩.
- (e) The section t ↦ [(ζ = 1, t), 0] extends holomorphically over t = 0.

*Proof.*

*Step 0: constants.* Shrinking ε₀, we may assume the data are holomorphic near the closed disc D̄_{ε₀}. Then
there are c₀, C₀, C₁ > 0 such that for |t| ≤ ε₀, real a, b and u ∈ ℂ:
- |a ē₂(t) + b ē₃(t)| ≥ c₀(|a| + |b|);
- |α_j(t)| + |ē₄(t)| ≤ C₀|t|;
- |φ_t(u)| ≤ C₁|t||u|. Here φ_t: ℂ → ℝ is the ℝ-linear functional with φ_t(ē_j(t)) = Re α_j(t) for
  j = 2, 3. It exists by (H2), and the estimate follows from (H1).

For t ≠ 0 put P(t) = log|ν| + log|t| − φ_t(ē₄(t)). Then P(t)/log|t| → 1 as t → 0. Choose ε < ε₀ such that
P(t)/log|t| ∈ [½, 3/2] for 0 < |t| < ε. This ε depends only on the data, not on any compact set. For the two
cusps used below, P is computed exactly in Remark 4.3, and ε = 10⁻² suffices.

*Step 1: the action.*
- The maps are holomorphic and preserve t, since the ℂ*-action and t are holomorphic and t ∘ ψ = t.
- They commute, since ψ commutes with the ℂ*-action.
- On t⁻¹(0) × ℂ, G₂^a G₃^b G₄^k acts by (ν^k ψ^k, w + a ē₂(0) + b ē₃(0)). This is the identity only for
  a = b = k = 0, so Γ ≅ ℤ³.

*Step 2: the open part.* For t ≠ 0 the fibre of 𝒳 is ℂ*_ζ, and ψ acts by ζ ↦ tζ, which gives the formulas in
(d). By Lemma 2.2 with ℓ = log(νt), the lattice there has determinant −2π P(t) ω(ē₂, ē₃) ≠ 0.

*Step 3: a normalised height.* On torus points of Ũ put H(p, w) = log|ζ(p)| − φ_t(w). Since
log|e^{α_j}| = Re α_j = φ_t(ē_j),

    H ∘ G₂ = H,   H ∘ G₃ = H,   H ∘ G₄ = H + P(t).

*Step 4: freeness.* Suppose γ = G₂^a G₃^b G₄^k fixes (p, w).
- If t(p) = 0, then a ē₂(0) + b ē₃(0) = 0, so a = b = 0. And ν^k ψ^k(p) = p forces k = 0, because ν^k ψ^k
  maps D_n to D_{n+k} and P_n to P_{n+k}.
- If t(p) ≠ 0, a fixed point means a e₂ + b e₃ + k e₄ ∈ ℤ e₁ for the lattice of Step 2. This is impossible
  unless a = b = k = 0.

*Step 5: proper discontinuity.* Let K ⊂ Ũ be compact. Then K ⊂ (⋃_{|n|≤N} W_n(2)) × {|w| ≤ R} for some N, R
(Lemma 2.5(1)). Let γ = G₂^a G₃^b G₄^k and q = (p, w) ∈ K with γq ∈ K.
- *The w-coordinate.* Comparing w-coordinates gives |a ē₂ + b ē₃| ≤ 2R + |k| C₀ ε, so
  |a| + |b| ≤ c₀⁻¹(2R + |k| C₀ ε). It therefore suffices to bound |k|.
- *If t(p) = 0.* q and γq lie on components D_m with |m| ≤ N + 1, and γ shifts component indices by k. So
  |k| ≤ 2N + 1.
- *If t(p) ≠ 0.* By Lemma 2.5(5) and Step 0, |H(q′)/log|t|| ≤ N + 2 + C₁ R ε / log(1/ε) =: B_K for
  q′ ∈ {q, γq}. By Step 3, H(γq) − H(q) = k P(t), and P(t)/log|t| ≥ ½. Hence |k| ≤ 4 B_K.

*Step 6: properness of f.* Fix ε′ < ε and let Π_t = {s₂ ē₂(t) + s₃ ē₃(t) : s_j ∈ [0, 1]} and
R₀ = max_{|t|≤ε′} max_{u∈Π_t} |u|. Let q = (p, w) with |t(p)| ≤ ε′.
- If t ≠ 0:
  - replace q by G₄^k q with H ∈ (P(t), 0] (Step 3);
  - then replace it by G₂^a G₃^b q with w ∈ Π_t; this leaves H unchanged.

  Now log|ζ| = H + φ_t(w) ∈ (log|t| − C₂, C₂) with C₂ = |log|ν|| + C₀C₁ε² + C₁R₀ε. So in the chart U₀,
  |y₀| = |ζ| < e^{C₂} and |x₀| = |t|/|ζ| < e^{C₂}.
- If t = 0: use G₄^k to move p into D₀, then G₂, G₃ to put w ∈ Π₀.

So every Γ-orbit in f⁻¹(D̄_{ε′}) meets the compact set

    K₀ = { (p, w) : p ∈ ({|x₀|, |y₀| ≤ e^{C₂}} ∪ {|x₋₁|, |y₋₁| ≤ 1}),  |t(p)| ≤ ε′,  w ∈ Π_{t(p)} }.

This uses D₀ ⊂ {|x₀| = 0, |y₀| ≤ 1} ∪ {y₋₁ = 0, |x₋₁| ≤ 1}. So f⁻¹(D̄_{ε′}) is the image of K₀, hence compact.

*Step 7: central fibre and local form.* At t = 0, by (H1), G₂ = (id, w + ē₂(0)), G₃ = (id, w + ē₃(0)) and
G₄ = (νψ, w). So f⁻¹(0) = (t⁻¹(0)/⟨νψ⟩) × E.
- The quotient of the chain by νψ is D₀ with P₋₁ ~ P₀, i.e. I₁.
- In the chart U_n × ℂ, f = x_n y_n. This proves (c).

(b) follows from (a), Step 6 and the smoothness of Ũ.

(e) The point (ζ = 1, t) has x₀ = t and y₀ = 1. It tends to the point (0, 1) ∈ D₀ ∩ U₀. ∎

### 4.1 The cusp at t = 0

Take the data

    α₂ = t,  α₃ = −it,  ē₂ = 1 − t²,  ē₃ = i(1 + t²),  ē₄ = 2t,  ν = 1.

Then (H1) holds, and (H2) holds with ω(ē₂, ē₃) = 1 − |t|⁴. Over D*_ε the generators of Proposition 4.1(d) are
exactly g₂, g₃, g₄.
- U₀ → D_ε denotes the filling.
- ι₀: U₀|_{D*_ε} → Y°|_{D*_ε} is the identity of (ℂ* × ℂ × D*_ε)/Γ°.

### 4.2 The cusp at t = ∞

**Lemma 4.2** [checks 3, 9]. On ℂ* × ℂ × {|t| > 1} the map

    Θ(ζ, w, t) = (ζ e^{w/t}, w/t², 1/t) =: (ζ′, w′, s)

is a biholomorphism onto ℂ* × ℂ × D*₁, with inverse (ζ′, w′, s) ↦ (ζ′ e^{−w′/s}, w′/s², 1/s). It conjugates

    g₂      ↦ h₂ = (e^s ζ′,        w′ + s² − 1),
    g₃      ↦ h₃ = (e^{is} ζ′,     w′ + i(1 + s²)),
    g₄⁻¹    ↦ h₄ = (e^{−2} s ζ′,   w′ − 2s).

So Y°|_{|t|>1} ≅ (ℂ* × ℂ × D*₁)/⟨h₂, h₃, h₄⟩. The data

    α₂ = s,  α₃ = is,  ē₂ = s² − 1,  ē₃ = i(1 + s²),  ē₄ = −2s,  ν = e^{−2}

satisfy (H1) and (H2), with ω(ē₂, ē₃) = |s|⁴ − 1 < 0.

*Proof.* A map (ζ, w) ↦ (mζ, w + u), with m and u depending only on t, conjugates to
ζe^{w/t} ↦ m e^{u/t} · ζe^{w/t} and w/t² ↦ w/t² + u/t². Insert (m, u) = (e^t, 1 − t²), (e^{−it}, i(1 + t²))
and (t, 2t). The value of ω is a direct computation. ∎

- U_∞ → D_ε (s-coordinate) is the filling for these data, with ε shrunk so that both fillings exist.
- ι_∞: U_∞|_{D*_ε} → Y°|_{|t|>1/ε} is Θ⁻¹.
- The central fibre is I₁ × E_∞, where E_∞ = ℂ/(ℤ(−1) + ℤi) is ℂ/(ℤ + iℤ) with the marked basis (−1, i).
  This is oriented oppositely to the basis (1, i) at 0.
- The zero sections correspond, since Θ(1, 0, t) = (1, 0, 1/t).

**Remark 4.3 (the heights P exactly)** [check 11]. By Lemma 2.2, P = −det_ℝ/(2π ω(ē₂, ē₃)), where det_ℝ is the
determinant of the lattice in Proposition 4.1(d).
- **At t = 0** this lattice is Λ_t and ω = 1 − r⁴ with r = |t|. Lemma 3.2 gives

      P₀(r) = log r − 2r²/(1 − r²).

- **At ∞** put ρ = |s|.
  - In additive coordinates (z″, w′) with ζ′ = e^{iz″}, Θ is the ℂ-linear map (z, w) ↦ (z − iw/t, w/t²).
    Its real determinant is |t|⁻⁴ = ρ⁴.
  - It maps (e₁, e₂, e₃, −e₄) to a basis of the lattice at ∞, up to adding multiples of e₁.
  - Hence det_ℝ = −ρ⁴ · 2π(1 + r²)(2r² − (1 − r²) log r) with r = 1/ρ, which equals
    −2π(1 + ρ²)(2 − (1 − ρ²) log ρ). With ω = ρ⁴ − 1 this gives

        P_∞(ρ) = log ρ − 2/(1 − ρ²).

So, for 0 < r < 1,

    P₀(r)/log r = 1 + 2r²/((1 − r²)|log r|),     P_∞(r)/log r = 1 + 2/((1 − r²)|log r|).

Both are increasing in r and tend to 1 as r → 0.
- For 0 < r < 10⁻² they lie in (1, 1.000044) and (1, 1.434339). So ε = 10⁻² satisfies Step 0 at both cusps.
- At r = 1/16 the second ratio is 1.724… > 3/2, so ε = 1/16 does **not** work at ∞.

---

## 5. The logarithmic transform and the threefold Y (manuscript §5)

The following map is a logarithmic transform of multiplicity one in the sense of Kodaira [Kodaira63]: a
regluing by a holomorphic section over a punctured disc. Its topological type is the nonzero class [e₄] in the
coinvariants Λ/(T − 1)Λ.

**Lemma 5.1** [check 4]. On Y°|_{D*_ε} define

    Φ[ζ, w, t] = [ e^{iπ s(s−1)} ζ,  w + 2ts,  t ],     s = log t / (2πi).

1. Φ is a well-defined biholomorphism over the identity of D*_ε and commutes with g₂, g₃, g₄.
2. In additive coordinates Φ is the fibrewise translation by σ̃(s) = (π s(s − 1), 2ts), and
   σ̃(s + 1) − σ̃(s) = e₄(s), where e₄(s) = (2πs, 2t).

*Proof.* Replacing s by s + 1 multiplies e^{iπs(s−1)} by e^{iπ((s+1)s − s(s−1))} = e^{2πis} = t, and adds 2t to
w + 2ts. So the two values differ by g₄, and Φ is well defined on the quotient.
- Φ is locally given by a holomorphic formula.
- It commutes with the g_j, because all these maps have the form (ζ, w) ↦ (mζ, w + u) with m, u depending
  only on t.
- Its inverse is given by the same formula with the signs reversed.
- (2) is a direct computation, using −i log t = 2πs. ∎

**Construction 5.2.** Fix ε as above and put Y_mid = Y°|_{ε/2<|t|<2/ε}. Define

    Y = (U₀ ⊔ Y_mid ⊔ U_∞) / ∼

with the identifications
- q ∼ Φ(ι₀(q)) for q ∈ U₀ with ε/2 < |t(q)| < ε;
- q ∼ ι_∞(q) for q ∈ U_∞ with ε/2 < |s(q)| < ε.

The maps t and s = 1/t glue to f: Y → ℙ¹.

**Proposition 5.3.** Y is a compact connected complex manifold of dimension 3, and f is a proper surjective
holomorphic map. Over ℂ* its fibres are the tori X_t, and f⁻¹(0) ≅ I₁ × E₀ and f⁻¹(∞) ≅ I₁ × E_∞.

*Proof.* The gluing maps are biholomorphisms between open subsets, compatible with f and equal to the identity
on the base.
- *Hausdorff.* Each piece is the full preimage of an open subset of ℙ¹ and maps injectively onto an open subset
  of Y. Since the gluings cover the identity of the base, every fibre of f lies entirely in at least one piece.
  Let q₁ ≠ q₂.
  - If f(q₁) ≠ f(q₂), separate them by preimages of disjoint discs.
  - If f(q₁) = f(q₂) = t₀, choose a piece containing the whole fibre over t₀. It contains both points, and
    they can be separated inside it.
- *Compact.* Y is the union of the compact sets f⁻¹(|t| ≤ 3ε/4), f⁻¹(3ε/4 ≤ |t| ≤ 4/(3ε)) and
  f⁻¹(|s| ≤ 3ε/4), by Proposition 4.1(b) and Lemma 3.3.
- *Connected.* Each piece is connected and the pieces overlap. ∎

---

## 6. The torus action (manuscript §6)

**Proposition 6.1** [check 5]. For c = (c₂, c₃) ∈ ℝ² define, on each piece,

    on Y°:   A_c[ζ, w, t]   = [e^{(c₂ − ic₃)t} ζ,   w + c₂(1 − t²) + ic₃(1 + t²),   t],
    on Ũ₀:   Ã_c(p, w)      = (e^{(c₂ − ic₃)t} · p,  w + c₂ ē₂(t) + c₃ ē₃(t)),
    on Ũ_∞:  Ã′_c(p, w′)    = (e^{(c₂ + ic₃)s} · p,  w′ + c₂(s² − 1) + ic₃(1 + s²)).

These define an action of T² = ℝ²/ℤ² on Y by biholomorphisms preserving f. The action is free and
real-analytic. Allowing c ∈ ℂ² gives a holomorphic action of (ℂ*)² = ℂ²/ℤ².

Write m_c for the multiplier of the 𝒳-factor: m_c = e^{(c₂ − ic₃)t} on Ũ₀ and m_c = e^{(c₂ + ic₃)s} on Ũ_∞.
Then on both pieces Ã_c(p, w) = (m_c · p, w + c₂ ē₂ + c₃ ē₃), with the data ē₂, ē₃ of the respective cusp.

*Proof.*

*Well defined on each piece.*
- Each formula is additive in c. The values c = (1, 0) and (0, 1) give g₂, g₃ (resp. G₂, G₃ and h₂, h₃), so
  ℤ² acts trivially on the quotients.
- Each formula is a map (ζ, w) ↦ (mζ, w + u) with m, u depending only on the base. So it commutes with g₄,
  G₄, h₄ and with Φ.

*Compatible with the gluings.*
- On U₀|_{D*} the formula is that of Y°, since ι₀ is the identity.
- The conjugate of A_c by Θ is Ã′_c. This is the computation of Lemma 4.2 with
  (m, u) = (e^{(c₂ − ic₃)t}, c₂(1 − t²) + ic₃(1 + t²)):
  m e^{u/t} = e^{(c₂ + ic₃)/t} and u/t² = c₂(s² − 1) + ic₃(1 + s²).

*Free.*
- On X_t with t ∈ ℂ*, A_c is translation by c₂e₂(t) + c₃e₃(t). This lies in Λ_t only if c ∈ ℤ², by
  ℝ-independence of e₁, …, e₄.
- On f⁻¹(0) = I₁ × E₀ the action is (id, u ↦ u + c₂ + ic₃), since e^{(c₂ − ic₃)·0} = 1. This is free on T².
- On f⁻¹(∞) it is translation by −c₂ + ic₃ on E_∞.

*Holomorphic.* Each A_c is holomorphic in all variables for fixed c, also for complex c. ∎

**Corollary 6.2.** W := Y/T² is a smooth, closed, connected, orientable 4-manifold, and π: Y → W is a smooth
principal T²-bundle. f factors as f = f_W ∘ π with f_W: W → ℙ¹ smooth and proper.

*Proof.* A free smooth action of a compact Lie group gives a smooth principal bundle. The quotient of an
oriented manifold by a free action of a connected group is orientable. f is T²-invariant. ∎

---

## 7. Topology of W (manuscript §7)

**The cover of W.** Let N₀ = π(U₀), M = π(Y_mid) and N_∞ = π(U_∞).
- These are open in W and cover it, and N₀ ∩ N_∞ = ∅.
- O₀ = N₀ ∩ M lies over the annulus ε/2 < |t| < ε, and O_∞ = N_∞ ∩ M over 1/ε < |t| < 2/ε.

### 7.1 The regular part

Let e₁*, …, e₄* be the real dual basis of (e₁, …, e₄) on ℂ² = ℝ⁴.
- Under e₄ ↦ e₄ + e₁ the dual basis becomes (e₁* − e₄*, e₂*, e₃*, e₄*). So e₄* is a well-defined continuous
  function of t ∈ ℂ*.
- e₄*(Λ_t) ⊂ ℤ, and e₄*(e₂) = e₄*(e₃) = 0.

**Lemma 7.1.**
1. For t ∈ ℂ*, f_W⁻¹(t) = X_t/T² = ℂ²/(Λ_t + ℝe₂ + ℝe₃) is a 2-torus with H₁ = ℤa ⊕ ℤb, where a = [e₁] and
   b = [e₄].
2. M → {ε/2 < |t| < 2/ε} is a smooth T²-bundle with monodromy a ↦ a, b ↦ b + a. Take the base point on the
   image of the zero section, and let x be the loop θ ↦ [1, 0, r₀e^{iθ}]. Then

       π₁(M) = ⟨ a, b, x | [a, b], [x, a], x b x⁻¹ = b a^{±1} ⟩.

3. The map Ψ: M → S¹ × S¹, Ψ[v, t] = (arg t, e₄*(v) mod ℤ), is well defined and continuous. It induces an
   isomorphism π₁(M)/⟨⟨a⟩⟩ ≅ ℤ² with x ↦ (1, 0) and b ↦ (0, 1).

*Proof.*
- (1) e₁, …, e₄ is a real basis, so ℝ⁴/(ℝe₂ + ℝe₃) ≅ ℝe₁ ⊕ ℝe₄.
- (2) Ehresmann's theorem applies, since f_W is a proper submersion over ℂ*. The monodromy is that of
  Lemma 3.3, and the zero section gives the semidirect product description.
- (3) Ψ is well defined by the properties of e₄*. Modulo the normal closure of a, π₁(M) is free abelian on
  b and x, and Ψ_* sends these to the standard basis. ∎

### 7.2 The local pieces

**Lemma 7.2.** Let j ∈ {0, ∞}.
1. N_j ≅ 𝒳^ε/⟨ψ̂⟩, where ψ̂(p) = κ(t(p)) · ν ψ(p).
   - κ = exp((γ₂ − iγ₃)t) at 0, and κ = exp((γ₂ + iγ₃)s) at ∞.
   - (γ₂, γ₃): D_ε → ℝ² is the real-analytic solution of γ₂ ē₂ + γ₃ ē₃ = −ē₄.
   - The action of ⟨ψ̂⟩ is free and properly discontinuous.
2. π₁(N_j) ≅ ℤ, generated by ψ̂, and H_q(N_j; ℤ) ≅ ℤ, ℤ, ℤ, 0, 0 for q = 0, …, 4. In particular χ(N_j) = 1.
3. Let a_j, b_j, x_j generate π₁(O_j):
   - a_j = [e₁] and b_j = [e₄] (resp. [−e₄]) are fibre classes;
   - x_j is the loop around the cusp on the zero section of U_j.

   Under O_j ⊂ N_j these map to 1, ψ̂ and 1.

*Proof.* (1) ℝ² acts on Ũ_j by Ã_c, and Γ ⊂ ℝ² × ⟨G₄⟩. So N_j = U_j/T² = Ũ_j/(ℝ² × ⟨G₄⟩).
- Since ē₂(t), ē₃(t) is a real basis of ℂ, the slice {w = 0} meets every ℝ²-orbit exactly once.
- The map (p, w) ↦ m_c · p, with c = c(t, w) solving w + c₂ ē₂ + c₃ ē₃ = 0, is a real-analytic submersion. It
  identifies Ũ_j/ℝ² with 𝒳^ε.
- G₄ maps (p, 0) to (νψ(p), ē₄), which lies on the ℝ²-orbit of (ψ̂(p), 0).

Since ℝ² is central in ℝ² × ⟨G₄⟩, N_j = (Ũ_j/ℝ²)/⟨G₄⟩ = 𝒳^ε/⟨ψ̂⟩.
- The slice inclusion p ↦ (p, 0) and the slice map (p, w) ↦ m_{c(t,w)} · p are smooth and equivariant.
- They induce mutually inverse smooth maps 𝒳^ε/⟨ψ̂⟩ ↔ U_j/T², since U_j → U_j/T² is a smooth principal bundle.
- So these spaces are diffeomorphic, and under the diffeomorphism f_W corresponds to t.

*Free.* Suppose ψ̂^k fixes p.
- Then Ã_c G₄^k fixes (p, 0) for some c.
- Its image in U_j is fixed by the class [c] ∈ T², so c ∈ ℤ² by freeness of the T²-action.
- Then Ã_c G₄^k ∈ Γ, and k = 0 by freeness of Γ.

*Properly discontinuous.* Let K ⊂ 𝒳^ε be compact, and put K̃ = Ã_{[0,1]²}(K × {0}), which is compact.
- Suppose ψ̂^k K ∩ K ≠ ∅. Write Ã_c G₄^k = G₂^{⌊c₂⌋} G₃^{⌊c₃⌋} G₄^k ∘ Ã_{{c}} with {c} ∈ [0, 1)².
- The element G₂^{⌊c₂⌋} G₃^{⌊c₃⌋} G₄^k ∈ Γ then moves K̃ to meet K̃.
- By Proposition 4.1(a) only finitely many k occur.

(2) By Lemma 2.6, 𝒳^ε is simply connected. So it is the universal cover of N_j, and π₁(N_j) = ⟨ψ̂⟩ ≅ ℤ.
- ψ̂ = κ′ · ψ with κ′ = ν κ(t), which has the continuous logarithm log ν + (γ₂ ∓ iγ₃)t (sign − at 0, + at ∞,
  with s in place of t). So Lemma 2.6 gives ψ̂_* = ψ_* = τ on homology.
- The Wang sequence (F2) gives short exact sequences

      0 → coker(τ_* − 1 | H_q(𝒳^ε)) → H_q(N_j) → ker(τ_* − 1 | H_{q−1}(𝒳^ε)) → 0.

- τ_* = 1 on H₀ = ℤ, H₁ = 0, and τ − 1 is injective with cokernel ℤ on ℤ[τ^{±1}]. These give:
  - H₀(N_j) = ℤ;
  - H₁(N_j) = ker(τ_* − 1 | H₀) = ℤ;
  - H₂(N_j) = coker(τ − 1 | ℤ[τ^{±1}]) = ℤ;
  - H₃(N_j) = ker(τ − 1 | ℤ[τ^{±1}]) = 0, and H₄(N_j) = 0.

(3) The loops a_j (the ζ-circle) and x_j (the loop t = re^{iθ}, ζ = 1, w = 0, on which c = 0 in the slice
map) lift to closed loops in 𝒳^ε. The fibre loop b_j lifts to a path from p to ψ̂(p). ∎

### 7.3 The effect of the logarithmic transform

**Lemma 7.3.** Let γ be the image in M, via Φ ∘ ι₀, of the loop x₀. Then Ψ_*(γ) = (1, 1), hence
γ ≡ xb (mod ⟨⟨a⟩⟩). On the other hand, x_∞ maps to x⁻¹ under ι_∞.

*Proof.* In additive coordinates the loop is θ ↦ [σ̃(s(θ)), t(θ)], with t(θ) = re^{iθ} and s(θ) a continuous
branch.
- arg t winds once.
- Since e₄* is continuous and single-valued, the change of e₄* along the loop is
  e₄*(σ̃(s(0) + 1) − σ̃(s(0))) = e₄*(e₄) = 1, by Lemma 5.1(2).
- The loop x has e₄* ≡ 0 along it.
- At ∞: ι_∞ maps zero section to zero section (§4.2), and a counterclockwise loop in s is clockwise in t. ∎

**Proposition 7.4.** π₁(W) = 1.

*Proof.* All pieces and intersections are path-connected, so van Kampen applies to W = (N₀ ∪ M) ∪ N_∞.
- By Lemma 7.2(3), gluing N_∞ imposes a = 1 and x_∞ = 1, i.e. x = 1 by Lemma 7.3. It also identifies the
  generator of π₁(N_∞) with b^{∓1}, up to a factor in ⟨⟨a⟩⟩.
- Gluing N₀ along O₀ via Φ ∘ ι₀ imposes a = 1 and γ = 1. It also identifies the generator of π₁(N₀) with b.
  (Base point changes only conjugate these relations.)
- By Lemma 7.3, γ = xbn with n ∈ ⟨⟨a⟩⟩. So 1 = γ = b in the quotient.

Thus all generators a, b, x of π₁(M), and the generators of π₁(N_j), are trivial. ∎

**Proposition 7.5.** χ(W) = 2.

*Proof.* All pieces have finitely generated homology. So the Mayer–Vietoris sequence of the open cover
{N₀, M, N_∞} gives

    χ(W) = χ(N₀) + χ(M) + χ(N_∞) − χ(O₀) − χ(O_∞).

M, O₀ and O_∞ are T²-bundles over annuli, so their Euler characteristics vanish. χ(N_j) = 1 by
Lemma 7.2(2). ∎

**Corollary 7.6.** W is a smooth homotopy 4-sphere.

*Proof.* W is closed, oriented and simply connected.
- Hence H₁ = 0, H₃ ≅ H¹ = 0, H₄ = ℤ, and H₂ ≅ H² ≅ Hom(H₂, ℤ) is free.
- Then χ(W) = 2 + b₂ = 2, so H₂ = 0.
- The degree-one collapse map W → S⁴ is a homology isomorphism of simply connected spaces, hence a homotopy
  equivalence by Whitehead's theorem. ∎

---

## 8. Proof of the Main Theorem (manuscript §8)

*Proof of the Main Theorem.*
- **The bundle is trivial.** By Corollary 6.2, Y → W is a principal T²-bundle. By Corollary 7.6,
  H²(W; ℤ²) = 0. By (F1) the bundle is smoothly trivial, so Y ≅ W × T².
- **h-cobordism.** By (F3), W is h-cobordant to S⁴: both are simply connected with zero intersection form.
  If H⁵ is such an h-cobordism, then H × T² is a 7-dimensional h-cobordism between W × T² and S⁴ × T², with
  fundamental group ℤ².
- **s-cobordism.** Its Whitehead torsion lies in Wh(ℤ²) = 0 (F4). By the s-cobordism theorem (F4),
  Y ≅ W × T² ≅ S⁴ × T².
- **Transport.** Let F: S⁴ × T² → Y be a diffeomorphism and J_Y the complex structure of Y. Then F*J_Y is an
  almost complex structure on S⁴ × T². It is integrable, since F is a biholomorphism from
  (S⁴ × T², F*J_Y) to Y; equivalently, N_{F*J_Y} = F*N_{J_Y} = 0. ∎

### 8.1 The quotient W is the standard S⁴

**Definitions** (following [Matsumoto85, §2]).
- A *(good) torus fibration* is a proper surjective smooth map f: M⁴ → B² of oriented manifolds with two
  properties:
  - near every point it is given, in local complex coordinates, by ξ = z₁^m z₂^n or ξ = z̄₁^m z₂^n;
  - it is a smooth T²-bundle off a discrete set of critical values.
- A singular fibre of type I₁^± is an immersed 2-sphere with one transverse self-intersection of sign ±1.
- A *twin* singular fibre (Tw) consists of two embedded 2-spheres meeting transversely in two points, of
  opposite signs.

**Proposition 8.1.** W is diffeomorphic to S⁴. Consequently Y ≅ S⁴ × T² as principal T²-bundles, and the
diffeomorphism of the Main Theorem can be chosen T²-equivariant.

*Proof.* Fix an orientation of W.

*Step 1: f_W is a torus fibration with two singular fibres, each of type I₁⁺ or I₁⁻.*
- Over ℂ*, f_W is a smooth T²-bundle (Lemma 7.1).
- By the proof of Lemma 7.2(1), there is a diffeomorphism N_j ≅ 𝒳^ε/⟨ψ̂⟩ under which f_W corresponds to t
  (to s at ∞). The quotient map 𝒳^ε → 𝒳^ε/⟨ψ̂⟩ is a covering. So the holomorphic chart (x₀, y₀) of U₀ near
  P₀ gives smooth complex coordinates on W near the singular point of f_W⁻¹(0), in which f_W = x₀y₀.
- If these coordinates are positively oriented for W, then f_W = z₁z₂ with (z₁, z₂) = (x₀, y₀). Otherwise
  f_W = z̄₁z₂ with the positively oriented coordinates (z₁, z₂) = (x̄₀, y₀).
- Elsewhere on the singular fibre dt ≠ 0, so f_W is locally ξ = z₁.
- The singular fibre is the image of the chain t⁻¹(0) modulo ψ̂, i.e. D₀ with P₋₁ and P₀ identified. This is
  an immersed 2-sphere with one transverse self-intersection. It is not multiple, since f_W = x₀y₀ is
  reduced. The same holds at ∞.

*Step 2: merging the nodes.* Let D ⊂ ℙ¹ be a closed disc containing 0 and ∞ in its interior and no other
critical value, e.g. the complement of a small open disc around t = 1.
- The monodromy around ∂D is trivial, since ∂D also bounds a disc without critical values.
- By [Matsumoto85, Thm. 4.1], f_W|_{f_W⁻¹(D)} can be deformed, without changing it near the boundary, into a
  torus fibration with a single singular fibre, of type Tw with divisor R + S.
- This gives a torus fibration W → S² whose only singular fibre is a non-multiple twin fibre.

*Step 3: Iwase's theorem* (as stated in [Matsumoto85, §1]).
> Let M → S² be a torus fibration whose singular fibres are non-multiple twin fibres, and assume
> e(M) ≠ 0. Then the diffeomorphism type of M is determined by π₁(M), e(M), w₂(M) and the parity of the
> intersection form on H₂(M; ℤ).

- S⁴ is the total space of such a fibration, with a single twin fibre [Matsumoto82], [Matsumoto85, §2, §4].
- Both W and S⁴ have π₁ = 1 and e = 2 (Propositions 7.4, 7.5).
- Both have H²(·; ℤ) = 0 and H²(·; ℤ/2) = 0 (Corollary 7.6), so w₂ = 0 and the intersection form is zero,
  hence even.
- Therefore W ≅ S⁴.

The last statement follows from Y ≅ W × T² (§8). ∎

**Remark 8.2.** Proposition 8.1 identifies the total space W with S⁴. It does **not** show that f_W is
equivalent, as a fibration, to Matsumoto's fibration of S⁴, although the two have the same singular fibres,
monodromy and twist class. It gives a second proof of the Main Theorem that replaces (F3)–(F4) by Iwase's
classification. The proof in §8 uses only (F1)–(F4).

*Open point recorded in the ledger (2026-10-01-TSG-6, flag 1).* Step 1 shows only that each node is of type
I₁⁺ or I₁⁻. That the two nodes have **opposite signs and the same vanishing cycle** a = [e₁] follows from the
trivial monodromy around ∂D, since the local monodromies are the transvections T and T⁻¹ along a. The opposite
signs alone also follow from σ(W) = 0 (Remark 10.2). The flag suggests stating this before invoking [Matsumoto85, Thm. 4.1],
whose exact hypotheses the reviewer could not consult.

---

## 9. Further properties: proof of Theorem B (manuscript §9)

**(1)** This is Proposition 5.3.

**(2)** The T²-statements are Proposition 6.1 and §8. For (z₀, w₀) ∈ ℂ², the formulas

    on Y°:    [ζ, w, t]  ↦ [e^{iz₀} ζ, w + w₀, t],
    on Ũ₀:    (p, w)     ↦ (e^{iz₀} · p, w + w₀),
    on Ũ_∞:   (p, w′)    ↦ (e^{iz₀ + w₀ s} · p, w′ + w₀ s²)

commute with all generators and with Φ, and are compatible with Θ (the computation of Lemma 4.2). So they
define a holomorphic ℂ²-action, which is transitive on each X_t.
- *The kernel.* If (z₀, w₀) acts trivially, then (z₀, w₀) ∈ Λ_t for every t. Write (z₀, w₀) = Σ n_j e_j(t)
  with locally constant integers. The second coordinate n₂(1 − t²) + in₃(1 + t²) + 2n₄t must be constant in t,
  so n₂ = n₃ = n₄ = 0.
- Hence the kernel is ℤ(2π, 0), and Aut⁰(Y) ⊃ ℂ²/ℤ(2π, 0) ≅ ℂ* × ℂ.

**(3) c₁ = 0.** Let O ⊂ X_t be a T²-orbit. Under Y ≅ W × T² it is {pt} × T², which generates
H₂(Y; ℤ) ≅ ℤ.
- There is an exact sequence of holomorphic vector bundles on X_t,

      0 → TX_t → TY|_{X_t} → N_{X_t} → 0,     N_{X_t} = f*T_tℙ¹|_{X_t} ≅ 𝒪_{X_t}.

- TX_t ≅ 𝒪² is trivialised by translations. By additivity of c₁, c₁(TY|_{X_t}) = 0.
- So c₁(Y) · [O] = 0, and c₁(Y) = 0 because H²(Y; ℤ) ≅ ℤ is free.

The sequence splits as a sequence of smooth complex vector bundles, but not holomorphically [check 13]. Its
extension class is the Kodaira–Spencer class of f at t. That class is nonzero, because the period ratio
p₁₃/p₁₂ = i(1 + t²)/(1 − t²) of Lemma 3.5 has derivative 4it/(1 − t²)² ≠ 0; at t = ±1 use p₁₂/p₁₃ instead.
Only c₁ is used above.

**(3) The canonical bundle** [check 6]. Consider the meromorphic 3-form Ω given on Y° by
dz ∧ dw ∧ dt = (1/i)(dζ/ζ) ∧ dw ∧ dt.
- It is invariant under the g_j and Φ, which have the form (z, w, t) ↦ (z + μ(t), w + u(t), t).
- On Ũ₀, since (dζ/ζ) ∧ dw ∧ dt = −(dζ/ζ) ∧ dt ∧ dw, it equals (1/i) dx_n ∧ dy_n ∧ dw (Lemma 2.4(6)). This is
  holomorphic and nowhere zero.
- Near ∞, Ω = −s⁻⁴ dz″ ∧ dw′ ∧ ds with ζ′ = e^{iz″}, and dz″ ∧ dw′ ∧ ds is a nowhere-zero holomorphic form on
  Ũ_∞.

Hence div(Ω) = −4F_∞ and K_Y ≅ 𝒪_Y(−4F_∞) = f*𝒪(−4), because F_∞ = f*[∞] is reduced. K_Y is not torsion:
- if f*𝒪(−4m) were trivial for some m ≥ 1, it would have a nonzero section;
- but f_*𝒪_Y = 𝒪_{ℙ¹} (connected fibres) gives H⁰(Y, f*𝒪(−4m)) = H⁰(ℙ¹, 𝒪(−4m)) = 0.

**(4)** Under Y ≅ W × T², a smooth fibre is π⁻¹(f_W⁻¹(t)) = f_W⁻¹(t) × T². Since H₂(W) = 0, [F] = 0 in
H₄(Y; ℤ).
- *Not balanced.* A balanced metric ω_Y (dω_Y² = 0) would give 0 < ∫_F ω_Y² = ⟨[ω_Y²], [F]⟩ = 0.
- *Not Kähler.* Kähler metrics are balanced. Alternatively, H²(Y; ℝ) = ℝa with a² = 0.
- Consistency check: c₁(Y) = −c₁(f*𝒪(−4)) = 4 PD[F] = 0.

**(5) ℳ(Y) = f*ℂ(ℙ¹).** Pulling back by f gives f*ℂ(ℙ¹) ⊂ ℳ(Y). Conversely, let h ∈ ℳ(Y), and let ℰ be the
countable set of Proposition 3.6.
- *Over a disc.* Let B ⊂ ℂ* be an open disc. Pull h back along the holomorphic covering ℂ² × B → f⁻¹(B) to a
  meromorphic function H, invariant under Λ_t.
- *Restriction to fibres.* The irreducible components of the polar set of H form a locally finite family.
  Every fibre ℂ² × {t} meets the compact set {0} × K for t in a compact K ⊂ B. So only a discrete set of
  fibres lies in the polar set.
- For every other t, H(·, t) is a Λ_t-periodic meromorphic function on ℂ², i.e. a meromorphic function on X_t.
  For t ∉ ℰ it is constant (Proposition 3.6).
- *H depends only on t.* The set of such t is dense in B. So ∂_zH and ∂_wH vanish on a dense subset of the
  open set where they are holomorphic, hence identically.
- Since the complement of the polar set is connected, H = c_B ∘ pr_B for a meromorphic function c_B on B.
  (Take c_B = H(z₀, w₀, ·) for (z₀, w₀) such that {(z₀, w₀)} × B is not contained in the polar set.)
- *Gluing.* The c_B glue to a meromorphic function c on ℂ* with h = f*c on Y°.
- *The ends.* Let σ₀ be the holomorphic section of U₀ → D_ε of Proposition 4.1(e).
  - Off a discrete set, σ₀(t) is a point where h = f*c is holomorphic. So σ₀(D_ε) is not contained in the
    polar set of h, and h ∘ σ₀ is meromorphic on D_ε.
  - It equals c on D*_ε, so c is meromorphic at 0.
  - The same argument with the section of U_∞ applies at ∞.

Thus c ∈ ℂ(ℙ¹), and h = f*c, since both sides are meromorphic on Y and agree on the dense open set Y°. Hence
ℳ(Y) = f*ℂ(ℙ¹) and a(Y) = 1. ∎

**Remark 9.1 (deformations)** [checks 7, 14]. Replace g₂, g₃, g₄ by

    (e^{λt} ζ, w + 1 − t²),   (e^{−iλt} ζ, w + i(1 + t²)),   (e^{iβ} t ζ, w + 2t),

with λ > 0 and β ∈ ℂ. The determinant of Lemma 3.2 becomes 2π(1 + r²) D(r), where r = |t|, b = Im β and

    D(r) = 2λ r² + (b − ln r)(1 − r²).

*The admissible region.* D is positive for every r > 0 if and only if

    |Im β − λ| < √(λ(λ+1)) + arsinh √λ.                                              (9.1)

- D(1) = 2λ > 0.
- For 0 < r < 1, D(r) > 0 iff b > g(r), where g(r) = ln r − 2λr²/(1 − r²).
  - g tends to −∞ at both ends of (0, 1).
  - g′(r) = 1/r − 4λr/(1 − r²)² vanishes only where 1 − r² = 2√λ r, i.e. at r_λ = √(λ+1) − √λ.
  - There r_λ²/(1 − r_λ²) = r_λ/(2√λ), so

        S_λ := max_{(0,1)} g = ln r_λ − √λ r_λ = λ − √(λ(λ+1)) − arsinh √λ.

- For r > 1 substitute r = 1/ρ. Since 1/(1 − ρ²) = 1 + ρ²/(1 − ρ²), D(1/ρ) > 0 iff b < 2λ − g(ρ).
- So D > 0 on (0, ∞) iff S_λ < b < 2λ − S_λ, which is (9.1).

The region (9.1) contains Im β = 0 for every λ > 0, since √(λ(λ+1)) > λ. It is not all of {λ > 0} × ℂ: at
λ = 1, β = −2i one has D(½) = −1 + ¾ ln 2 < 0. So D vanishes for some r ∈ (½, 1), and the four period vectors
fail to span a lattice there. Outside (9.1) the construction breaks down.

*Modified proofs.* For (λ, β) satisfying (9.1), the proofs go through with these changes:
- use Θ_λ(ζ, w, t) = (ζ e^{λw/t}, w/t², 1/t) at ∞;
- there the filling data are α₂ = λs, α₃ = iλs, and ν = e^{−iβ−2λ};
- the twist is translation by (πs(s − 1) + βs, 2ts);
- the T²-action has multipliers e^{λ(c₂ − ic₃)t} near 0 and e^{λ(c₂ + ic₃)s} near ∞.

This gives a family of complex structures on S⁴ × T², parametrised by the region (9.1). The original
construction is (λ, β) = (1, 0). It has not been determined which members are biholomorphic.

---

## 10. Remarks and consistency checks (manuscript §10)

**Remark 10.1 (Leray spectral sequence: a heuristic check; not used in the proof).** The Betti numbers follow
rigorously from Y ≅ W × T².
- Rationally, the monodromy-invariant covectors are e₂*, e₃*, e₄*.
- In the Leray spectral sequence of f, the twist class [e₄] gives d₂(e₄*) = ±ω_{ℙ¹}, while
  d₂(e₂*) = d₂(e₃*) = 0.
- This yields b₁(Y) = 2, and the class of H²(ℙ¹) dies. The surviving class e₂* ∧ e₃* accounts for b₂(Y) = 1.

A complete Leray computation would need the sheaves R^q f_*ℚ at the singular fibres; it is not carried out.

**Remark 10.2 (signature).** Matsumoto's local signature of a genus-one node is ∓⅔ for positive and negative
nodes. Since σ(W) = 0, the two nodes of f_W must have opposite signs, in accordance with Remark 3.4.

**Remark 10.3 (polarised families have no cusps over ℂ*).** Let {A_t}_{t∈ℂ*} be a holomorphic family of
polarised abelian surfaces.
- Its period map lifts to the universal cover, ℂ → 𝔥₂. Siegel space is biholomorphic to a bounded domain, so
  the lift is constant by Liouville's theorem.
- Hence the monodromy fixes a point of 𝔥₂ and has finite order. After a finite base change the family is
  trivial near 0 and ∞; in particular it has no cusp there.

So the family of §3, whose monodromy T has infinite order, admits no polarisation over ℂ*. Its period map
takes values in the space of complex structures on ℝ⁴, which is not hyperbolic, and is non-constant.
- This argument concerns polarisations of the whole family. It does not decide whether individual fibres are
  projective, and some are: X₁ is an abelian surface.
- Non-algebraicity of the very general fibre is a statement about integral Hodge classes, proved directly in
  Proposition 3.6.

**Remark 10.4 (a necessary condition).** Let f: Y → ℙ¹ be a fibration by complex 2-tori with finitely many
singular fibres and connected fibres, with Y ≃ S⁴ × T². Assume that the restriction
H⁰(ℙ¹, R¹f_*ℚ) → H¹(F_t; ℚ)^{π₁} to a smooth fibre is injective.
- A section that vanishes on one smooth fibre vanishes on the whole regular locus. So the assumption holds
  whenever each singular fibre has a neighbourhood N with H¹(N; ℚ) → H¹(F_t; ℚ) injective.
- That is the case for the fibres I₁ × E used here, since H₁(F_t) → H₁(I₁ × E) is onto.

The assumption is used in the second and third steps below.
1. The base contributes no H¹, so b₁(Y) = 2 means ker(d₂: H⁰(R¹) → H²(ℙ¹; ℚ)) has rank 2.
2. If d₂ = 0, then f*[pt] ≠ 0 in H². The cup product of the two classes in H¹ is then a second, independent
   class: it projects to ψ₁ ∧ ψ₂ ≠ 0 in E_∞^{0,2}. That gives b₂ ≥ 2, a contradiction.
3. The other ranks give b₁ ≠ 2, so rank H⁰(R¹) = 3. Hence every monodromy fixes a common rational hyperplane
   Π ⊂ H¹(F_t; ℚ), or acts trivially.

Monodromies of a holomorphic family preserve the complex orientation of the fibres, so they have
determinant 1.
- An element T of SL₄ fixing Π pointwise satisfies T* = I + u ⊗ ℓ with ℓ|_Π = 0 and det T = 1 + ℓ(u) = 1. So
  u ∈ Π, and T is a transvection.
- Dually, every monodromy is a transvection λ ↦ λ + u_γ(λ) λ₀ with a common λ₀.
- A transvection I + N with N² = 0 satisfies (I + N)^k = I + kN. So nontrivial monodromies of finite order,
  as in the (3,4,∞) family, cannot occur.

In this construction λ₀ = e₁.

**Remark 10.5 (questions).**
1. W ≅ S⁴ (Proposition 8.1). Is f_W equivalent to Matsumoto's fibration? Can one write an explicit
   diffeomorphism W ≅ S⁴, and hence an explicit integrable J on the round S⁴ × T²? (Editorial note, not in the
   manuscript: this is the goal of the symbolic campaign.)
2. Compute the Dolbeault cohomology and the Kuranishi space of Y.
3. Twist by n[e₄] instead of [e₄], and call the quotient W_n.
   - The arguments of §7 give χ(W_n) = 2, and π₁(W_n) = ℤ/n for n ≠ 0. So for n ≠ 0, W_n is a rational
     homology 4-sphere.
   - Untwisted, π₁(W₀) = ℤ, so b₁ = b₃ = 1 and b₂ = 2.
   - Steps 1–2 of the proof of Proposition 8.1 apply to every W_n. So Iwase's theorem determines W_n from
     π₁(W_n), w₂(W_n) and the parity of its intersection form. Which 4-manifolds are they?
4. Calabi's question for V² of genus ≥ 2 remains open as far as is known. For the products of spheres of
   dimensions 10 and 12 that remain open, see §11.7.

---

## 11. Products of spheres (manuscript §11)

Unlike the Main Theorem, **Theorem C and Corollary D depend on recent preprints.** Products of spheres carry
their standard smooth structures, S^n always has n ≥ 1, and products of spheres of dimension n correspond to
partitions of n.

### 11.1 Inputs

Besides the classical constructions (complex projective spaces, elliptic curves, Hopf and Calabi–Eckmann
manifolds, and products), the following inputs are used.
- **(I1)** S⁶ admits a complex structure [AC26, Engel26]. X_sph ≅ S⁶ denotes the threefold of [AC26].
- **(I2)** S² × S⁴ admits a complex structure [GFL26, LLZ26].
- **(I3)** *The package (LLZ).* There is a compact connected complex threefold X₊, a diffeomorphism
  X₊ ≅ S² × S⁴, and a holomorphic line bundle L₊ → X₊ whose first Chern class generates H²(X₊; ℤ) ≅ ℤ
  [LLZ26].

Dependencies among the inputs:
- (I3) implies (I2).
- (I3) assumes more than (I1). It starts from X_sph with:
  - Aut⁰(X_sph) ≅ ℂ*;
  - normal weights ±1 on the fixed curve of this action;
  - the structure of the cusp fibre, with its ℂ*-invariant double curves.

  See [AC26, Thm. 4.5, Prop. 4.6, Props. 9.23–9.24, Rem. 9.25].
- Engel's independent account [Engel26] verifies the construction of X_sph, but treats these auxiliary
  properties only in a remark.
- (I2) also rests on [AC26] by either route: [GFL26] reuse its analytic framework, and [LLZ26] start from
  X_sph.

So none of (I1)–(I3) is independent of [AC26].

### 11.2 Elliptic bundles over a complex base (Calabi–Eckmann)

**Proposition 11.1.** Let B be a compact connected complex manifold of dimension n, let L₁, L₂ → B be
holomorphic line bundles, and fix τ ∈ ℂ ∖ ℝ. Let 𝒫 = L₁^× ×_B L₂^× = {(v₁, v₂) : v₁ ≠ 0 ≠ v₂}, and let ℂ act by
t·(v₁, v₂) = (e^t v₁, e^{τt} v₂). Then:
1. Z = 𝒫/ℂ is a compact complex manifold of dimension n + 1;
2. Z → B is a holomorphic fibre bundle with fibre E_τ = ℂ/(ℤ + τℤ) and structure group E_τ acting by
   translations;
3. for any hermitian metrics, Z is diffeomorphic to the principal T²-bundle
   𝒬 = S(L₁) ×_B S(L₂) = {(v₁, v₂) ∈ 𝒫 : |v₁| = |v₂| = 1}.

This is the Calabi–Eckmann construction [CalabiEckmann53] over an arbitrary compact complex base; see
[Morimoto63, Hofer93, GoldsteinProkushkin04, GGP08]. The manuscript proves it in full, with G = (ℂ*)² and
H_τ = {(e^t, e^{τt})}, in three lemmas.
- **Lemma 11.2.**
  - t ↦ (e^t, e^{τt}) is injective, since τ ∉ ℝ.
  - ρ: G → E′_τ = ℂ/2πi(ℤ + τℤ), (w₁, w₂) ↦ [log w₂ − τ log w₁], is a well-defined surjective holomorphic
    homomorphism with kernel H_τ.
  - So H_τ is closed, and G/H_τ ≅ E′_τ ≅ E_τ via z ↦ z/2πi.
- **Lemma 11.3.** 𝒫 is a holomorphic principal G-bundle. Gluing the charts U_α × E′_τ by
  (b, e) ↦ (b, e + ρ(g_αβ(b))) gives a compact complex manifold Z′. The map 𝒫 → Z′ induced by ρ is a
  surjective submersion whose fibres are exactly the H_τ-orbits, so Z = 𝒫/H_τ ≅ Z′.
- **Lemma 11.4.** Let μ(v) = (log|v₁|, log|v₂|) and A_τ(t) = (Re t, Re(τt)), so that det A_τ = −Im τ ≠ 0.
  - Then μ(t·v) = μ(v) + A_τ(t) (11.1), so 𝒬 = μ⁻¹(0) is a compact submanifold and a principal T²-bundle.
  - 𝒬 → Z is injective by (11.1).
  - Its inverse is induced by v ↦ t(v)·v with t(v) = −A_τ⁻¹μ(v).

**Remark 11.5 (one line bundle).** If L₂ is trivial, then Z ≅ L₁^×/λ^ℤ with λ = e^{2πi/τ}, and
|λ| = exp(2π Im τ/|τ|²) ≠ 1. Conversely, every λ with |λ| ≠ 1 arises, from τ = 2πi/log λ. Then
Proposition 11.1(3) gives L₁^×/λ^ℤ ≅ S(L₁) × S¹. For L₁ = 𝒪_{ℙ^p}(−1) this is the Hopf manifold
(ℂ^{p+1} ∖ 0)/λ^ℤ ≅ S^{2p+1} × S¹.

### 11.3 Circle bundles

The standard facts used are as follows.
- **(C1)** Every complex line bundle has a hermitian metric. S(L) is a smooth principal U(1)-bundle whose
  Chern class is c₁(L), and rescaling the metric gives an isomorphic circle bundle.
- **(C2)** Smooth principal U(1)-bundles over a manifold are classified by their Chern class [Husemoller].
- **(C3)** Pullback of a bundle along a diffeomorphism of the base has a diffeomorphic total space.
- **(C4)** For 𝒪(−1) → ℙ^q with the induced metric, S(𝒪(−1)) ≅ S^{2q+1}.
  - For q = 1 this is the Hopf fibration, and over ℙ¹ S(𝒪(1)) is also S³.
  - For q = 0, S(𝒪(−1)) = S¹.

**Lemma 11.6.** For L_j = pr_j*L′_j on B₁ × B₂, S(L₁) ×_B S(L₂) ≅ S(L′₁) × S(L′₂).

**Lemma 11.7.** If X ≅ S² × S⁴ and c₁(L) generates H²(X; ℤ), then S(L) ≅ S³ × S⁴.
- *Proof.* Pull back along a diffeomorphism φ: S² × S⁴ → X. Then φ*c₁(L) = ±pr₁*c₁(𝒪_{ℙ¹}(1)). By (C2),
  φ*S(L) ≅ pr₁*S(𝒪(±1)), whose total space is S(𝒪(±1)) × S⁴ ≅ S³ × S⁴ by (C4). Conclude with (C3). ∎

### 11.4 The complex S² × S⁴ of Liang–Liu–Zhang

The construction of [LLZ26], step by step:
1. *The source.* X_sph has Aut⁰ ≅ ℂ*, a fixed smooth rational curve C with normal weights ±1, and a
   ℂ*-invariant double curve D₋ of the cusp fibre meeting C in two points Q, Q′ [LLZ26, §1].
2. *Blow-up and quotient.* Put X̃ = Bl_C(X_sph), with exceptional divisor ℰ.
   - The involution −1 ∈ ℂ* lifts to an involution ι̃, locally ι̃(z, x, t) = (z, −x, t). It fixes ℰ = {x = 0}
     pointwise and acts freely elsewhere.
   - ϖ: X̃ → X₋ = X̃/⟨ι̃⟩ is a double cover branched along ℱ = ϖ(ℰ), and ϖ|_ℰ: ℰ → ℱ is a biholomorphism
     [LLZ26, §2.1].
3. *The line bundle on X₋.* L₋⁻¹ is the anti-invariant part of ϖ_*𝒪_{X̃}.
   - L₋^{⊗2} ≅ 𝒪_{X₋}(ℱ) [LLZ26, §4.1].
   - ϖ*L₋ ≅ 𝒪_{X̃}(ℰ) (11.2) [LLZ26, (4)].
4. *The flop.* C₋ = ϖ(D̃₋) is a smooth rational curve with normal bundle 𝒪(−1) ⊕ 𝒪(−1) [LLZ26, Lemma 3.1].
   - The analytic Atiyah flop along C₋ gives X₊, with flopped curve C₊ [LLZ26, Def. 3.2].
   - X° = X₋ ∖ C₋ = X₊ ∖ C₊.
5. *The line bundle on X₊.* L₊ is the unique holomorphic line bundle with
   Bl_{C₊}*L₊ ≅ Bl_{C₋}*L₋ ⊗ 𝒪(R) on the common blow-up [LLZ26, (6)]. Hence L₊|_{X°} ≅ L₋|_{X°} (11.3).
6. *Topology.* X₊ is simply connected, with torsion-free homology and H²(X₊; ℤ) ≅ ℤ [LLZ26, Prop. 3.5]. By the
   Wall–Jupp classification [Wall66, Jupp73], X₊ ≅ S² × S⁴ orientation-preservingly
   [LLZ26, Thm. 1.2 and §5.2].

**Lemma 11.8.** c₁(L₊) generates H²(X₊; ℤ). More precisely, there is a smooth rational curve Σ ⊂ X° with
deg(L₊|_Σ) = −1.

*Proof.*
- Choose c ∈ C ∖ D₋, let Σ_c ≅ ℙ¹ be the fibre of ℰ → C over c, and put Σ = ϖ(Σ_c).
- Σ_c misses D̃₋, and also ι̃(D̃₋), since ι̃ fixes Σ_c pointwise. Since ϖ⁻¹(C₋) = D̃₋ ∪ ι̃(D̃₋), Σ ⊂ X°.
- By (11.3) and (11.2),
  deg(L₊|_Σ) = deg(L₋|_Σ) = deg((ϖ*L₋)|_{Σ_c}) = deg(𝒪_{X̃}(ℰ)|_{Σ_c}) = −1.
  The last equality holds because 𝒪(ℰ)|_ℰ = 𝒪_{ℙ(N_C)}(−1).
- So ⟨c₁(L₊), [Σ]⟩ = −1, and c₁(L₊) is ± a generator. ∎

This uses only items 2–6. It is consistent with [LLZ26, Lemma 4.1], which computes c₁(L₋)³ = 1.

**Remark 11.9 (the complex S² × S⁴ of Guo–Fang–Lu).** For [GFL26] no line bundle with primitive Chern class was
found. The cusp-fibre components satisfy [D₀] = 2x = −[D₁] (end of the proof of [GFL26, Prop. 2.2]). Applied
to 𝒪(D₀), the construction would give a lens-space factor L(2,1) in place of S³. Theorem C therefore uses
[LLZ26].

### 11.5–11.7 Proof of Theorem C, Corollary D, higher dimensions

*Proof of Theorem C.* Assume (I3) and fix τ ∈ ℂ ∖ ℝ, for instance τ = i.
- **S⁴ × S³ × S^{2q+1}, q ≥ 0.**
  - Take B = X₊ × ℙ^q, L₁ = pr₁*L₊ and L₂ = pr₂*𝒪_{ℙ^q}(−1); for q = 0, L₂ is trivial.
  - By Proposition 11.1, Z_q = (L₁^× ×_B L₂^×)/ℂ is a compact complex manifold of dimension q + 4.
  - It is diffeomorphic to S(L₁) ×_B S(L₂) ≅ S(L₊) × S(𝒪_{ℙ^q}(−1)) ≅ (S³ × S⁴) × S^{2q+1}, by Lemmas 11.6,
    11.7 (with 11.8) and (C4).
  - Dimension check: 3 + 4 + (2q + 1) = 2(q + 4).
- **(S⁴ × S³)².** Take B = X₊ × X₊ and L_j = pr_j*L₊. This gives a compact complex 7-fold diffeomorphic to
  S(L₊) × S(L₊) ≅ (S³ × S⁴)². ∎

**Remark 11.10.**
1. *The case q = 0.* By Remark 11.5, for τ = 2πi/log 2, Z₀ ≅ L₊^×/2^ℤ. (For τ = i, replace 2 by e^{2π}; ledger
   entry 2026-10-01-TSG-4.) The map v ↦ (v/|v|, log|v| mod log 2) identifies it with
   S(L₊) × S¹ ≅ S⁴ × S³ × S¹.
2. With X₊ replaced by ℙ^p and L₊ by 𝒪(−1), this is exactly the Calabi–Eckmann manifold.
3. Z_q is not Kähler, since H²(S⁴ × S³ × S^{2q+1}; ℝ) = 0.
4. K_{Z_q} ≅ π*K_{X₊×ℙ^q}, since Z_q is a holomorphic principal E_τ-bundle.
5. *Why this does not give S⁴ × T².* Every S⁴-factor here comes from X₊ ≅ S² × S⁴, paired with the S³ of a
   Hopf circle over the S²-factor. One would need a compact complex surface carrying a principal T²-bundle
   with total space S⁴ × T². None is known, and the obvious base S⁴ is not complex. So the Main Theorem is
   genuinely needed for S⁴ × T² and S⁴ × T⁴.

**Lemma 11.11.** S⁴, S⁸ and S⁴ × S⁴ admit no almost complex structure.

*Proof.* M is stably parallelisable, so all its Pontryagin classes vanish. Suppose M has an almost complex
structure J. Then:
- p₁ = c₁² − 2c₂ and p₂ = c₂² − 2c₁c₃ + 2c₄;
- c_m = e(TM) for the orientation of J, and ⟨e, [M]⟩ = χ(M), which is 2, 2 and 4 respectively.

The cases:
- *M = S⁴.* c₁ = 0 since H² = 0, so p₁ = 0 forces c₂ = 0, since H⁴(S⁴; ℤ) ≅ ℤ. But c₂ = e ≠ 0.
- *M = S⁸ or S⁴ × S⁴.* c₁ = c₃ = 0. Then p₁ = 0 forces c₂ = 0, since H⁴ is torsion-free, and p₂ = 0 then
  forces 2e = 0. But 2e ≠ 0 in H⁸(M; ℤ) ≅ ℤ. ∎

This recovers the cases of Ehresmann–Hopf (S⁴), Borel–Serre [BorelSerre53] (S⁸), and Sutherland [Sutherland65]
and Datta–Subramanian [DattaSubramanian90] (S⁴ × S⁴).

*Proof of Corollary D.*
- The products of dimension 2, 4, 6, 8 correspond to partitions: there are 2, 5, 11 and 22 of them, 40 in
  all.
- S⁴, S⁸ and S⁴ × S⁴ are not almost complex (Lemma 11.11).
- Each of the other 37 is diffeomorphic to the complex manifold in the table below. The notation:
  - 𝓗_k = (ℂ^k ∖ 0)/2^ℤ ≅ S^{2k−1} × S¹ (Hopf);
  - CE_{p,q} ≅ S^{2p+1} × S^{2q+1} (Calabi–Eckmann);
  - E is an elliptic curve;
  - X_{S²×S⁴} is either complex S² × S⁴ of (I2);
  - Y is the threefold of the Main Theorem;
  - Z₀ is from Theorem C. ∎

| dim | product | complex manifold | input |
|---|---|---|---|
| 2 | S² | ℙ¹ | classical |
| 2 | T² | E | classical |
| 4 | S⁴ | none (Lemma 11.11) | — |
| 4 | S³ × S¹ | 𝓗₂ | classical |
| 4 | S² × S², S² × T², T⁴ | (ℙ¹)², ℙ¹ × E, E² | classical |
| 6 | S⁵ × S¹, S³ × S³ | 𝓗₃, CE_{1,1} | classical |
| 6 | S³ × S² × S¹, S³ × T³ | 𝓗₂ × ℙ¹, 𝓗₂ × E | classical |
| 6 | (S²)³, (S²)² × T², S² × T⁴, T⁶ | products of ℙ¹ and E | classical |
| 6 | S⁶ | X_sph | (I1) |
| 6 | S⁴ × S² | X_{S²×S⁴} | (I2) |
| 6 | **S⁴ × T²** | **Y** | **Main Theorem** |
| 8 | S⁸, S⁴ × S⁴ | none (Lemma 11.11) | — |
| 8 | S⁷ × S¹, S⁵ × S³ | 𝓗₄, CE_{2,1} | classical |
| 8 | S⁵ × S² × S¹, S⁵ × T³ | 𝓗₃ × ℙ¹, 𝓗₃ × E | classical |
| 8 | S³ × S³ × S², S³ × S³ × T² | CE_{1,1} × ℙ¹, CE_{1,1} × E | classical |
| 8 | S³ × S² × S² × S¹, S³ × S² × T³, S³ × T⁵ | 𝓗₂ × (ℙ¹)², 𝓗₂ × ℙ¹ × E, 𝓗₂ × E² | classical |
| 8 | (S²)⁴, (S²)³ × T², (S²)² × T⁴, S² × T⁶, T⁸ | products of ℙ¹ and E | classical |
| 8 | S⁶ × S², S⁶ × T² | X_sph × ℙ¹, X_sph × E | (I1) |
| 8 | S⁴ × S² × S² | X_{S²×S⁴} × ℙ¹ | (I2) |
| 8 | S⁴ × S² × T² | X_{S²×S⁴} × E, or Y × ℙ¹ | (I2) or Main Theorem |
| 8 | S⁴ × S³ × S¹ | Z₀ | (I3), Theorem C |
| 8 | **S⁴ × T⁴** | **Y × E** | **Main Theorem** |

**Remark 11.12 (what was open).** Before 2026, the almost complex products of dimension ≤ 8 with no known
complex structure were, to the authors' knowledge, the nine products in the rows not marked "classical" or
"—". They are:
- S⁶, S⁴ × S², S⁴ × T²;
- S⁶ × S², S⁶ × T², S⁴ × S² × S², S⁴ × S² × T², S⁴ × S³ × S¹, S⁴ × T⁴.

They are almost complex: S⁶ and S² × S⁴ classically, and the others as products of almost complex manifolds or
by Staples' theorem [Staples67]. How they are settled:
- [AC26, Engel26, GFL26, LLZ26] settle S⁶, S² × S⁴ and their products with ℙ¹ and E.
- The Main Theorem settles S⁴ × T² and S⁴ × T⁴.
- Theorem C settles S⁴ × S³ × S¹.

Without the Main Theorem, S⁴ × T² would be the only open case in dimension ≤ 6, and S⁴ × T⁴ the only one in
dimension 8.

**Higher dimensions (manuscript §11.7).**
- Among products of even spheres, the almost complex ones are the products of S², S⁶ and S² × S⁴
  [Sutherland65, Thm. 3.1], as stated in [LLZ26]. Given (I1)–(I2), all of these are complex
  [LLZ26, Cor. 1.3].
- Products with an odd factor are parallelisable [Staples67], hence almost complex.
- The known complex blocks are:
  - S², S⁶, S² × S⁴;
  - Calabi–Eckmann S^{2p+1} × S^{2q+1}, which includes Hopf manifolds and T²;
  - S⁴ × S³ × S^{2q+1} and (S⁴ × S³)² (Theorem C);
  - S⁴ × T² (Main Theorem).
- No known block contains S^{2n} with n ≥ 4. Each S⁴ must be absorbed in one of four ways:
  - by an S²;
  - by two S¹'s;
  - by an S³ and one further odd sphere;
  - together with a second S⁴, by two S³'s.

The exhaustive search (`sphere_products.py`; Appendix A) finds these almost complex products with no complex
structure obtainable in this way:
- dimension 10: S⁸ × T², S⁵ × S⁴ × S¹, S⁴ × S⁴ × T²;
- dimension 12: S¹⁰ × T², S⁸ × S³ × S¹, S⁸ × S² × T², S⁸ × T⁴, S⁷ × S⁴ × S¹, S⁴ × S⁴ × S³ × S¹.

The simplest of these, S^{2n} × T² for n ≥ 4, are the higher-dimensional analogues of Calabi's question.

---

## 12. Relation to the literature (manuscript §12, summary)

To the authors' knowledge, no complex structure on S⁴ × T², or on V² × S⁴ with V of positive genus, was
previously known. Calabi's question for V = T² is recorded as open in [Calabi58] and [MO89798].

**Non-integrability results.** [EuhSekigawa13, EuhSekigawa16] prove non-integrability only for *orthogonal*
almost complex structures with respect to special metrics.

**The 2026 constructions.**
- [AC26, Engel26, GFL26, LLZ26], and the subsequent constructions the authors know of, all have finite
  fundamental group.
- None treats S⁴ × T², and none uses two cusps of opposite chirality or a multiplicity-one logarithmic
  transform.
- They are cited for their stated results and were not independently verified. The Main Theorem and
  Theorem B use none of them.

**Classical families.** The manuscript names, for each family, the invariant that rules it out. For
comparison, π₁(S⁴ × T²) = ℤ², π₂ = 0, and the Betti numbers are 1, 2, 1, 0, 1, 2, 1.

| family | why it cannot contain S⁴ × T² | source |
|---|---|---|
| Calabi–Eckmann, Hopf, Brieskorn–Van de Ven | in real dimension 6 these are products of two odd (homotopy) spheres, so b₂ = 0 | — |
| LVM manifolds | π₁ = ℤ² would force N ≅ X_B × T² with X_B a closed 2-connected 4-manifold, so χ(X_B) = 2; but every link has χ = 0 | [BosioMeersseman06, §12, Lemmas 0.9, 0.10] |
| maximal torus actions | the standard T⁴-action on S⁴ × T² does not occur, since two of its characteristic submanifolds meet in a disconnected set; this excludes only that action | [Ishida19] |
| twistor spaces | Leray–Hirsch gives b₃(Z) = b₃(M) + b₁(M), so b₁(Z) = 2 forces b₃(Z) ≥ 2, but b₃(S⁴ × T²) = 0 | — |
| Vaisman manifolds | b₁ is odd | [Vaisman82] |
| Kato manifolds | π₁ ≅ ℤ | [Kato77], [IOP19, §1] |
| compact homogeneous | fibres over a rational homogeneous Q with compact complex parallelisable fibre F; if Q is a point, the universal cover of the Lie group is contractible or has H³ ≠ 0; otherwise F is a torus and the exact sequence of homotopy groups gives a contradiction | [Tits62, Wang54] |
| symplectic | H² is generated by a with a² = 0, so any complex structure is non-Kähler | — |

**Closest antecedent.** [DKZ17] make the positive half of the Matsumoto–Fukaya fibration holomorphic on
ℝ⁴ ⊂ S⁴. The present construction realises *both* nodes, after fibring by T², with very generally
non-algebraic 2-tori, in complex dimension three.

---

## 13. Status

- **Proof status.** A full review of the manuscript, with every computation re-derived by hand, found no
  mathematical error (ledger entry 2026-10-01-TSG-6). The three local imprecisions it found were fixed in
  TSG-3 to TSG-5.
- **Open flags from that review.**
  1. Proposition 8.1, Steps 1–2: state explicitly that the nodes have opposite signs and a common vanishing
     cycle (see §8.1 above).
  2. A notation clash: Φ and Ψ in Lemma 11.4 versus the log transform Φ and the map Ψ of Lemma 7.1.
  3. PDF bookmarks.
  4. Pinpoint citations that the reviewer could not check against their sources.
  5. The open points of TSG-1: the unused macro `\hflogo`, and `amsmath` and `multirow` loaded twice.
- **Dependencies.**
  - The Main Theorem and Theorem B are unconditional. They use only (F1)–(F4) and the computations above.
  - Proposition 8.1 additionally uses [Matsumoto85, Iwase84].
  - Theorem C and Corollary D depend on the 2026 preprints through (I1)–(I3).
- **Constructivity.** The diffeomorphism Y ≅ S⁴ × T² is obtained by classification, not by a formula. The
  explicit diffeomorphism and the transported J on the round S⁴ × T² are the subject of the symbolic
  campaign.

---

## Appendix A. Verification map

All checks use exact `sympy` arithmetic, with no floating-point tolerance. The few inequalities are
decided exactly at the stated points; the arguments that extend them to all r are by hand, as noted below. The two scripts live in
`manuscript_draft/code/`. `extra/s4t2_analytic/` has the environment, a runner and the reference outputs.

| check in `verify_construction.py` | statement checked | manuscript |
|---|---|---|
| 1 | det_ℝ(e₁, …, e₄) = 2π(1 + r²)(2r² − (1 − r²) ln r), independent of arg t; the sign of (1 − r²) ln r is spot-checked at r = 1, ½, 2, while the general sign argument for positivity is by hand | Lemma 3.2 |
| 2 | monodromy e₄ ↦ e₄ + e₁ | Lemma 3.3 |
| 3 | conjugation by Θ: g₂ ↦ h₂, g₃ ↦ h₃, g₄ ↦ h₄⁻¹ | Lemma 4.2 |
| 4 | σ̃(s + 1) − σ̃(s) = e₄(s) | Lemma 5.1(2) |
| 5 | T²-multiplier e^{(c₂+ic₃)s} → 1 and shift → −c₂ + ic₃ at ∞ | Proposition 6.1 |
| 6 | dz″ ∧ dw′ = s² dz ∧ dw on fibres; dz ∧ dw ∧ dt = −s⁻⁴ dz″ ∧ dw′ ∧ ds | Theorem B(3), canonical bundle |
| 7 | deformed determinant 2π(1 + r²)(2λr² + (Im β − ln r)(1 − r²)) | Remark 9.1 |
| 8 | toric charts: t = xy, ζ = xⁿyⁿ⁺¹, the shift ψ, chart overlaps, (dζ/ζ) ∧ dt = −dx ∧ dy | Lemma 2.4(1),(5),(6) |
| 9 | ω(ē₂, ē₃) = 1 − \|t\|⁴ near 0 and \|s\|⁴ − 1 near ∞ | Remark 3.4, Lemma 4.2 |
| 10 | det_ℝ = −2πPω(ē₂, ē₃) for arbitrary symbolic data | Lemma 2.2 |
| 11 | exact heights P₀, P_∞; the derivative formula of the common denominator; bounds 1.000044 and 1.434339 at r = 10⁻²; failure at r = 1/16. The sign of the derivative, hence monotonicity and the bounds for all r < 10⁻², is by hand | Remark 4.3 |
| 12 | the six periods; no nonzero integral relation; X₁ as a product lattice | Lemma 3.5, Proposition 3.6(2) |
| 13 | d/dt (p₁₃/p₁₂) = 4it/(1 − t²)² | Theorem B(3), Kodaira–Spencer class |
| 14 | critical point r_λ, maximum S_λ, counterexample D(½) = −1 + ¾ ln 2 at λ = 1, β = −2i | Remark 9.1, (9.1) |
| 15 | dual bases (−1, n+1), (1, −n) and Aᵀm_i^{(n+1)} = m_i^{(n)} | Lemma 2.4(1),(5) |

`sphere_products.py` (manuscript Appendix A.1):
- enumerates the products of spheres of dimensions 2–12;
- decides almost complexity (Sutherland when all factors are even, Staples otherwise);
- searches all decompositions into the blocks of §11.7;
- records the inclusion-minimal sets of inputs;
- asserts the dimension-8 classification of Corollary D, including that S⁴ × S³ × S¹ needs (I3) and that
  S⁴ × T⁴ needs the Main Theorem;
- produces the lists of §11.7 for dimensions 10 and 12.

**Running the checks**, from the repository root:

```bash
python3 manuscript_draft/code/verify_construction.py     # ends with "all checks passed"
(cd manuscript_draft/code && python3 sphere_products.py)  # ends with the dimension-8 summary
extra/s4t2_analytic/run_all.sh                            # both, plus regenerating the reference outputs
```
