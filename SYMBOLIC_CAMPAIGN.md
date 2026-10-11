# The symbolic campaign: a rigorous, constructive, symbolic-elementary proof that the round S⁴ × T² carries the paper's complex structure

*Entry point for collaborators and their agents. Opened 2026-10-07 (TSG). Rules: `AGENTS.md`, in particular §8.
Read this file, then `docs/s4t2_symbolic_programme.md` (the programme) and, for the construction itself,
`docs/s4t2_explicit_diffeo.md` §12.*

## 0. The goal, and what "symbolic-elementary" means

The paper proves that the compact complex threefold Y is diffeomorphic to S⁴ × T² non-constructively. The campaign
replaces the diffeomorphism by an explicit one, F: Y → S⁴ × T² with explicit inverse G, so that the complex
structure J = G*J_Y lives on the standard round S⁴ × T² ⊂ ℝ⁵ × ℝ² × ℝ², given by formulas. The target is a
**certificate** of the statement

> *J is a smooth (1,1)-tensor field on all of S⁴ × T² with J² = −1 and N_J = 0, hence (Newlander–Nirenberg) a
> complex structure,*

whose correctness a reader can check by (i) verifying that the formulas are what the documents say, (ii) running
short scripts, and (iii) reading short elementary arguments. Concretely, every step of the proof must be one of:

* **(a) an identity**, verified by a sympy script that prints `ALL PASS`;
* **(b) a manifest inequality**: a sum of squares, a convex combination of terms of known sign, a one-variable
  monotonicity with the derivative written out;
* **(c) an elementary inequality with explicit constants**, written step by step, every constant stated, every
  series bounded by a series with positive terms (model: `docs/s4t2_symbolic_programme.md` §4.2);
* **(d) a classical theorem, cited** (Newlander–Nirenberg; the zero count of a degree-1 theta function; the
  implicit function theorem).

No interval arithmetic, no proof assistant, no appeal to numerics. Numerical checks remain as a **safety net**
and as the source of margins, never as proof. This is the strongest form of certificate the statement admits:
the theorem says that an explicit map is a diffeomorphism, which is quantitative, so some inequalities are
unavoidable; the campaign reduces them to type (b)/(c) and keeps them few and short.

**Logical structure (fixed, see `docs/s4t2_explicit_diffeo.md` §13.1 and `docs/s4t2_symbolic_programme.md` §7).**
For J = G*J_Y the identities J² = −1 and N_J = 0 hold for *every* smooth G with invertible differential
(naturality). They are consistency checks of an implementation, not of the theorem, and are verified once,
generically. **The entire content is that G is well defined**: smooth on each piece, DG invertible (certified by an
explicit left inverse F with F∘G = id), pieces agreeing on overlaps (identities through Y's chart changes), domains
covering S⁴ × T² (manifest unions). That is the lemma list of §3 below.

## 1. State of the art (2026-10-07)

### 1.1 The construction exists and is verified numerically

F is T³-equivariant (the paper's T² plus the circle ζ ↦ e^{iθ}ζ, which is a holomorphic action on Y with two fixed
points on W = Y/T²). F = (Ψ_W ∘ π, φ): Ψ_W: W → S⁴ is an equivariant lift of an orbit-space map Λ: B = W/S¹ → S⁴/S¹,
and φ: Y → T² is a trivialisation. Everything is in `docs/s4t2_explicit_diffeo.md` §12 (authoritative), with the
code in `extra/s4t2_explicit/`. Verified: deg Ψ_W = 1.000000 with det > 0 at all 18.4·10⁶ grid nodes; det DF of one
sign in all charts; F⁻¹ to 10⁻¹³; J on the round S⁴ × T² with |N|/scale ≤ 10⁻¹⁰; J is T³-invariant. **This is
evidence, not proof** (§12.11(4) of that file states precisely what is and is not established).

### 1.2 The symbolic programme has started

**Update (EH-18):** programme §8 now gives a global construction **by exact one-dimensional
quadratures**, with an explicit inverse and the tensor formula (8.14)–(8.15). It replaces the fitted
lift and takes R=Id, so I3–I5 are absent from this new construction. This extends the original
expression class; a finite elementary expression has not been obtained. Its base and connection
blocks are implemented and tested, while a complete numerical evaluator of the quadrature tensor
is still to be implemented. The global statement and expression-class extension require human review.
The table below retains the status of the original fitted-lift architecture.

**Closed-form continuation (EH-23):** programme §9 eliminates both scalar profile
inverses, certifies the redesigned shells with contraction ratio <0.7, and supplies
an algebraic tensor on two whole pole neighbourhoods after an explicit torus gauge
change. These blocks and tensors are implemented and independently checked. The
global height inverse, shell inverse and circle transport/phase quadratures remain.
There is a proof excluding algebraic deck-equivariant heights, **not** a proof that
an elementary global tensor is impossible. The full closed-form goal remains open.

**Literature continuation (2026-10-08-EH-1):** programme §10 gives a finite elementary
cusp height with certified node margin >14/25, implemented in a separate base and
connection option. A different, elliptic height has a closed inverse in standard
special functions, including the complete cusp-to-middle collar and exact central
axes. An algebraic round-S⁴ fibration with separated critical values is also
implemented and checked. Finite phase-logarithm and projection operations are
proved and implemented, conditional on explicit derivative/nonvanishing bounds.
Their concrete bundle embedding and bounds, and the elliptic block's pole matching,
remain unresolved. There is no completed finite global lift or tensor.

**Finite global tensor (2026-10-08-EH-2):** programme §11 gives a separate,
elementary inverse-first construction. It supplies the full cusp inverse/collar,
open central-fibre charts, a smooth pole patch with a finite elementary inverse
in both directions (programme (11.12g)–(11.12j)),
a concrete C¹⁰ bundle projector with proved bound L=10¹¹⁴, and finite projection
and phase formulas with fixed counts N=4L+1, K=8·32(N+1)L+1. Together with the
finite forward height for the torus trivialization, (11.18)–(11.19) specify a
global elementary G_cl and J_cl, and prove invertibility, gluing, J²=−Id and N_J=0.
This is a new mathematical statement for human review. The counts have 115 and
232 digits: the full expression has not been numerically expanded or evaluated.
Block, overlap, central-fibre and reduced-count tensor checks pass; the legacy
safety net also stays green. **The tensor expression is closed, but §0's stronger
finite elementary inverse-pair goal remains open:** F_cl exists smoothly, but
the global height transition (11.3) may still require implicit inversion.
The pole patch no longer requires one. F_cl is not used to evaluate J_cl. This result
does not resolve the old evaluator's I3–I5/T14 flags.

`docs/s4t2_symbolic_programme.md` gives the block-by-block symbolic specification of F and G (§3) and the lemma
list with status (§6.1). Summary:

| lemma | statement | status |
|---|---|---|
| A | Y's charts, gluings, lattice coordinates, Φ, Θ, ψ̂, T³-actions | identities (paper's `verify_construction.py`; `sym/chart_identities.py`) |
| B | lattice system nonsingular: det = (\|t\|²+1)(log\|t\|·(\|t\|²−1) + 2\|t\|²) > 0 | manifest |
| C | cusp height μ̃ strictly monotone along fibres (dμ̃/dρ = −2 Var), range ℝ, closed-form on the central fibre | manifest |
| D | fibre-coordinate blend a₁ = (1−λ)μ̃ + λ(x₄^slice − ½) monotone, no `wrap` needed | manifest |
| E | node correction monotone (**I1**) | **certified, corrected rational bounds** (§4.2): common margin −a₁′/(2rad) ≥ 7/20; earlier displayed constants replaced |
| F | base profile \|T\| = r·k(r) strictly increasing (needs s′ ≤ 2) | proven (`sym/lemmas_1d.py` L1, L4) |
| G | radial profile ν strictly increasing | manifest (L3) |
| H | fibre twist R, moebius, P_final: closed inverses | identities (L5, L6) |
| I | exact linear models M₀ = 2πQ₀, M_∞ = πQ_∞ at the fixed points | identities (`symbolic_identities.py`) |
| J | linear-blend shells of Λ_mod nonsingular (**I2**) | **certified** (§4.3, I2a–n): σ_min > 2 / > 1, including the infinity profile transition (EH-11) |
| K | the lifts land in S⁴, cover Λ, are equivariant and linear in the balls | identities (L7) |
| L | **I3** E ≠ 0 off the tube; **I4** G_A ≠ 0 on the tube; **I5** \|Arg g\| < π on the shell | **open**; the proposed single-lift replacement also needs smooth division, not only matching linearisations (programme L8) |
| M | T²-trivialisation φ equivariant, invariant, closed inverse | identities (C6) |
| N | F∘G = id on each piece; pieces of G agree on overlaps | C8 verifies the closed linear-ball pieces; C9 gives the complete tensor there; other pieces and overlaps remain open, including T14 |
| O | DG invertible ⇒ J smooth, J² = −1, N = 0 ⇒ complex structure | C10 checks all six-dimensional N components with arbitrary symmetric second jets; the diffeomorphism hypothesis is essential |

Exact facts found on the way: φ(t) = 2|t|²/(1 − |t|²) at both cusps; the lattice determinant above; in node
variables the theta series has terms of total degree j²; its two-term truncation is accurate to fourth order; the
closed-form inverse of the smooth step; the central-fibre heights are logistic functions.

### 1.3 Candidate architecture (open design): "single compatible lift"

`docs/s4t2_symbolic_programme.md` §5–§7. With the current architecture (two lifts glued by a phase) the
conditions I4, I5 are intrinsically global. With a **single lift through one weight-1 function E whose zero set
is manifest and for which (v₁+iv₂)/E extends smoothly and never vanishes near the whole arc Γ_E = Λ⁻¹(meridian)**, lift B alone
is smooth off the fixed points (its only denominator vanishes on Γ_E; it is smooth across the antipodal arc), the
gauge, its normaliser ψ, the shell and both Chebyshev fits disappear, and what remains quantitative is short:
I1 and I2 (certified), E ≠ 0 off the arc in the cusps for |t| ≤ 0.3 (reducible to type (c)), and one
sign condition along the arc. The ingredients are in §7 of the programme file and in the task list below.
**Review correction:** matching transverse first derivatives does not imply this smooth-quotient condition;
programme L8 gives a symbolic counterexample and a sufficient exact-factorisation criterion. T4 must still
construct the quotient and extend the local modification globally.

## 2. Files, and how to run everything

| path | content |
|---|---|
| `docs/s4t2_explicit_diffeo.md` | construction log; **§12 authoritative numerical construction** (12.0 results, 12.1 files, 12.2 E, 12.3 Λ, 12.4 lift, 12.5–12.6 verification, 12.7 J and F⁻¹, 12.8 pitfalls, 12.11 J_round and the precise meaning of the result, 12.12–12.13 rigour analysis); **§13 the agreed proof scheme** |
| `docs/s4t2_symbolic_programme.md` | **the programme**: §2 one-variable lemmas, §3 symbolic specification of F and G, §4 fixed points (I1 certified, I2 structured), §5 global conditions, §6 status/lemma list, §7 the certificate goal and the single-lift architecture |
| `extra/s4t2_explicit/model.py, cusp.py, glue.py, hopf.py` | numpy reference implementation (independent of torch) |
| `extra/s4t2_explicit/middle.py, curve.py, lam.py` | superseded numpy modules, kept as a record (docstrings say so) |
| `extra/s4t2_explicit/tcore.py` | torch port of the building blocks (1:1 with the numpy reference; `check_tcore.py` compares) |
| `extra/s4t2_explicit/tlift.py` | E, Λ, the gauge, the lift Ψ_W |
| `extra/s4t2_explicit/tfull.py` | φ and F in the three holomorphic charts |
| `extra/s4t2_explicit/tinverse.py` | F⁻¹ (explicit chain + Newton polish) |
| `extra/s4t2_explicit/jround.py` | `J_round(q)`: J and DJ at any point of the round S⁴ × T² |
| `extra/s4t2_explicit/fit_arc.py → arc_fit.json`, `fit_psi.py → psi_fit.json` | the two fitted one-variable functions (deterministic; to be retired by the single-lift design) |
| `extra/s4t2_explicit/check_*.py`, `nijenhuis_round.py` | the numerical safety net (§12.0 table) |
| `extra/s4t2_explicit/symbolic_identities.py`, `check_nblind.py`, `toy_symbolic_blindness.py` | sympy: dμ̃/dρ = −2 Var and DΛ(b₀) = 2πQ₀; the N-check is blind to F (numerically and in a sympy toy) |
| `extra/s4t2_explicit/sym/lemmas_1d.py` | L1–L7: smooth step s′ ≤ 2 and its inverse, ν′ > 0, \|T\| increasing, P_final zones, moebius, the lift identities |
| `extra/s4t2_explicit/sym/chart_identities.py` | C1–C6: lattice (Cramer, determinant), actions, ψ̂, central-fibre heights, dμ̃/dρ, φ₀ |
| `extra/s4t2_explicit/sym/node_remainder.py` | N1–N5: φ(t) exactly, the node-variable series, μ̃ − m₀ structure, the I1 derivative identity, the I2 Jacobian identity, sampled margins |
| `extra/s4t2_explicit/sym/i2_margins.py` | I2 margins on the linear-blend shells (torch, 2 s) |
| `src/geometry/ambient.py`, `src/loss/ambient_loss.py` | vendored verbatim from `no_overlaps` @ `dcf2bab` (the ambient Nijenhuis tensor); frozen |

**Environment.** Python 3.10, torch 2.10.0 (CPU), sympy 1.14.0, numpy 2.2.6. Everything runs on a laptop.
**Run from the repository root.**

| command | time | expected |
|---|---|---|
| `python3 extra/s4t2_explicit/sym/lemmas_1d.py` | 20 s | `ALL PASS` |
| `python3 extra/s4t2_explicit/sym/chart_identities.py` | 10 s | `ALL PASS` |
| `python3 extra/s4t2_explicit/sym/node_remainder.py` | 40 s | `ALL PASS`, prints the margins of §4.2 |
| `python3 extra/s4t2_explicit/sym/i2_margins.py` | seconds | `ALL PASS`; proven bounds 2 / 1, measured 5.64 / 2.82 |
| `python3 extra/s4t2_explicit/symbolic_identities.py` | 10 s | passes |
| `python3 extra/s4t2_explicit/check_tcore.py`, `check_tlift.py`, `check_full.py` | seconds each | `ALL PASS` |
| `python3 extra/s4t2_explicit/check_step1.py`, `check_step2.py`, `check_glue.py`, `check_hopf.py` | ≈ 3 min total | PASS (numpy reference) |
| `python3 extra/s4t2_explicit/check_lift.py` | 2 min | PASS (§12.5 numbers) |
| `python3 extra/s4t2_explicit/check_degree.py` | 2.5 min (finest grid) | deg = 1.000000 |
| `python3 extra/s4t2_explicit/nijenhuis_round.py 2000` | 1.7 min | §12.7 numbers; `nijenhuis_round.py 60` is a 10 s smoke test |
| `python3 extra/s4t2_explicit/check_inverse.py 3000` | < 1 min | PASS, residual ≤ 10⁻¹³ (**run with 3000**: with small sample sizes, e.g. 200, torch 2.10 raises an internal forward-AD assertion in `lattice_coords`; pre-existing on the code branch, see T13) |
| `python3 extra/s4t2_explicit/check_jround.py 300` | 14 min | PASS (§12.11 table) |

`check_middle.py`, `check_curve.py`, `check_lam.py` test the superseded modules (`check_lam` item 3 aborts by
design, §11.8). Any change to `tcore.py` must be followed by `check_tcore.py`; any change to E, Λ or the gauge by
`fit_arc.py`, then `fit_psi.py`, then all checks (until the single-lift design retires the fits).

## 3. Task list

Each task names its inputs, its deliverable and its acceptance criterion. Claim a task by a ledger entry of type
`flag` saying "working on T<n>" (so that parallel work is visible), deliver it in one or more commits, each with a
ledger entry. Tasks marked **(paper)** can be done by a chat assistant from the documents alone (it writes the
argument and the script; the collaborator runs the script). Tasks marked **(code)** need the repository.

* **T1 (paper, small). Write the I2 certificate.** Input: programme §4.3, `sym/i2_margins.py`, `sym/node_remainder.py`
  N4. Deliverable: an explicit bound |Λ − M₀v_q| ≤ C_Λ|v_q|² and |DΛ − M₀| ≤ C′_Λ|v_q| on rad ≤ 0.008, for the
  closed-form chain h_B → moebius → P_final in the variables (t, a₁), with C_Λ ≤ 60, C′_Λ ≤ 120 (cusp 0; half at ∞),
  by second-derivative bounds of explicit functions. Acceptance: type (c) argument, constants ≥ measured (39.5, 79.0;
  19.7, 39.5) and below the thresholds; a sympy script verifying the second derivatives used.
  **Review note (2026-10-07-EH-7):** include the nonconstant base profile k on part of the infinity shell;
  the measured near-halving of the constants is not an exact scaling identity (programme §4.3).
  **Completed (2026-10-07-EH-11):** a direct derivative estimate for a common rational–trigonometric
  local map gives (C_Λ,C′_Λ)=(89,178) and (89/2,89), including the infinity profile. These replace the
  proposed smaller targets above and still prove σ_min > 2 and > 1. Exact identities and rational endpoint
  checks: `sym/i2_margins.py` I2a–n, ALL PASS; numerical formula/derivative agreement on 4000 points per cusp.
* **T2 (paper, medium). The aligned-phase lemma (I3 at the cusps).** Input: programme §5.2 (the exact ratio
  w_{n+1}/w_n = (c + |m|²|y|²t)/(m̄(|y|² + ct̄))), `cusp.py`/`tcore.py` (χ with inner 0.4, outer 0.9). Deliverable:
  (i) a choice of χ-supports such that on every overlap of two χ-supports |y₀|² > c|t| manifestly for |t| ≤ 0.3,
  proven from the series (type (c): μ̃ ≤ −δ forces Y ≥ Y_min(δ, |t|)); (ii) explicit bounds on the three phases
  giving |arg(w_{n+1}/w_n)| < π; (iii) the conclusion E ≠ 0 off the arc on |t| ≤ 0.3 at both cusps. Acceptance: a
  written lemma with constants; a script checking the series inequalities used; `check_step2.py` still passing with
  the new χ (code change, see T6).
* **T3 (paper, medium). The theta section with r-dependent modulus.** Input: programme §5.2 (f = e^{−πs²/κ}ϑ(−x +
  is/κ; i/κ)), §7 (iii). Deliverable: an explicit family E_mid = e^{2πix₁}·ϑ_{τ(r)} on 0.3 ≤ |t| ≤ 1/0.3 with the
  twisted periodicity f(s+1, x) = e^{2πix}f(s, x), one simple zero per cell at a prescribed, explicit point, and
  with the complex structure z_r = x + c(r)s matching the cusp structure z_E = log|y₀| + (i/2)arg t on the gluing
  annuli (closed-form c(r) from §3.3 of the programme and the slice formulas x₄ = log|ζ|/P + ξ₄). Acceptance:
  sympy verification of the periodicities and of the zero; a torch implementation agreeing with the formulas.
* **T4 (paper, hard). Compatibility of Λ with E along the arc.** Input: programme §5.2 (the sufficient linear-ball
  design and its limits), §7 (iv). Deliverable: (i) the induced complex structure on the torus T_r at the arc point from
  the closed-form chain, as an explicit Beltrami coefficient μ_Λ(r); (ii) either a choice of τ(r) in T3 making the
  theta's structure equal to it (then no twist), or the explicit anti-conformal twist (A, B)(λ) in the w-chart with
  |B| < |A| proven along the arc (one-variable, type (b)/(c)); (iii) the exact statement and proof of the
  "compatibility ⇒ smoothness of lift B across the arc" lemma (the Fact of `docs/s4t2_explicit_diffeo.md` §11.8),
  written out. Acceptance: a written lemma; a script evaluating μ_Λ(r) and |B/A|(r) on a fine grid with the margin.
  **Partial result (2026-10-07-EH-1):** (iii) corrected: exact smooth factorisation is sufficient, and
  first-order compatibility alone is false (programme L8, `sym/lemmas_1d.py`). Parts (i)–(ii), the required
  smooth quotient, and the global extension remain open.
* **T5 (code, medium). The single-lift implementation.** Input: T3, T4. Deliverable: `tlift.py` with E_mid the theta
  section, R from the explicit zero (no `fit_arc`), the compatibility twist, lift B alone, no gauge/ψ/shell; the
  log-ratio gluing of the cusp E to E_mid on the annuli with the windings computed symbolically. Acceptance:
  `check_lift`, `check_degree` (= 1.000000), `check_full`, `nijenhuis_round`, `check_inverse`, `check_jround` all
  pass; `check_tlift` updated; the smoothness test of §12.5 (jump/h independent of h) passes across the arc.
* **T6 (code, small). Narrowed χ-supports** (from T2) in `cusp.py`/`tcore.py`; rerun `check_step2`, `check_tcore`.
* **T7 (code, small). Remove `wrap` from the λ-blend** (`tcore.a1_c0`, `a2_ci`; programme §3.3) and verify that
  nothing changes (`check_hopf`, `check_tcore` against an updated numpy reference).
  **Implemented (2026-10-07-EH-2):** both numpy and torch use the unwrapped convex blend; C7 proves exact
  equality for the infinite series via discrete-Gaussian reflection and monotonicity. Regression checks cover
  10725 points per cusp, including half-cell and blend boundaries.
* **T8 (paper+code, medium). F∘G = id per piece, symbolically.** Input: programme §3 (every block's inverse).
  Deliverable: sympy scripts, one per piece of G, treating the certified-implicit inverses as symbols with their
  defining identities, proving F∘G = id and DF(G)·DG = 1, and the overlap agreement of the pieces through Φ, Θ and
  the toric chart changes. Acceptance: `ALL PASS`; a table in the programme file listing each piece, its domain and
  the script.
  **Partial result (2026-10-07-EH-3):** C8 verifies the closed linear-ball inverse compositions at both poles.
  The other pieces, derivative identities and overlap table remain open.
  **Further result (2026-10-07-EH-13):** C9 gives a complete symbolic coframe and tensor throughout
  both caps, with a nonzero explicit determinant. `J_linear_ball` evaluates it independently of the
  forward-map derivatives; checks cover 256 points per cap and finite differences at poles and axes.
  **Quadrature construction (EH-18):** programme §8.5 gives both inverse compositions globally,
  with phase, radial, chart and torus identities C11–C15. Its shell inverse is the contraction (8.1),
  ratio <17/25, and its base inverse includes the nonnodal central fibres. This does not repair the
  legacy `Finv` limitation T14 or certify the original fitted lift.
* **T9 (paper, small). The generic naturality lemma.** A sympy script with G = four (or six) undetermined
  functions showing N(A J₀ A⁻¹) ≡ 0 and (A J₀ A⁻¹)² ≡ −1 for A = DG, and a paragraph citing the textbook argument.
  (`toy_symbolic_blindness.py` is an instance.)
  **Completed (2026-10-07-EH-15):** C10 checks all 216 components in dimension six with an arbitrary
  symmetric second jet; programme §6.1 gives the coordinate-normalization and commuting-frame proof.
* **T10 (paper, small). Domain covering.** A written list of the pieces of G with their domains (explicit sets in
  the (Z, X, φ) coordinates of S⁴ × T² and in B_S), and the one-line proof that they cover; the same for the
  charts of Y used.
  **Completed for G_Q (EH-18):** programme §§8.3, 8.5–8.6 list the free/cap domains, source phase
  charts, holomorphic charts and round coordinate patches, with their covering and overlap arguments.
* **T11 (paper, medium). The certificate document.** Assemble, in a new file `docs/s4t2_certificate.md`, the
  definition of J piece by piece and the complete lemma list A–O with, for each lemma, the type (a)–(d), the
  argument or the script, the constants and the margin. The document must be readable in one sitting.
  **Quadrature certificate supplied in the existing programme §8 (EH-18).** No new file was created;
  it states the expanded expression class and separates the exact construction from implemented blocks.
* **T12 (review, ongoing). Adversarial review.** Try to break every type (c) argument: look for a missing case, a
  constant that is not justified, a series bound that does not cover all terms. Record findings as `flag` entries.

* **T13 (code, small). `check_inverse.py` at small sample sizes.** With n = 200 (not with 3000) torch 2.10 raises
  `INTERNAL ASSERT FAILED … has_fw_view()` in `tcore.lattice_coords` during the Newton polish (forward-mode AD
  through `.imag` of a view). Find the batch-layout trigger, make the script robust for any n, keep the results at
  3000 unchanged. Not a mathematical issue.
  **Implemented (2026-10-07-EH-6):** reproduced also with torch 2.12.1; contiguous copies of the lattice
  vectors before real/imaginary views remove the forward-AD layout assertion, without changing any formula.

* **T14 (code+paper). Inverse on the central fibres.** The review found that the polar inverse chain is singular
  at the exact poles and can fail on non-nodal central-fibre points, which uniform target sampling misses.
  **Pole caps implemented (2026-10-07-EH-3):** the exact linear-ball inverse is now used before polar inversion,
  with C8 identities, chart-coverage checks and pole J/DJ regressions. **Open:** extend a chart-based inverse
  across the central fibres outside those caps; review reproduction and limits are in 2026-10-07-EH-4.

* **T15 (closed-form continuation).** Programme §9 completes algebraic replacements
  of k and P_final, with new I2 bounds, and cancels the height from the pole tensor.
  **Open:** a finite elementary height/inverse resolving the central fibres; a finite
  inverse for the modified shells; a closed meridian lift and compatible real phase.
  The algebraic-height obstruction in §9.4 is narrower than these open questions.
  Retain the quadrature construction until a replacement has a full gluing proof.
  **Further progress (2026-10-08-EH-1):** §10.2–3 eliminates the infinite Gaussian
  series from a redesigned forward height and certifies the same local shells.
  §10.5 solves a different height and its entire collar inverse using Jacobi amplitude
  and standard elliptic F. This is not elementary, and appending the old node
  correction would destroy the closed inverse on its transition. §10.6 implements
  finite phase/projection operations and identifies the remaining inverse-first
  shell, concrete bundle embedding and derivative-bound work. The complete
  finite global expression and inverse pair have not been obtained.
  **Tensor expression closed (2026-10-08-EH-2):** §11 completes the inverse-first
  pole construction, the concrete bundle embedding, global derivative bounds,
  finite circle lift and phase gluing, and the global coefficient formula for J.
  It uses an elementary height inverse instead of elliptic functions. The closed
  formula has fixed but impractically large counts; reduced-count tests are not
  a certificate for smaller counts. **Still open:** a finite elementary expression
  for F_cl, a compact practical expression, and human review of the new global
  statement. Those are distinguished from the now supplied finite formula for J.

Dependencies: T5 needs T3, T4 (and T2 for the cusp part); T8 and T11 need T5; T1, T2, T3, T6, T7, T9, T10, T13
can start now, in parallel.

## 4. How to contribute

1. Read `AGENTS.md` (all of it; §8 is the campaign section), this file, and the programme file.
2. Synchronise: `git fetch origin && git switch agentic_ledger && git pull --ff-only origin agentic_ledger`.
3. Run the fast checks (table in §2) to know the starting state.
4. Work on one task. Write results **into the documents** (`docs/s4t2_symbolic_programme.md` for lemmas and
   designs, with a new subsection per result; `docs/s4t2_explicit_diffeo.md` only for corrections of the
   construction record) and **into scripts** (`extra/s4t2_explicit/sym/<topic>.py`, printing `ALL PASS`).
5. Label every claim (a)–(d). Record margins next to bounds. Never present a numerical observation as a proof.
6. One logical change per commit; ledger entry in the same commit (`LEDGER.md` template; `Where` = the campaign
   files; `Verification` = the scripts run and their results, or `not run`).
7. Keep the safety net green (AGENTS.md §8.3). Keep this file's task list current: mark a task done with the ID of
   the ledger entry that closed it.
8. Notation: follow the two notes files (Unicode in Markdown; ASCII in code). Points of Y in the three charts:
   `mid` (z, w, t), `c0` (x, y, w), `ci` (x′, y′, w′). The fixed points p₀ (t = 0) and p_∞ map to S (X = −1) and
   N (X = +1). B_S = {(v, X): |v| + X² = 1}, v = Hopf(Z) = (2Z₁Z̄₂, |Z₁|² − |Z₂|²).

## 5. Pitfalls (binding; details in `docs/s4t2_explicit_diffeo.md` §11.10 and §12.8)

* `torch.linalg.solve` and `linalg.det` give wrong second derivatives when the inner derivative is forward-mode
  (torch 2.10); use Cramer's rule or `inv`/`pinv`. `z.abs()**2` has a wrong second derivative at z = 0; use
  Re² + Im² (`tcore.abs2`).
* Never log-blend two sections across a common zero unless their linearisations are compatible (the §12.2
  discontinuity).
* The slice coordinate log|ζ|/P collapses the central fibre onto the node; only the theta height resolves it.
* A ψ̂-sum of local models y + c·x̄ has its zero arc at arg t ∈ (π/2, 3π/2); the sign of c is forced.
* The cusp slice formulas degenerate near |t| = 1 (|F| ~ 10¹³ at 0.97): never evaluate them beyond |t| = 0.6
  except at the plateau radii.
* The N-check cannot detect a wrong F (naturality); the cross-chart agreement of J and DJ can.

## 6. Status log of this file

* 2026-10-07 (TSG): campaign opened; files migrated from the code branch; lemma list A–O; tasks T1–T12 open.
* 2026-10-07 (EH/Codex, worktree review): T7 and T13 implemented; C8 adds the closed pole-cap inverse
  and symbolic pole values of J; I1 receives a complete rational certificate. L8 corrects the single-lift
  smoothness criterion and L9 shows exact orbit-map linearity is sufficient, not necessary. I2 remains open,
  including its infinity base-profile transition. New T14 records the remaining central-fibre inverse issue.
  See ledger entries 2026-10-07-EH-1–9; EH-9 records the complete regression results. Statement changes
  require human review. No commit or push requested.
* 2026-10-07 (EH/Codex, continuation): **T1 completed** with a rational I2 certificate, including the infinity
  profile transition; proven singular-value margins >2 / >1. **T9 completed** by C10's generic six-dimensional
  integrability identity and a non-integrable control. **T8 advanced**: C9 gives the complete symbolic tensor
  throughout both pole caps, with an independent evaluator and verified determinant. See EH-10–16 for proofs,
  numerical comparisons and test scope. I3–I5 and the remaining global inverse/domain/overlap work are open.
  The I2 status change needs human review; all changes remain uncommitted and unpushed.
* 2026-10-07 (EH/Codex, global continuation): programme §8 specifies **G_Q, F_Q and J_Q globally by
  exact quadratures**. The new base has a contraction inverse (<0.68); the connection has normalization
  >1/2; a two-chart integral gives the north mismatch a global real logarithm. This removes I3–I5 from
  the new map without asserting them for the original one. C11–C15 check the additional identities;
  the new block tests include 4560 inverse targets, exact axes and quadrature gauge changes. The
  complete quadrature evaluator and a finite elementary formula remain separate further tasks.
  See EH-17 onward for review qualifications and final regression results. No commit or push requested.
