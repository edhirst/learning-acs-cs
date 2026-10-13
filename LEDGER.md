# Ledger of changes to `manuscript_draft/`

This is an append-only record of every change made to the manuscript on this branch:
what changed, when, by whom, why, and how it was verified. Git keeps the exact
diffs and every past version (`git log`, `git show`, `git blame`). The ledger keeps
the reasons and the verification, which git does not record. The binding rules are
in `AGENTS.md` §5.

## Rules

1. One entry per commit that touches `manuscript_draft/`, plus one per `flag` (an
   issue found but not fixed). Entries are appended at the end, in chronological
   order.
2. **Never edit or delete an entry.** Corrections and reversals are new entries that
   cite the earlier entry's ID.
3. **Entry ID:** `YYYY-MM-DD-<initials>-<n>`.
   - The date is the day of the commit.
   - The initials are those of the collaborator who deployed the agent:
     EH = Edward Hirst, TSG = Tancredi Schettini Gherardini, TSRS = Tomás S. R. Silva,
     JS = Jakob Stein.
   - n = 1, 2, … counts that collaborator's entries on that day.

   IDs are therefore unique without coordination. The commit message starts with
   `[<ID>]`, so `git log --grep '<ID>'` finds the commit.
4. **Types:**
   - `edit`: correction or clarification; meaning unchanged.
   - `statement`: changes a mathematical statement; needs human review.
   - `addition`: new content.
   - `structure`: files, sections or numbering.
   - `revert`: undoes an earlier entry.
   - `flag`: an issue reported, nothing changed.
   - `infra`: the workflow files; humans only.

## Template

```markdown
### <ID> — <short title>
- **Date:** YYYY-MM-DD
- **Collaborator:** <name of the human who deployed the agent>
- **Agent:** <tool and model, e.g. "Claude Code (Claude Opus 5.5)", "Codex (<model>)", or "ChatGPT (<model>), via chat; applied by <collaborator>">
- **Type:** edit | statement | addition | structure | revert | flag | infra
- **Based on:** <short hash of the commit the work started from>
- **Where:** <file and location, e.g. "main.tex, proof of Lemma~\ref{lem:slabs}(3)">
- **Change:** <precisely what changed; for small edits quote old → new>
- **Why:** <why this is an improvement; for mathematical changes, the derivation or argument>
- **Verification:** <how it was checked; LaTeX build / verify_construction.py / sphere_products.py: pass | fail | not run>
- **Open points:** <anything left for the humans, or "none">
```

## Entries

### 2026-10-01-TSG-1 — Set-up of the workflow and baseline import of the manuscript
- **Date:** 2026-10-01
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Claude Code (Claude Opus 5.5)
- **Type:** infra
- **Based on:** `f94b92e` (the emptied branch)
- **Where:** new files `README.md`, `AGENTS.md`, `CLAUDE.md`, `LEDGER.md` and `.gitignore`.
  `manuscript_draft/` was imported **unchanged**, exactly as provided by the
  collaborators: `main.tex`, `references.bib`, `utphys.bst`, `figures/github.png`,
  `code/verify_construction.py` and `code/sphere_products.py`.
- **Change:** set up the agentic-ledger workflow. The manuscript was not modified. This
  commit is the **baseline** against which all later changes are measured.
- **Why:** requested by TSG, so that agents deployed by the collaborators can improve
  the manuscript in a controlled, recorded and reversible way.
- **Verification (baseline state of the manuscript):**
  - `latexmk -pdf -interaction=nonstopmode main.tex`: 34 pages, 0 errors, 0 undefined
    references or citations, 0 BibTeX warnings. The only warnings are 37 hyperref
    "Token not allowed in a PDF string" warnings, caused by mathematics in section
    titles. No overfull or underfull boxes.
  - `python3 code/verify_construction.py`: `all checks passed`, all 15 checks
    (sympy 1.14.0).
  - `python3 code/sphere_products.py`: `Dimension 8 check passed`.
  - Bibliography: 41 `\cite` keys, every key defined and every entry cited.
- **Open points:** two observations, deliberately left unchanged.
  1. The macro `\hflogo` (preamble) refers to `figures/huggingface-color.png`, which is
     not in `figures/`. This is harmless while the macro is unused, but the build
     would fail if it were used.
  2. `amsmath` and `multirow` are each loaded twice in the preamble. This is harmless.

### 2026-10-01-TSG-2 — Cover ChatGPT and other chat assistants without repository access
- **Date:** 2026-10-01
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Claude Code (Claude Opus 5.5)
- **Type:** infra
- **Based on:** `fd02c6d`
- **Where:** `AGENTS.md` (introduction, §1.5, new §7); `README.md` (the `AGENTS.md` row, a Codex
  note, and a new section "Using ChatGPT, or another chat assistant without access to the
  repository"); `LEDGER.md` (examples in the template's Agent field).
- **Change:** a chat assistant cannot run git, the build or the checks. So it now works only
  from the exact files it is given. It returns changes as unified diffs or exact OLD/NEW
  pairs, never as rewritten text. It shows its verification, and it drafts the ledger entry
  and the commit message. The collaborator applies the proposal, checks `git diff`, runs the
  build and both scripts, completes the entry and commits. The README also covers Codex,
  which works in the repository and reads `AGENTS.md`, and data settings before uploading
  the unpublished manuscript to a chat service. §1.5 now allows the AI service the
  collaborator has chosen.
- **Why:** requested by TSG, since some collaborators use ChatGPT.
- **Verification:** `manuscript_draft/` is untouched; this commit changes only `AGENTS.md`,
  `README.md` and `LEDGER.md`.
- **Open points:** none.

### 2026-10-01-TSG-3 — Lemma 2.6, proof: name the space obtained by attaching one slab
- **Date:** 2026-10-01
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Claude Code (Claude Opus 5.5)
- **Type:** edit
- **Based on:** `0205e49`
- **Where:** `main.tex`, proof of Lemma~\ref{lem:Xtop}, second bullet ("Mayer--Vietoris …").
- **Change:** `the boundary map $H_2(Z_{n+1})\to H_1(W_n\cap W_{n+1})$` →
  `the boundary map $H_2(Z_n\cup W_{n+1})\to H_1(W_n\cap W_{n+1})$`.
- **Why:** $Z_N=\bigcup_{|n|\le N}W_n(2)$ is symmetric. So attaching $W_{n+1}$ to $Z_n$
  gives $Z_n\cup W_{n+1}$, which is not $Z_{n+1}$: the latter also contains $W_{-n-1}$.
  The Mayer–Vietoris sequence of $Z_n\cup W_{n+1}=Z_n\cup_{W_n\cap W_{n+1}}W_{n+1}$
  reads $0=H_2(W_n\cap W_{n+1})\to H_2(Z_n)\to H_2(Z_n\cup W_{n+1})\xrightarrow{\partial}
  H_1(W_n\cap W_{n+1})\cong\Z\to H_1(Z_n)\oplus H_1(W_{n+1})=0$. This is exactly the
  statement made, for the corrected space. Here $W_{n+1}\cap Z_n=W_{n+1}\cap W_n$ by
  Lemma~\ref{lem:slabs}(2), and $H_1(Z_n)=0$ because $Z_n$ is simply connected. The meaning
  of the proof is unchanged.
- **Verification:** sequence re-derived as above. LaTeX build: pass (34 pages, 0 errors,
  0 undefined references, same 37 hyperref warnings as the baseline; the PDF text shows the
  new formula). verify_construction.py: pass. sphere_products.py: pass.
- **Open points:** none.

### 2026-10-01-TSG-4 — Remark 11.8(1): $Z_0\cong L_+^\times/2^{\Z}$ holds for a specific $\tau$
- **Date:** 2026-10-01
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Claude Code (Claude Opus 5.5)
- **Type:** edit
- **Based on:** `8f9c0e8`
- **Where:** `main.tex`, Remark~\ref{rem:Zq}(1) ("The case $q=0$").
- **Change:** "By Remark~\ref{rem:onebundle}, $Z_0\cong L_+^\times/2^{\Z}$, the quotient … by
  $v\mapsto2v$." → "By Remark~\ref{rem:onebundle}, for $\tau=2\pi i/\log2$ one has
  $Z_0\cong L_+^\times/2^{\Z}$, the quotient … by $v\mapsto2v$ (for $\tau=i$, replace $2$ by
  $e^{2\pi}$)."
- **Why:** $Z_0$ depends on $\tau$. By Remark~\ref{rem:onebundle}, $Z_0\cong L_+^\times/\lambda^{\Z}$
  with $\lambda=e^{2\pi i/\tau}$. The proof of Theorem C fixes "for instance $\tau=i$", which
  gives $\lambda=e^{2\pi}$, not $2$. The value $\lambda=2$ corresponds to
  $\tau=2\pi i/\log 2$. The conclusion $Z_0\cong S(L_+)\times S^1\cong S^4\times S^3\times S^1$ is
  unaffected, since it holds for every $\tau$. For any real $\lambda>1$, such as both values
  here, the explicit map $v\mapsto(v/|v|,\log|v|\bmod\log\lambda)$ is well defined.
- **Verification:** $\lambda=e^{2\pi i/\tau}$ evaluated for both values of $\tau$ (gives
  $2.0000000000$ and $535.4916555248=e^{2\pi}$). This agrees with
  $|\lambda|=\exp(2\pi\Ima\tau/|\tau|^2)$ from Remark~\ref{rem:onebundle}. LaTeX build: pass
  (34 pages, 0 errors, 0 undefined references, same 37 hyperref warnings). verify_construction.py:
  pass. sphere_products.py: pass.
- **Open points:** none.

### 2026-10-01-TSG-5 — Table 11.1: count products, not rows
- **Date:** 2026-10-01
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Claude Code (Claude Opus 5.5)
- **Type:** edit
- **Based on:** `abe16ce`
- **Where:** `main.tex`, proof of Corollary D ("The list.") and Remark~\ref{rem:wasopen}.
- **Change:** "they are exactly the rows of Table~\ref{tab:spheres}" → "they are exactly the
  products listed in Table~\ref{tab:spheres}". "the nine rows of Table~\ref{tab:spheres} not
  marked …" → "the nine products in the rows of Table~\ref{tab:spheres} not marked …".
- **Why:** the table groups several products per row. Counted from the LaTeX source, it has
  22 rows holding the 40 products. The rows not marked "classical" or "---" are 8, holding 9
  products: $S^6$, $S^4\times S^2$, $S^4\times T^2$, $S^6\times S^2$, $S^6\times T^2$,
  $S^4\times S^2\times S^2$, $S^4\times S^2\times T^2$, $S^4\times S^3\times S^1$,
  $S^4\times T^4$. Both numbers, 40 and nine, are correct counts of products, so only the noun
  was wrong. The nine also agree with the script, which lists 3 products of dimension 6 and 6
  of dimension 8 needing a 2026 input.
- **Verification:** rows and products counted by parsing the `tabular` in `main.tex`;
  `sphere_products.py` output compared. LaTeX build: pass (34 pages, 0 errors, 0 undefined
  references, same 37 hyperref warnings). verify_construction.py: pass. sphere_products.py: pass.
- **Open points:** none.

### 2026-10-01-TSG-6 — Full review of the manuscript: scope, outcome, open points
- **Date:** 2026-10-01
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Claude Code (Claude Opus 5.5)
- **Type:** flag
- **Based on:** `74362a6`
- **Where:** all of `main.tex`, plus `code/`.
- **Change:** none. This entry records the review behind TSG-3 to TSG-5 and lists what was
  deliberately left to the authors.
- **Scope and outcome:** every section was read, and every computation was re-derived by hand.
  - Lemmas 2.1–2.6: the determinant identities; the toric charts, dual bases, shift and log
    form; Hausdorffness; the slab bounds; the Mayer–Vietoris induction.
  - §3: the lattice determinant, the null conic, the periods and the no-relation argument,
    the change of basis at $t=1$.
  - §4: every step of the Mumford–Tate filling, including the bounds in Steps 5–6; the
    conjugation at $\infty$; both exact heights and the bounds at $r=10^{-2}$ and $r=1/16$.
  - §5–§7: the twist; the gluing; the $T^2$-action; $e_4^*$; the slice argument; the Wang
    sequences; $\pi_1(W)=1$; $\chi(W)=2$.
  - §8–§9: the $h$- and $s$-cobordism argument; Theorem B(1)–(5); the deformation region,
    critical point and counterexample.
  - §10–§12: the remarks; Proposition~\ref{prop:CE} and its lemmas; Theorem C; the
    obstruction lemma; Table 11.1 checked against the partitions and their complex models;
    the higher-dimensional lists against `sphere_products.py`; the literature arguments that
    the manuscript itself proves (LVM, twistor, homogeneous manifolds).
  - Appendix A matches the script's 15 checks item by item.

  **No mathematical error was found.** TSG-3 to TSG-5 fix three local imprecisions.
- **Verification:** LaTeX build, `verify_construction.py` and `sphere_products.py` all pass
  after every commit. The build matches the baseline: 34 pages, 0 errors, 0 undefined
  references, the same 37 hyperref warnings.
- **Open points (flags, not changed):**
  1. *Proof of Proposition~\ref{prop:Wstandard}, Steps 1–2.* Step 1 shows only that each node
     is of type $I_1^+$ or $I_1^-$. Step 2 then merges them by [Matsumoto85, Thm.~4.1]. That
     the two nodes have **opposite signs and the same vanishing cycle** $a=[e_1]$ follows
     from the trivial monodromy around $\partial D$: the local monodromies are the
     transvections $T$ and $T^{-1}$ along $a$. It also follows from the signature remark
     ($\sigma(W)=0$). Consider stating it explicitly before invoking Thm.~4.1. The reviewer
     could not consult [Matsumoto85] to confirm that theorem's exact hypotheses.
  2. *Notation clash.* In Lemma~\ref{lem:Q}, $\Phi$ and $\Psi$ denote the maps
     $\mathcal Q\to Z$ and $Z\to\mathcal Q$. Elsewhere, $\Phi$ is the logarithmic transform
     (Lemma~\ref{lem:Phi}, Construction~\ref{con:Y}, \S\ref{sec:topology}) and $\Psi$ is the
     map $M\to S^1\times S^1$ (Lemma~\ref{lem:regular}). Consider renaming the two maps of
     Lemma~\ref{lem:Q}.
  3. *PDF bookmarks.* The 37 baseline hyperref warnings come from mathematics in 10 section
     and subsection titles: 26 math shifts, 8 superscripts, 2 `\times`, 1 `\infty`. The
     bookmarks drop these symbols. Using `\texorpdfstring{$…$}{…}` in those 10 titles would
     remove all 37 warnings without changing the typeset text.
  4. *Pinpoint citations to check against the sources.* The reviewer could not consult these:
     - LLZ26: (4), (6), Lemma 3.1, Def.~3.2, Prop.~3.5, Lemma 4.1, Thm.~1.2, Cor.~1.3,
       §§1, 2.1, 4.1, 5.2;
     - AC26: Thm.~4.5, Prop.~4.6, Props.~9.23–9.24, Rem.~9.25;
     - GFL26: proof of Prop.~2.2;
     - Engel26: Rem.~8.5;
     - Matsumoto85: §§1, 2, 4, Thms.~4.1–4.2;
     - BosioMeersseman06: §12, Lemmas 0.9–0.10;
     - Calabi58: p.~408;
     - GompfStipsicz99: Ex.~8.4.7;
     - Sutherland65: Thm.~3.1.

     Per `AGENTS.md` §3.2 they were not changed.
  5. The open points of TSG-1 (the unused `\hflogo`; `amsmath` and `multirow` loaded twice)
     remain open.

### 2026-10-07-TSG-1 — Extend the workflow to the symbolic campaign
- **Date:** 2026-10-07
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Claude Code (Claude Fable 5.1)
- **Type:** infra
- **Based on:** `7235e94`
- **Where:** `AGENTS.md` (§1.1 scope; new §8 "The symbolic campaign"), `README.md` (table rows for the campaign
  files, a paragraph in the introduction, a section "The symbolic campaign", one line in "Rules in brief"),
  `.gitignore` (campaign scratch directory).
- **Change:** by TSG's explicit request, the branch becomes the shared workspace of the *symbolic campaign* (a
  rigorous, constructive, symbolic-elementary proof that the round S⁴ × T² carries the paper's complex structure,
  through an explicit diffeomorphism and its inverse). §8 of `AGENTS.md` extends the agents' scope to the campaign
  files and states the campaign's rules: no claim without a check (four claim types (a)–(d), each labelled); the
  numerical safety net stays green; margins recorded; the pitfalls are binding; ledger entry per commit;
  reviewability first. The manuscript rules (§0–§7) are unchanged; manuscript tasks and campaign tasks do not
  touch each other's files.
- **Why:** TSG asked to migrate the campaign here so that all collaborators and their agents can contribute to the
  goal; the previous scope rule would have forbidden agents to touch the campaign files.
- **Verification:** workflow files only; LaTeX build / `verify_construction.py` / `sphere_products.py`: not run
  (the manuscript and its scripts are untouched by this commit).
- **Open points:** none.

### 2026-10-07-TSG-2 — Open the symbolic campaign: migrate the construction and the symbolic programme
- **Date:** 2026-10-07
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Claude Code (Claude Fable 5.1)
- **Type:** addition
- **Based on:** `78999e4`
- **Where:** new files `SYMBOLIC_CAMPAIGN.md` (entry point: goal, state of the art, file map, how to run, task list
  T1–T13, contribution rules, pitfalls); `docs/s4t2_explicit_diffeo.md` (construction log, §12 authoritative
  numerical construction, §13 proof scheme); `docs/s4t2_symbolic_programme.md` (the symbolic programme, §§0–7);
  `extra/s4t2_explicit/` (numpy reference, torch implementation, numerical checks, fits, sympy certificates in
  `sym/`, and the package marker `extra/__init__.py`); `src/geometry/ambient.py`, `src/loss/ambient_loss.py` vendored verbatim from the code branch
  `no_overlaps` @ `dcf2bab` with stub `__init__.py` files (frozen, AGENTS.md §8.1).
- **Change:** migration of everything produced on the code branch between 2026-10-06 and 2026-10-07 for the
  explicit diffeomorphism F: Y → S⁴ × T² and the symbolic description of F and G = F⁻¹, unchanged in content.
  The two notes files and the code are copied byte for byte; the only additions are the campaign entry point and
  the two stub package files. Mathematical state recorded in `SYMBOLIC_CAMPAIGN.md` §1: the construction is
  verified numerically (deg 1, det ≠ 0, N ≈ 0, F⁻¹); lemmas A–O with status (B, C, D, G manifest; E certified with
  explicit constants; F proven; A, H, I, K, M identities; J reduced to a second-derivative bound; L = the global
  conditions I3–I5 open; N, O to be assembled); the single-compatible-lift architecture as the route to the goal.
- **Why:** TSG's request: a shared workspace for the collaborators and their agents to reach a rigorous,
  constructive, symbolic-elementary proof that is clear and easy to review.
- **Verification (run on this branch, in this worktree, 2026-10-07):** `sym/lemmas_1d.py`, `sym/chart_identities.py`,
  `sym/node_remainder.py`: ALL PASS; `sym/i2_margins.py`: certificate bounds 4.42 / 2.20 vs measured 5.64 / 2.82;
  `check_tcore.py`, `check_tlift.py`, `check_full.py`: ALL PASS; `nijenhuis_round.py 60`: max |N|/scale
  6.3·10⁻¹¹, cross-chart controls 10⁻¹³, non-integrable control detected; `check_inverse.py 3000`: PASS. The
  vendored modules are identical to the code branch (`diff -q`). LaTeX build / `verify_construction.py` /
  `sphere_products.py`: not run (manuscript untouched).
- **Open points:** T1–T13 of `SYMBOLIC_CAMPAIGN.md` §3; in particular `check_inverse.py` with a small sample size
  (200) hits a torch 2.10 forward-AD internal assertion (pre-existing on the code branch; T13).

### 2026-10-07-EH-1 — Correct the proposed smoothness criterion for the single lift
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** statement
- **Based on:** `d4aa9dc`
- **Where:** `docs/s4t2_explicit_diffeo.md` §11.8; `docs/s4t2_symbolic_programme.md` L8 and §7;
  `extra/s4t2_explicit/sym/lemmas_1d.py` L8; `SYMBOLIC_CAMPAIGN.md` status and T4.
- **Change:** withdraw the claim that matching transverse first derivatives suffices for a smooth lift.
  Supply the sufficient exact condition f = H E with H smooth and nonvanishing, and the extended formulas
  Z₁ = |H|E/√(2d), Z₂ = (H̄/|H|)√(d/2), d = |v|−v₃ > 0. No runtime lift is changed.
- **Why:** E(z)=z+z̄² and v=(2 Re(z)√(1−|z|²),2 Im(z)√(1−|z|²),2|z|²−1) have matching first jets.
  Nevertheless lift B's Z₂ has directional derivatives 0,i,−i along 1,i,1+i, and Z₁'s second derivatives
  violate the parallelogram identity. This is a counterexample to the criterion imported in TSG-2.
- **Verification:** `sym/lemmas_1d.py`: ALL PASS, including the factorisation identities and both
  counterexample derivative calculations. Manuscript `verify_construction.py` and `sphere_products.py`: pass.
  LaTeX build: not run successfully (`latexmk` unavailable; fallback `pdflatex` stops at missing `bbm.sty`).
- **Open points:** human review of this statement change. T4 still needs smooth division and a globally
  invertible extension of its proposed local twist; agreement of first derivatives does not close T4.

### 2026-10-07-EH-2 — Prove and implement the unwrapped fibre blend (T7)
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** addition
- **Based on:** `d4aa9dc`
- **Where:** `docs/s4t2_symbolic_programme.md` §3.3; `extra/s4t2_explicit/hopf.py`, `tcore.py`,
  `check_hopf.py`, `sym/chart_identities.py` C7; `SYMBOLIC_CAMPAIGN.md` T7.
- **Change:** replace `mn + lam*wrap(sigma-mn)` by `mn + lam*(sigma-mn)` in both cusp blends, in numpy and torch.
  Keep `wrap` for its other uses. Correct the clutching equality to equality modulo one and the infinity
  node-support estimate to include 1/√a₀.
- **Why:** completing the square gives a discrete Gaussian centred at σ=ρ/P−1/2. Its mean is strictly
  increasing and equals each half-integer centre by reflection; hence |mean−σ|<1/2. Pairing j with −j and
  with 1−j also proves the bound for the implemented symmetric finite truncation. The node correction's
  support has |t|≤e²/40<0.2, disjoint from the active blend |t|>0.3. Thus removal changes no exact value.
- **Verification:** C7: ALL PASS. `check_tcore.py`: ALL PASS. `check_hopf.py`: old/new error exactly zero
  on 10725 points per cusp (half-cell boundaries, translations and blend endpoints included); shift errors
  3.6e−15 / 8.0e−15. Numpy `check_step1`, `check_step2`, `check_glue`, `check_hopf`: completed with the
  expected small errors and positive sampled fibre derivatives. Manuscript scripts: pass. LaTeX: not run
  successfully (environment limitations in EH-1).
- **Open points:** the construction's other proof obligations remain open; no new files or restructuring.

### 2026-10-07-EH-3 — Implement the closed inverse in both fixed-point balls
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** addition
- **Based on:** `d4aa9dc`
- **Where:** `extra/s4t2_explicit/tinverse.py`, `check_inverse.py`, `check_jround.py`,
  `sym/chart_identities.py` C8; programme §3.4 and lemma N; construction log §12.11; campaign T8/T14.
- **Change:** invert the existing unitary linear lift directly on each target cap, recover C=φ+Hc̃,
  and undo the T² slice. Dispatch these caps before the polar inverse, merge chart indices correctly,
  and reject incomplete chart coverage. Add exact-pole, nearby-axis and cap-boundary checks. Record the
  elementary pole values of J and test those values independently of the inverse residual.
- **Why:** the original `Finv` returned `{}` at both exact poles, so `J_round` silently returned J=0 and
  residual zero there. The specified linear-ball inverse was documented but not implemented. The new
  inverse is closed and avoids Newton in these caps. C8 verifies its algebraic compositions.
- **Verification:** C8: ALL PASS. `check_inverse.py 3000`: PASS, max residual 8.4e−14; 96 added pole,
  near-pole and boundary points: complete coverage, max residual 2.0e−14. Sixteen pole/axis J checks:
  max residual 8.9e−16, |J²+P|≤5.7e−15, |N|/scale≤2.6e−15, finite-difference DJ relative error 3.7e−9;
  the independently derived south/north pole matrices agree. Manuscript scripts: pass. LaTeX: not run
  successfully (EH-1).
- **Open points:** inverse evaluation on central fibres outside these caps still fails in a targeted test
  (EH-4). This is a local improvement, not a claim of global certification. No forward-map formula changed.

### 2026-10-07-EH-4 — Flag inverse failures on the non-nodal central fibres
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** flag
- **Based on:** `d4aa9dc`
- **Where:** `extra/s4t2_explicit/tinverse.py`, central-fibre inputs outside the linear caps; campaign T14.
- **Change:** record the limitation found during T12 review; no inverse for this region is supplied here.
- **Why:** uniform target sampling does not test the central fibres, and source-chart evaluations of N
  do not test whether the target-to-source inverse handles them. The polar chain degenerates there.
- **Verification (numerical, eight points per cusp, torch 2.12.1):** for each C=CUSP0,CUSPINF, let p be an
  8×6 zero tensor. For r in [0.005,0.02,0.05,0.5], put x=√(r/C.a0),y=0 in one row and x=0,y=√r in the
  next row (both real, w=0). Set q=TI.CHART_F[key](p) and call TI.Finv(q). Both batches raise an SVD
  convergence error at batch element 2. Repeated after the pole fix and the contiguous-copy AD fix;
  the failure remains. This is distinct from the `has_fw_view()` assertion of T13.
  Manuscript scripts: pass; LaTeX: not run successfully (EH-1).
- **Open points:** T14 must use chart coordinates through the central fibres, with branch choices and
  overlap agreement certified. Do not describe `J_round` as a successfully checked evaluator at every point.

### 2026-10-07-EH-5 — Replace incomplete I1 estimates by a rational certificate
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** statement
- **Based on:** `d4aa9dc`
- **Where:** programme §4.2 and lemma E; `extra/s4t2_explicit/sym/node_remainder.py`; campaign lemma E.
- **Change:** replace the unspecified epsilon and incomplete tail/variance estimates by explicit common
  constants C=27/10 and C₂=4, proving −a₁′/(2rad)≥7/20 for 0<rad≤1/20 at both cusps. Verify the exact
  phi formula and both slice-determinant signs, and make a failed symbolic check exit unsuccessfully.
- **Why:** the displayed estimates imported in TSG-2 did not establish C≤2.63 or C₂≤2.8 at infinity;
  the N5 sample checks never tested that proof. The replacement uses e⁴<55 with a rational Taylor-tail
  bound, explicit positive geometric tails, X̂²Y≤4rad³/27, and a retained X̂ factor in the variance estimate.
  The scaling-flow version includes both axes without using log coordinates there.
- **Verification:** `sym/node_remainder.py`: ALL PASS, including all rational comparisons and N3c on
  the axes. Proven margin 0.35; measured minima on 200000 points per cusp 0.8521 / 0.8261. No numerical
  assertion or tolerance was weakened, and no construction formula changed. Manuscript scripts: pass.
  LaTeX: not run successfully (EH-1).
- **Open points:** human review of the changed quantitative statements. I2 remains open.

### 2026-10-07-EH-6 — Fix the small-batch forward-AD layout failure (T13)
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** edit
- **Based on:** `d4aa9dc`
- **Where:** `extra/s4t2_explicit/tcore.py`, `e_vectors`; campaign T13.
- **Change:** return contiguous owned copies of the lattice vectors before consumers take real/imaginary
  views. This changes storage only, preserving values and all mathematical derivatives.
- **Why:** `check_inverse.py 200` reproduced the documented `has_fw_view()` assertion even on torch 2.12.1.
  An isolated wrapper copying these vectors removed the failure; the same change is now in the implementation.
- **Verification:** after the edit, `check_tcore.py`: ALL PASS; `check_inverse.py 200`: PASS, max residual
  3.9e−14; `check_inverse.py 3000`: PASS, max residual 8.4e−14, unchanged from the earlier run. Both include
  the existing transported-J calculation and the new 96-point cap test. The pole J/DJ finite-difference
  check also passes with the copies. Manuscript scripts: pass. LaTeX: not run successfully (EH-1).
- **Open points:** the central-fibre inverse failure of EH-4 is separate and remains open.

### 2026-10-07-EH-7 — Correct the interpretation of the I2 shell measurements
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** statement
- **Based on:** `d4aa9dc`
- **Where:** programme §4.3 and lemma J; `extra/s4t2_explicit/sym/i2_margins.py`; campaign T1/lemma J.
- **Change:** label the sampled lower estimates as numerical evidence, withdraw the unsupported exact
  scaling description, and point out the base-profile transition in the infinity shell.
- **Why:** at the outer equator |s|=0.004e², so |t|=250e⁻²<e⁴ (the positive Taylor sum for e⁶ through
  degree 7 exceeds 250). Thus k is not yet on its constant plateau, and a rigorous infinity estimate
  cannot be obtained simply by halving the cusp-0 constants. The shell's global nondegeneracy remains open.
- **Verification:** `sym/i2_margins.py`, 4000 points per cusp: measured minimum singular values
  5.6446 / 2.8197, sample-based estimates 4.4162 / 2.1969, unchanged. The rational partial-sum assertion
  passes; at the specified single point k/KINF−1=4.871e−7. Manuscript scripts: pass. LaTeX: not run
  successfully (EH-1).
- **Open points:** human review of the corrected status; T1 still needs a bound over the whole shell.

### 2026-10-07-EH-8 — Distinguish a sufficient linear model from a necessary condition
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** statement
- **Based on:** `d4aa9dc`
- **Where:** programme §5.2; construction log §13.3(1); `sym/lemmas_1d.py` L9; campaign T4 input description.
- **Change:** withdraw the general assertion that smoothness forces an exactly linear orbit map near a fixed
  point. Retain the linear balls as the sufficient choice used by the actual implementation.
- **Why:** for w=Uq, r=|w|², the map Z=√scale·√(1+r²)w with
  X=±√(1−scale·r(1+r²)) is smooth near the node and induces the nonlinear orbit map
  v=(1+|v_q|²)Mv_q. Its positive radial factor preserves the required divisibility. Thus a generic
  problematic remainder is not an obstruction to every nonlinear model.
- **Verification:** `sym/lemmas_1d.py` L9: ALL PASS (Hopf homogeneity, norm identity, and the explicit
  height-radicand bound >9/10 for r≤1/100 and scale≤7). This is a local example only, not an implemented
  replacement for the global map. Manuscript scripts: pass. LaTeX: not run successfully (EH-1).
- **Open points:** human review of the statement correction. A nonlinear model may be useful for a future
  simplification, but no claim about its global extension or the removal of I2 is made here.

### 2026-10-07-EH-9 — Record the campaign review and complete regression results
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** addition
- **Based on:** `d4aa9dc`
- **Where:** `SYMBOLIC_CAMPAIGN.md` task/status summaries; this ledger; changes detailed in EH-1–8.
- **Change:** record the completed review and verification, distinguishing proven identities and bounds,
  numerical evidence, implementation repairs, and outstanding proof obligations. All changes remain in the
  working tree under the standing no-commit/no-push instruction (construction notes §13.5).
- **Why:** random testing missed both the exact-pole coverage defect and the remaining central-fibre issue.
  The new targeted checks supplement the original safety net; they do not certify a global diffeomorphism.
- **Verification:** synchronized `agentic_ledger` with origin before work. Python 3.10, sympy 1.14.0,
  numpy 2.2.6, torch 2.12.1 from the existing uv environment. Baseline and post-change runs completed for
  `check_tcore.py`, `check_tlift.py`, `check_lift.py`, `check_full.py`: ALL PASS. Post-change
  `check_degree.py 1121 128 128`: degree 1.000000, positive determinant at all 18366464 grid nodes.
  `nijenhuis_round.py 2000` after the storage fix: maximum normalized N 2.59e−10; the finite-difference
  derivative and cross-chart comparisons completed, and the non-integrable controls were detected.
  `check_inverse.py 3000` after the storage fix: PASS, maximum residual 8.4e−14; the 200-point run also
  passes, replacing the reproduced layout failure. The added 96 cap/boundary cases pass (EH-3).
  `check_jround.py 300`: ALL PASS, including 100 points on each chart switch, finite-difference DJ,
  torus invariance and infinity/middle agreement. This process loaded the cap/wrap changes before the
  final storage-only change; `check_jround.py 12` subsequently tested that final runtime: ALL PASS.
  The strengthened 16-point pole/axis check, including the symbolic pole-matrix comparison, was also run
  separately on the final code and passed. NumPy `check_step1.py`, `check_step2.py`, `check_glue.py`,
  `check_hopf.py`: baseline and post-change runs completed with the expected diagnostic values;
  the new 10725-point-per-cusp blend regression passed. `sym/lemmas_1d.py`, `sym/chart_identities.py`,
  `sym/node_remainder.py` on the final scripts: ALL PASS. `sym/i2_margins.py`: numerical measurements
  and the exact profile-transition comparison completed (EH-7); no global certificate is claimed.
  The unchanged `symbolic_identities.py` passed during baseline review. Both manuscript scripts ran
  before and after the work: `all checks passed` and `Dimension 8 check passed`. LaTeX: not run
  successfully (`latexmk` unavailable; the fallback `pdflatex` attempt stopped at missing `bbm.sty`).
  `git diff --check`: pass. No manuscript, frozen `src/`, protected workflow file or existing untracked
  collaborator file was modified; no file was added or removed from the repository.
- **Open points:** human review of EH-1, EH-5, EH-7 and EH-8. The global I2–I5 bounds, smooth quotient
  and global extension for the proposed single lift, remaining inverse/overlap identities, and T14 remain
  open. This audit does not establish correctness of every unproved claim or replace a manuscript review.

### 2026-10-07-EH-10 — Working on T1: certify the transition shells
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** flag
- **Based on:** `d4aa9dc`, with the uncommitted changes EH-1–9 retained.
- **Where:** campaign T1; programme §4.3; `extra/s4t2_explicit/sym/i2_margins.py`.
- **Change:** claim T1 for this continuation. Seek a direct derivative bound for a common rational–trigonometric
  local model, including the nonconstant infinity profile identified in EH-7.
- **Why:** I2 is an outstanding local certificate independent of the global lift redesign.
- **Verification:** synchronized with origin; existing changes retained. No new mathematical claim is made
  by this entry. LaTeX: not run (`latexmk` unavailable). Results will be appended after verification.
- **Open points:** derive and check the formula, its domain, the derivative estimates and their constants.

### 2026-10-07-EH-11 — Certify both I2 transition shells by a common local formula
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** statement
- **Based on:** `d4aa9dc`, continuing the uncommitted work recorded in EH-1–10.
- **Where:** programme §4.3 and status; campaign lemma J/T1; `sym/i2_margins.py` I2a–n.
- **Change:** change I2 from open to certified. Reduce the unblended local orbit maps to B at 0 and an
  orthogonal transform of B/2 at infinity, with the actual base-profile factor retained. Prove common
  normalized derivative bound 9|z|/2; after scaling, C′=178 / 89 and C=C′/2 suffice. The blended maps
  have minimum singular value >2 / >1 and are injective on their respective node balls.
- **Why:** the derivative of the rational–trigonometric quotient can be bounded directly, avoiding a long
  Hessian expansion. The infinity profile correction flagged in EH-7 has |1−β|<1/500000 and logarithmic
  derivative <1/25000, proven by positive Taylor bounds. The earlier target constants 60/120 and 30/60
  are replaced by larger proven constants that still meet the actual nondegeneracy requirement.
- **Verification:** `sym/i2_margins.py`: ALL PASS for the exact reduction identities, quotient derivative,
  rational endpoint comparisons, exponential/π bounds and the positive margins. Independent comparison
  with the unchanged forward implementation on 4000 shell points per cusp: value error ≤7.0e−16,
  first-derivative error ≤1.0e−15. Measured minimum singular values remain 5.6446 / 2.8197 beside proven
  bounds 2 / 1. Baseline `check_tcore.py`, `check_tlift.py`, `check_full.py` and the three other sym scripts:
  ALL PASS. Both manuscript scripts: pass. LaTeX: not run (`latexmk` unavailable). No forward-map or
  inverse implementation changed for this result; the new reduced formula is independently evaluated
  inside the certificate script. Existing assertions and tolerances were not weakened.
- **Open points:** human review of the new proof and status change. I3–I5 and the remaining inverse/gluing
  work remain open; local injectivity does not assert global injectivity or prove smoothness of the lift.

### 2026-10-07-EH-12 — Working on an explicit tensor formula in the pole caps
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** flag
- **Based on:** `d4aa9dc`, with EH-1–11 retained in the working tree.
- **Where:** campaign T8; programme §3.4; `sym/chart_identities.py`; `jround.py` and its check.
- **Change:** derive the complete local tensor J from an explicit coframe, starting with the closed cap
  inverse C8, and verify it independently of automatic differentiation of the forward map.
- **Why:** closed inverse coordinates should yield a directly evaluable symbolic complex structure on
  these domains, rather than only its values at the two poles.
- **Verification:** this entry records work in progress, not a new mathematical assertion.
- **Open points:** derive the coframe, determinant, tensor and torus-periodicity identities; compare with
  the existing evaluator, including axes and cap boundaries.

### 2026-10-07-EH-13 — Give the complete symbolic tensor on both pole caps
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** addition
- **Based on:** `d4aa9dc`, with the earlier entries' work retained.
- **Where:** programme §3.4; `sym/chart_identities.py` C9; `jround.py`; `check_jround.py`; campaign T8.
- **Change:** simplify the closed inverse to A=tW+Hf and w=σ(Wbar−Wt²+2Ht), write all three complex
  coframe forms and all their coefficients, including the explicit series for dH. Obtain the full tensor
  J_cap=adj(M)J_std M/det(M) with det(M)=−σ(1−|t|⁴)/(scale²a₀), and its ambient tangent extension.
  Implement this independently as `J_linear_ball`, retaining the existing `J_round` path as a reference.
- **Why:** this gives every tensor coefficient throughout the caps, not only at their poles. Holomorphic
  torus translations remove the angle branches; the evaluation needs neither automatic differentiation
  nor Newton. The symbolic expression uses the convergent infinite height series; numerical evaluation
  uses the existing finite sum j=−8,…,8, with no new truncation error certificate claimed.
- **Verification:** `sym/chart_identities.py` C9a–l: ALL PASS. The determinant follows from the verified
  unit-determinant complex normal block, triangular real shear, slice determinant and real linear
  Z-to-(x_s,y_s) map. Positivity of 1−|t|⁴ follows from the earlier cap bound |t|<0.016.
  `check_symbolic_linear_balls`: 256 points per cap, including poles, both source axes and boundaries;
  relative error against the independent forward-map tensor ≤5.8e−15, relative determinant error ≤9.6e−16.
  Finite differences on 12 interior/pole/axis points per cap: relative DJ error ≤1.2e−10 and normalized
  N ≤5.4e−11. J²=−P and tangency assertions pass. The new evaluator is separate; F, G and the original
  global evaluator are unchanged. Manuscript checks passed at the continuation baseline; LaTeX not run
  (`latexmk` unavailable). Final regression results will be appended with the continuation's closing entry.
- **Open points:** this is a complete local formula on the stated caps, not yet a global formula with all
  pieces certified. I3–I5, the remaining inverse domains/overlaps and T14 still require work.

### 2026-10-07-EH-14 — Working on T9: the six-dimensional naturality identity
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** flag
- **Based on:** `d4aa9dc`, continuing EH-1–13 in the working tree.
- **Where:** campaign T9; `sym/chart_identities.py`; programme lemma O.
- **Change:** claim T9: verify the generic integrability identity with arbitrary symmetric second derivatives
  of a six-dimensional coordinate map, and state the local-diffeomorphism hypothesis explicitly.
- **Why:** this isolates the automatic integrability of an actual coordinate pullback from the unresolved
  global diffeomorphism and smoothness obligations.
- **Verification:** work in progress; no new identity claimed in this entry.
- **Open points:** check all components, the derivative convention, and the coordinate-normalization argument.

### 2026-10-07-EH-15 — Verify the generic six-dimensional integrability identity (T9)
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** addition
- **Based on:** `d4aa9dc`, with EH-1–14 retained.
- **Where:** `sym/chart_identities.py` C10; programme §6.1; campaign T9 and lemma O.
- **Change:** verify N(AJ_std A⁻¹)=0 for A the Jacobian of a smooth local diffeomorphism, using an
  arbitrary symmetric second jet in six dimensions. Explain the linear normalization A=Id at a point
  and the equivalent commuting-frame proof. Keep the local-diffeomorphism hypothesis explicit.
- **Why:** this supplies the requested algebraic certificate for integrability independently of sampled
  N values, while leaving global smoothness, invertibility and overlap obligations separate.
- **Verification:** `sym/chart_identities.py`: ALL PASS. C10a checks the differentiated J² identity;
  C10b expands all 216 N components using 126 independent symmetric Hessian entries, all zero.
  C10c uses A=Id+q₂E₀₀, an invertible frame matrix that is not a coordinate Jacobian: the same component
  formula gives N¹_{0,3}=−1 at zero, so the cancellation is not vacuous. No bibliography data added.
  Both manuscript scripts pass; LaTeX not run (`latexmk` unavailable).
- **Open points:** none for the local identity. This does not certify that the campaign's global F is a
  diffeomorphism; the remaining lift and inverse/domain conditions are unchanged.

### 2026-10-07-EH-16 — Record final verification for the symbolic continuation
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** addition
- **Based on:** `d4aa9dc`, continuing the earlier uncommitted work.
- **Where:** campaign task/status summaries; results EH-11, EH-13 and EH-15; this ledger.
- **Change:** record T1 and T9 as completed, and the complete cap tensor as further progress on T8.
  Preserve all earlier work and the existing numerical construction.
- **Why:** the new I2 proof and tensor formulas need a clear distinction between exact certificates,
  numerical evaluation and the still-open global problem.
- **Verification:** synchronized origin before work. At the continuation baseline, both manuscript scripts
  passed; `sym/lemmas_1d.py`, `sym/chart_identities.py`, `sym/node_remainder.py`, `check_tcore.py`,
  `check_tlift.py`, `check_full.py`: ALL PASS. Final changed symbolic scripts `sym/chart_identities.py`
  and `sym/i2_margins.py`: ALL PASS. Final `check_inverse.py 200`: PASS, maximum residual 3.9e−14;
  96 cap/boundary points also pass with maximum residual 2.0e−14. Final `check_fixed_point_balls` and
  `check_symbolic_linear_balls`: PASS, with the margins recorded in EH-13 and the original 16-point
  pole/axis derivative check unchanged (relative error 3.7e−9). Both manuscript scripts rerun after
  the edits: pass. LaTeX not run (`latexmk` unavailable). `git diff --check`: pass.
  The full lift, finest degree, 2000-point N and 300-point J suites were not rerun in this continuation:
  their completed runs are recorded in EH-9, and F, G, tcore and the existing J_round function body
  were not changed here. The new independent cap evaluator has its own targeted value, determinant,
  tangency, J², finite-difference DJ and N checks. No assertion or tolerance was weakened.
- **Open points:** human review of the I2 statement/status change in EH-11. A global certified symbolic
  expression still requires I3–I5 and the remaining inverse/domain/overlap work, including T14.
  All changes remain uncommitted and unpushed under the standing instruction. No manuscript, frozen
  src file, protected workflow file or existing untracked collaborator file was modified.

### 2026-10-07-EH-17 — Investigate a global connection-transport lift
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** flag
- **Based on:** `d4aa9dc`, retaining the uncommitted results EH-1–16.
- **Where:** campaign global lift task; programme §5; proposed additions to `sym/chart_identities.py`.
- **Change:** investigate replacing the fitted phase-gluing construction by an explicitly specified
  circle-bundle connection, parallel transport over the certified base map, and a north-cap phase correction.
  This route would use exact one-dimensional quadratures; its admissibility as the requested symbolic
  expression has been asked as a preference, without stopping the derivation.
- **Why:** the current I3–I5 are still unproved. A connection-based lift may make nonvanishing and
  invertibility structural, while agreeing exactly with both existing linear models near the poles.
- **Verification:** synchronized origin, preserving all earlier work. Continuation baseline: both manuscript
  scripts pass; `sym/chart_identities.py`, `sym/i2_margins.py`, `check_tcore.py`, `check_full.py`: ALL PASS.
  LaTeX not run (`latexmk` unavailable). This entry makes no claim that the proposed lift is complete.
- **Open points:** connection descent, explicit chart transitions, transport/inverse formulas, north-cap
  correction, complete domain covering and independent verification. Existing runtime construction unchanged.

### 2026-10-07-EH-18 — Specify a global complex structure by exact quadratures
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** statement
- **Based on:** `d4aa9dc`, with EH-1–17 retained in the working tree.
- **Where:** programme §8 and status line; campaign task/status updates; notes §13.6;
  `tlift.py` connection/base additions; `tinverse.py` independent base inverse;
  `sym/chart_identities.py` C11–C15; `sym/i2_margins.py` I2o–p;
  `check_tlift.py`, `check_inverse.py` new block checks.
- **Change:** give the complete definition of a new parametrization G_Q, its inverse F_Q,
  and J_Q in coordinate and ambient matrix form. The base uses R=Id and the existing
  certified node modifications. A projected-defect contraction supplies the shell inverse;
  the central axes have a separate monotone scalar inverse. An explicit circle connection,
  scalar transport integrals, and a two-chart real logarithm of the north mismatch supply
  the global lift. Prove chart independence, both inverse compositions, domain covering,
  smooth agreement with the linear models on the reduced caps |Z|²≤1/100, and integrability.
  Correct the claim that the middle circle coordinates form a global product: their log-branch
  monodromy is essential. Existing source holomorphic charts and numerical evaluator are preserved.
- **Why:** this removes the unproved I3–I5 from the new construction, rather than assuming
  them or attempting to conclude global correctness from numerical Nijenhuis tests.
- **Verification:** derive the connection from the weighted Hopf form and reindex its deck
  average; S(H)>1/2. Derive the inverse contraction from the proven I2 derivative bounds;
  κ<17/25 by exact rational Machin bounds. The raw and modified extended maps agree outside
  the ball, giving identical ball images and a global inverse-domain proof. Derive the
  logarithm by integrating the closed form Im(conj(g)dg) on two stereographic planes; the
  two real primitives agree on their connected overlap. Derive both inverse compositions
  by transport reversal and phase cancellation. The tensor is the pullback through this
  smooth diffeomorphism; C10 supplies the generic six-dimensional N=0 identity.
  C11–C14 and I2a–p: ALL PASS at first completed runs; C15 and the final expanded overlap
  checks are included in the closing verification entry to follow. New numerical base test:
  4560 targets including both poles, shells and central axes; shell error 3.6e−17, full
  inverse residual 1.1e−14. Measured defect derivative/scale on 4000 points per shell:
  0.184870 / 0.184762, against the certified bound <0.68. Connection and quadrature
  measurements and their sample sizes are listed in programme §8.7. These are evidence
  for those blocks, not a numerical evaluation of the entire G_Q or J_Q.
- **Open points:** human review is required for this global statement and the expression-class
  extension. The exact expression permits one-dimensional quadratures, in addition to series
  and certified inverses; no finite elementary closed form is asserted. The optional question
  about accepting quadratures received no answer before this construction was written; this
  choice is stated explicitly, not treated as approval. The complete numerical quadrature
  evaluator is not implemented. Legacy Finv/J_round remain separate, with T14 and I3–I5
  still unresolved for that original parametrization. Full regression results follow when
  all launched suites finish. No manuscript or frozen source change, commit or push.

### 2026-10-07-EH-19 — Flag the legacy stereographic rationalization at its regular antipode
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** flag
- **Based on:** `d4aa9dc`, with the preceding working-tree additions.
- **Where:** legacy `tcore.moebius`; the new base test and `tlift.transport_stereographic`.
- **Change:** record that the legacy denominator
  (P₁²+P₂²+P₃²)/(1−P₀) is 0/0 at P=(1,0,0,0), although the stereographic map
  from (−1,0,0,0) is regular there. This point occurs on the infinity central fibre
  at the height half-cell. The new construction uses the direct denominator 1+P₀;
  its singular point is covered by the exact local node formula. The legacy function
  is left unchanged, so earlier evaluator limitations are not silently withdrawn.
- **Why:** random checks can miss this exact axis point; a complete domain specification
  must distinguish an actual chart singularity from a removable numerical one.
- **Verification:** the first 4560-target new base test produced nonfinite data at this
  point with the legacy rationalization. The direct formula passes the same test with
  residual 1.1e−14. Its inverse identities are L6; P₀=1 gives denominator 2 exactly.
- **Open points:** a future legacy-code repair should preserve stability near P₀=−1,
  for example by choosing the two algebraically identical denominator formulas on
  appropriate hemispheres. This flag is distinct from the broader legacy inverse issue T14.

### 2026-10-07-EH-20 — Clarify the domain and pairing of the legacy phase formula
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** statement
- **Based on:** `d4aa9dc`, with EH-1–19 retained.
- **Where:** programme §3.2, S¹-phase paragraph.
- **Change:** restrict the phase formula to Z≠0 and the indicated regular phase patches;
  state its pairing as sum Z_j conjugate(W_j), matching the implementation. Direct
  fixed-point inputs to the already verified cap inverse. Remove the verbatim duplicate
  fragment of that paragraph, which carried no additional content.
- **Why:** a phase is undefined at a fixed point; its pairing convention must not be
  left ambiguous, and a y-positive reference requires y≠0.
- **Verification:** for common nonzero Hopf data, C15b–c determine the same rank-one
  complex projector and squared norm, hence Z=λW and the numerator is λ|W|²≠0.
  The implemented phase in `tinverse.initial_guess` is exactly sum Z_j conjugate(W_j).
  No code behavior changed. C8 already verifies the cap inverse used at Z=0.
- **Open points:** this explicit domain clarification is a statement change for human
  review. It does not resolve the old lift's I3–I5 or Finv's T14 limitation.

### 2026-10-07-EH-21 — Audit the series domain and the new base's pole derivatives
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** statement
- **Based on:** `d4aa9dc`, with EH-1–20 retained.
- **Where:** programme §§3.2, 8.2, 8.6–8.7; new base functions in `tlift.py`;
  `check_inverse.check_transport_base`; `sym/i2_margins.py` comparisons.
- **Change:** replace the overbroad word "entire" for the height series with normal
  convergence on compact subsets of its cusp domain |xy|<1. Add tests of every
  second-derivative direction at both exact poles for the new base. Use the real
  quadratic Hopf polynomial, the inner polynomial radius and a cutoff guard, so
  nested derivatives do not pass through a norm at zero or a scaled complex-conjugate
  view. This changes only independent additions for the quadrature construction.
  Clarify the bundle coordinates used for the determinant argument and the T²
  invariance simplification of the final coframe formula.
- **Why:** correct values at a pole are insufficient; differentiability and the
  convergence domain must agree with the exact mathematical construction.
- **Verification:** the tail terms at index ±k are bounded by M^k r₀^{k(k−1)} on
  |x|²,|y|²≤M and |xy|≤r₀<1, with the same quadratic exponential decay for local
  derivatives. No claim of convergence beyond the cusp domain is needed or made.
  `check_transport_base`: all 16 second-derivative directions at each pole agree
  with an independent real polynomial model; 4560 value/inverse targets still pass
  with residual 1.1e−14. The first derivative-test attempt exposed a torch singleton
  layout assertion after a scaled conjugate; the real polynomial avoids that operation.
  `sym/i2_margins.py`: ALL PASS, now also comparing the new raw/blended base with
  the independent original shell on 4000 points per cusp: value errors ≤6.2e−16,
  derivative errors ≤4.7e−14. The mathematical formula and all tolerances are unchanged.
  Final `sym/chart_identities.py` C1–C15 and `check_tlift.py`: ALL PASS; the latter
  now includes both poles of the phase-log test and 64 real-overlap comparisons
  (error 1.2e−14, not merely equality of exponentials).
- **Open points:** human review of the convergence-domain clarification. The full
  quadrature evaluator remains unimplemented; the mathematical prescription is in §8.
  The legacy numerical J suite is still running as this entry is appended; its result
  will be recorded separately, rather than reported prematurely.

### 2026-10-07-EH-22 — Record the completed global-construction regression run
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** addition
- **Based on:** `d4aa9dc`, with all preceding working-tree results preserved.
- **Where:** verification of EH-18–21; programme §8; campaign status; this ledger.
- **Change:** close the verification record for the exact quadrature construction
  and its independent code additions. Preserve the distinction between the new
  analytic proof and the original numerical construction's regression suite.
- **Why:** a completed run, with exact scope and known limitations, is required
  before handing the new global statement to the collaborators for review.
- **Verification:** final `sym/lemmas_1d.py`, `sym/chart_identities.py` C1–C15,
  `sym/node_remainder.py`, `sym/i2_margins.py` I2a–p: ALL PASS. Final
  `check_tcore.py`, `check_tlift.py`, `check_full.py`, `check_lift.py`: ALL PASS.
  Numpy `check_step1.py`, `check_step2.py`, `check_glue.py`, `check_hopf.py`:
  completed successfully, with no failing checks. `check_inverse.py 3000`: PASS,
  maximum residual 8.4e−14; 96 additional cap/boundary targets pass to 2.0e−14.
  Its new 4560-target base check passes to 1.1e−14; the final separate rerun adds
  and passes all 32 pole second-derivative comparisons. The final I2 rerun compares
  the new base with the independent original shell: value ≤6.2e−16, derivative
  ≤4.7e−14 on 4000 points per cusp. `check_degree.py 1121 128 128`: degree
  1.000000, positive determinant at all 18,366,464 nodes. `nijenhuis_round.py 2000`:
  maximum normalized N 2.59e−10, with the independent chart and non-integrable controls.
  `check_jround.py 300`: ALL PASS, including finite-difference DJ, both chart switches,
  T³-invariance and the non-integrable controls; the 16 original pole/axis checks and
  both 256-point symbolic-cap comparisons also pass. The full legacy J run completed
  in 1553 seconds under concurrent verification load; it was not interrupted or weakened.
  Both manuscript scripts pass. LaTeX **not run** (`latexmk` unavailable), with no
  manuscript edits. `git diff --check`: pass. No assertion or tolerance was weakened.
- **Open points:** human review of EH-18's new global statement and expression-class
  extension, and of the domain clarifications in EH-20–21. J_Q is specified exactly
  by series, certified inverses and one-dimensional quadratures; a finite elementary
  expression is not established. Its base and connection are implemented; the whole
  quadrature map/tensor is not yet a numerical API. Legacy J_round retains its earlier
  I3–I5/T14 qualifications and the specific flag EH-19. All work remains uncommitted
  and unpushed. No manuscript, frozen src, protected workflow file, or pre-existing
  untracked collaborator file was changed; no new repository file was created.

### 2026-10-07-EH-23 — Eliminate scalar profile inverses and make both pole tensors algebraic
- **Date:** 2026-10-07
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** statement
- **Based on:** `d4aa9dc`, with EH-1–22 and the pre-existing untracked collaborator
  files preserved. Fetch/switch/fast-forward pull completed; origin and HEAD agree.
- **Where:** `docs/s4t2_symbolic_programme.md` §9 and status; campaign status/T15;
  `docs/s4t2_explicit_diffeo.md` §13.7; `extra/s4t2_explicit/{tlift,tinverse,jround}.py`,
  their three check scripts, and `sym/{lemmas_1d,chart_identities,i2_margins}.py`.
- **Change:** define a further global parametrization G_E using algebraic replacements
  for both scalar profiles and an explicit torus-coordinate translation. Supply and
  implement their closed inverses, a revised shell certificate, and the complete
  finite algebraic cap tensor. Retain all old defaults. Define the stricter requested
  expression class explicitly; distinguish §8's exact prescription from a finite
  closed form. Compress the height value using standard theta notation, and prove
  that an algebraic deck-equivariant height (or algebraic Y(height)) is impossible.
  Record the remaining global closed-form problem as open, with no claim of general
  impossibility. No existing theorem or manuscript passage is changed.
- **Why:** the user's request to remove the non-closed operations requires actual
  elimination, not new names for integrals or inverses. The balanced profile makes
  both local cusp corrections identical and admits a quadratic inverse. The radial
  compactification has an algebraic inverse. The torus translation cancels H,dH
  exactly on the caps, where even the previous C9 formula still contained them.
  The restricted negative result rules out finite rational height truncations but
  must not be overstated as an obstruction to elementary J itself.
- **Verification:** baseline and final manuscript scripts: `all checks passed` and
  `Dimension 8 check passed`. LaTeX **not run** (`latexmk` unavailable); no manuscript
  edit. All four symbolic suites: ALL PASS, now including L10–L12, C16–C17 and I2q–s.
  The new shell proof gives derivative constants 182/91, value constants 91/(91/2),
  singular-value bounds >1.9/>0.95 and contraction <0.7. Old bounds/assertions are
  retained for the old map. On 4000 independent samples per cusp, measured minimum
  singular values are 5.636589/2.815100 and contraction defects .183538/.185246;
  direct-composition/derivative discrepancies are at most 8.6e−15/2.4e−13.
  Both scalar inverse blocks have 4097 checks over 24 orders of magnitude, relative
  residuals 2.4e−15/5.8e−16. The revised base has 4560 random/shell/pole/axis checks
  (residual 1.3e−14) and all 32 pole second-jet checks. Its cusp profiles also pass
  all 14 derivatives through order three per central axis against independent Taylor
  polynomials. Squared moduli use real polynomials at zero, respecting §12.8's AD pitfall.
  `check_tlift.py`: ALL PASS, including 1024 new chart/middle/deck/branch comparisons
  per cusp for both the revised base and torus correction, exact inner cancellation,
  and 256 exact pole/axis values per cusp. `check_elementary_caps`: 256 points per
  cap, independent DG comparison ≤8.9e−16, determinant ≤8.0e−16, J²=−P, tangency,
  poles, axes and integer torus changes all pass. Separate finite-difference N checks
  on 44 points per cap are <1.3e−10 after division by max(1,cancellation scale);
  eight-point nonintegrable controls give median .26/.22. The original two 256-point
  cap suites also pass, including their independent DJ/N checks. C16's 25-point,
  60-digit theta/series comparison has error 9.3e−61; the exact identity is the
  termwise exponent/derivative calculation, with DLMF 20.2.3 checked at its source.
  The nonalgebraic-height proof uses polynomial growth versus q^n, with explicit
  constants and threshold in (9.16); it does not infer nonexistence from a CAS failure.
  Final `check_tcore.py`, `check_tlift.py`, `check_full.py`, `check_lift.py`: ALL PASS.
  Numpy step1/step2/glue/hopf suites complete successfully. `check_inverse.py 3000`
  passes (legacy residual 8.4e−14, 96 cap/boundary targets ≤2.0e−14); all three base
  variants pass. The final revised-base/axis-jet/cap checks were also run separately
  after the last code edits. Python syntax checks and `git diff --check` pass.
  The unchanged full degree/N/J regression results remain recorded in EH-22; those
  three expensive whole-domain suites were **not rerun in this continuation**.
- **Open points:** human review of the new parametrization and global substitutions.
  The global height inverse, modified-shell inverse and transport/phase quadratures
  remain: no full finite elementary global expression, and no proof that one is
  impossible, has been obtained. A complete numerical API for G_E/J_E is still absent;
  the new implemented tensor is restricted to the two pole neighbourhoods. The
  algebraic-height obstruction concerns Y=|zeta|² at fixed t, not all elementary
  functions, all choices of gauge, or all possible complex structures. This does not
  resolve the legacy I3–I5/T14 flags. No assertion/tolerance was weakened, no files
  were created/deleted, and no commit or push was made.

### 2026-10-08-EH-1 — Literature-led finite heights, elliptic inversion and finite bundle operations
- **Date:** 2026-10-08
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** statement
- **Based on:** `d4aa9dc`, retaining the existing uncommitted campaign work through EH-23.
- **Where:** `docs/s4t2_symbolic_programme.md` §10 and status;
  `docs/s4t2_explicit_diffeo.md` §13.8; `SYMBOLIC_CAMPAIGN.md` status/T15;
  `extra/s4t2_explicit/{cusp,tcore,tlift,tinverse,check_tcore,check_tlift,check_inverse}.py`;
  `extra/s4t2_explicit/sym/{lemmas_1d,node_remainder,chart_identities}.py`.
- **Change:** continued the requested literature investigation and closed-form campaign.
  Primary sources inspected are linked with their specific scope in §10.1: Di Scala–Kasuya–Zuddas
  on the Hopf-suspension construction and almost-complex embeddings; Zung and
  Bolsinov–Izosimov on singular-fibration gluing; Massamba–Thompson on universal
  connections; DLMF on theta, Jacobi amplitude, elliptic F and the beta/reflection formulas.
  No unpublished manuscript text was uploaded, and no manuscript/bibliography file changed.
  Added a finite elementary height H_f, evaluated by four seed values and an exact
  geometric tail on a reduced cell, including smooth axes and representative seams.
  Proved the cubic remainder bound 35 and corrected node slope >14/25. Implemented
  a separate finite-height option throughout the independent base/inverse/connection/
  torus-correction blocks; the original numerical construction keeps its defaults.
  A second height uses Jacobi amplitude with a closed inverse in standard elliptic F.
  Modulus interpolation closes its entire cusp-to-middle collar inverse. Its central
  extension is proved by a normally convergent independent arctangent series; its
  normalized cubic remainder is <38 and a hypothetical old node correction has
  slope >13/25. That correction is NOT appended to this closed inverse design.
  Added scalar reference evaluators, an algebraic S4->S2 target with two distinct
  critical values, and finite phase-logarithm/projection operations under explicit
  derivative/nonvanishing hypotheses. These are new mathematical statements for
  human review. The new section expressly does not claim a complete global closed form.
- **Why:** exact linear tails remove the infinite Gaussian sum in one forward design.
  Fourier comparison with DLMF's dn series identifies the alternative periodized
  arctangent height with Jacobi amplitude, whose standard inverse is elliptic F;
  changing the modulus instead of convexly blending heights preserves a closed inverse.
  The explicit S3 boost removes an unspecified map from the literature's target model.
  For a circle-valued g with ||Dg||<=L, N=floor(8L)+1 finite stereographic steps
  put every phase ratio in Re>7/8. Summed elementary arguments give a smooth real
  logarithm with exact overlap agreement, without its former quadrature. This result
  is conditional on a supplied g and proved L. For nearby rank-one projections,
  normalized projection is nonzero and reversing the chain is an exact inverse.
  These operations reduce the required bundle construction to a concrete embedding,
  pole compatibility and quantitative derivative bounds; they do not supply them.
- **Verification:** synchronised the branch with origin, preserving preexisting edits
  and untracked files. Baseline and final manuscript `verify_construction.py`:
  `all checks passed`; `sphere_products.py`: `Dimension 8 check passed`.
  LaTeX: **not run**, `latexmk` is unavailable; the manuscript is untouched.
  All four final scripts in `sym/`: `ALL PASS`, including new L13, N6/N7 and C18–C20.
  Identities and explicit constants are checked symbolically; normal convergence,
  monotonicity, Fourier identification and finite gluing are derived in §10.
  Finite height: 3000 independent sum comparisons per cusp, deck translations,
  both axes and third jets pass; relative errors <=1.3e-15. N6 samples 200000
  points per cusp: measured cubic constants .528244 / 32.280107 versus 35,
  corrected slope minima .998628 / .903375 versus 14/25.
  Elliptic height: 120 independent series/value/inverse tests at 80 digits, including
  P=-100 and translated cells, max error 1.5e-80; 20 exact-axis inverses and the node;
  72 complete collar inverses, error 1.1e-81; 36 seam jets through order three,
  error 1.1e-81. A failing inverse test exposed premature mpf rounding outside the
  guarded precision context; conversion was moved inside it. The original 1e-65
  assertion is unchanged. N7 samples 200000 points per cusp: cubic constants
  2.289646 / 32.441924 versus 38, hypothetical corrected slopes .982202 / .899635
  versus 13/25. These samples use an independent eight-term series with tail <2e-27.
  Finite-height base: 4560 inverse targets and all 32 pole second jets pass,
  residual <=1.3e-14; 2048 collar/base/connection/torus gluing cases per cusp and
  405 deck-seam jets per cusp/order pass. Algebraic target: 4096 sphere/inverse
  tests, error <=6.7e-15; tangent-rank margin min sigma2/sqrt(1-X^2)=1.175380>1;
  both exact critical points checked. Finite phase: 259 values including the poles,
  real error <=7.1e-15; 65 overlaps, real error <=3.6e-15 and derivative error
  <=3.4e-14. Finite projections: 128 chains, all steps <.05 apart, reversal/image/
  norm/equivariance errors <=1.2e-15.
  Legacy safety net: `check_tcore`, `check_tlift`, `check_full`, `check_lift` all
  `ALL PASS`; numpy `check_step1`, `check_step2`, `check_glue`, `check_hopf` complete
  successfully. `check_inverse.py 3000`: PASS, residual 8.4e-14; all independent
  base variants also pass. `check_degree.py 1121 128 128`: degree 1.000000,
  determinant positive at all 18366464 nodes. `nijenhuis_round.py 2000`: max
  normalized N 2.59e-10. `check_jround.py 300`: ALL PASS, including caps, finite
  differences, chart switches and nonintegrable controls (1453 seconds); uniform
  finite-difference normalized N 2.4e-6 at h=1e-6 versus the existing 1e-5 tolerance.
  Those global regressions test the legacy evaluator, not an unimplemented new J.
  Modified Python syntax and `git diff --check`: PASS. No assertion/tolerance was
  weakened; the existing `check_tcore` failure list now also has a final assertion.
- **Open points:** the requested full finite global expression has NOT been obtained.
  The finite elementary height still has an implicit transition inverse. The elliptic
  design closes the scalar block in an expanded special-function class, but its
  smooth pole matching and global circle lift remain to be constructed. Finite
  projection transport still needs the concrete bundle embedding, meridian bound
  and pole agreement, then a certified mismatch derivative bound for the finite
  logarithm. The specified algebraic target is not yet identified with W.
  No general impossibility theorem is claimed; these results do not exhaust all
  strategies. The new finite-height global substitution is a further quadrature
  construction, with no full numerical G_f/J_f API. Human review of the new
  statements remains required. Legacy I3–I5/T14 flags are not resolved by these
  results. No repository files were created/deleted, no commit or push was made,
  and the preexisting untracked analysis files and `extra/good_models/` were preserved.

### 2026-10-08-EH-2 — Finite elementary global tensor by an inverse-first construction
- **Date:** 2026-10-08
- **Collaborator:** Edward Hirst
- **Agent:** Codex (GPT-6)
- **Type:** statement
- **Based on:** `d4aa9dc`, retaining the uncommitted campaign work through EH-1.
- **Where:** `docs/s4t2_symbolic_programme.md` §11 and status;
  `docs/s4t2_explicit_diffeo.md` §13.9; `SYMBOLIC_CAMPAIGN.md` status/T15;
  `extra/s4t2_explicit/{cusp,tlift,check_tlift,check_inverse}.py`;
  `extra/s4t2_explicit/sym/{chart_identities,i2_margins}.py`.
- **Change:** supplied a finite elementary global coefficient prescription for J,
  through a new inverse-first G_cl. This changes the tensor-expression status
  recorded as open in EH-1, with the argument and limitations stated explicitly.
  The new scalar inverse uses square roots, exp/log and the existing flat step;
  f_a is bounded between 8/9 and 14, its full collar is monotone, and every deck
  seam jet agrees. Smooth open axis formulas, including neighbouring height cells,
  replace polar value formulas. A closed inverse pole replacement has derivative
  defect at most 13200·10^-8<1/1000. The ten-component circle-bundle embedding has
  norm squared between 1/2 and 2, with exact deck/branch compatibility and exact
  homogeneous pole models. An explicit chain-rule certificate gives projector
  bound L=10^114 on the specified slab. Fixed counts N=4L+1,
  L_g=32(N+1)L, K=8L_g+1 make all normalized projections and finite phase arguments
  nonvanishing. The projection phase derivative accumulates linearly in N.
  The resulting circle map matches the closed inverse pole maps on open overlaps.
  Restore T2 using the already verified finite forward height H_f, then use
  (11.19). The seven-component first-jet arithmetic rules give a coefficient
  prescription without an unevaluated differentiation operator. The proof supplies
  global bijectivity, smoothness, invertible differential, chart agreement,
  J^2=-Id and N_J=0. This is a new mathematical statement requiring human review.
- **Why:** a closed pullback tensor needs the inverse parametrization and its first
  derivatives, but does not need to evaluate the forward inverse. Designing the
  height and pole smoothing in that direction removes the former scalar and shell
  solves from J. Finite products define exact circle-bundle isomorphisms, rather
  than limits or approximations to transport. A concrete embedding and explicit
  derivative bounds remove the remaining unspecified inputs in EH-1's finite
  bundle operations. The torus height can be chosen independently of the base
  height, so its four-seed forward formula finishes the pullback. This produces
  the same holomorphic Y with different coordinates, without changing the legacy
  fitted evaluator or introducing a preprint dependency.
- **Verification:** synchronised with origin on `agentic_ledger`; no remote change.
  Baseline and final manuscript scripts: `all checks passed` and
  `Dimension 8 check passed`. LaTeX **not run** (`latexmk` unavailable); no manuscript
  changes. All four symbolic suites: ALL PASS, including C21-C24 and I2t-v.
  C21/C22 check reindexing, branch phases, normalization, stereographic/Hopf
  identities and the horizontal inverse norm. C23 checks the quadratic height,
  the normalized-projection phase derivative with arbitrary two-line tangent jets,
  132 inverse/deck/collar cases, 108 seam jets (error 6.4e-86), and 12 open-axis
  inverse comparisons plus 12 seam jets (2.6e-86). C24 verifies the finite
  first-jet rules and elementary argument chart. I2t-v check the rational endpoint
  comparisons of the written derivative certificates. The f_a samples lie in
  [1,4.536928] versus the proved [8/9,14]. Sixty-four independent pole forward
  compositions have relative error 2.1e-63; measured value/derivative constants
  2.552/11.348 versus proved 200/5000; blend defect 3.2e-7 versus proved 1.32e-4.
  The torch embedding has 4096 tests for each old height (error <=4.7e-14), and
  per height/cusp 1024 chart comparisons, 256 exact-axis tangent comparisons
  (<=8.6e-10), and 128 exact homogeneous models. The elementary scalar code has
  70 cusp/collar branch/deck/circle/independent-torus checks (7.4e-56), eight exact
  axis representatives, 96 raw-base round trips (1.4e-47), and projector residual
  2.0e-56. Sixteen projector derivative samples give Frobenius norm <=32.8697,
  beside the deliberately coarse proved operator bound 1e114. Twelve finite
  phase-family real-overlap checks give error 2.5e-54 despite multiple phase turns.
  With explicitly reduced counts N=512,K=32, both full cap/slab overlaps pass
  (3.1e-60/1.2e-54); sampled projection norms are >=.05504 and phase increment
  real parts >=.91584. These are tests of different finite instances, not numerical
  evaluation of the certified counts; .05504 does not meet the latter's >1/2
  bound. A preliminary N=64 instance actually vanished, so no smaller global
  count is claimed. Six pole/round-axis tensors and a free N=512 tensor pass
  nonzero-DG and J^2 checks; the free determinant is 2217.99. Independent sphere
  charts and integer torus shifts agree. Both a=1/2 central-circle tensors include
  transverse derivatives and have |det DG|=2.67101/.0900846. Independent finite-
  difference N at two caps and a free point falls by a factor of four when h is
  halved from 2e-7 to 1e-7; maxima 4.73e-10/1.18e-10. Nonintegrable controls
  are .566/.0255. The initial single-complex-line shear was integrable and was
  replaced by a sphere/torus shear without weakening the assertion.
  The overlap checks found a genuine log-branch metadata error: the supplied
  height described a principal-log cusp but conversion used an unreduced log.
  Principal-log normalization preserving (z,w) fixes it; the same assertions pass.
  Guard precision was increased for a cancellation-prone independent north-pole
  check, retaining its tolerance. Exact endpoint fixtures and rounding of returned
  guard digits avoid treating rounding beyond a closed domain as a counterexample.
  The complete legacy safety net was rerun: check_tcore, check_tlift, check_full,
  check_lift, all four numpy checks, check_inverse 3000 (residual 8.4e-14),
  check_degree 1121 128 128 (degree 1.000000, all 18,366,464 determinants positive),
  nijenhuis_round 2000 (maximum normalized N 2.59e-10), and check_jround 300
  (ALL PASS, 1509 seconds). These protect the old evaluator, not the unevaluated
  huge new formula. Python syntax and git diff --check pass. No assertion or
  tolerance in the existing safety net was weakened.
- **Open points:** human review of the new global mathematical statement. The finite
  expression is not compact or computationally practical: N and K have 115 and
  232 digits, and their full products have not been expanded or numerically run.
  The scalar implementation has a work guard; the numerical derivative helper
  is restricted to interiors of source charts and the tested reduced-count maps.
  This closes the finite elementary **tensor-expression** question by supplying
  G_cl and its finite coefficient rules. It does not close §0's stronger finite
  elementary **inverse-pair** goal: F_cl exists smoothly, but the inverses of the
  transition height and inverse pole patch may remain implicit. They do not enter
  J_cl. No general impossibility theorem or globally valid small count is claimed.
  The old fitted evaluator's I3-I5/T14 flags remain open. No manuscript, frozen
  module, instruction or bibliography file was modified; no repository files were
  created/deleted. No commit or push was made, and the preexisting untracked
  analysis documents and extra/good_models/ were preserved.

### 2026-10-09-EH-1 — Full write-up of the analytic construction, with its environment and reference outputs
- **Date:** 2026-10-09
- **Collaborator:** Edward Hirst
- **Agent:** Claude Code (Claude Opus 5.5)
- **Type:** addition
- **Based on:** `6e31eb3`
- **Where:** new file `docs/s4t2_analytic_construction.md`; new folder `extra/s4t2_analytic/` (`README.md`,
  `run_all.sh`, `export_table.py`, `environment/README.md`, `environment/requirements.txt`,
  `data/verify_construction.txt`, `data/sphere_products.txt`, `data/sphere_products.csv`,
  `data/.gitattributes`). No existing file changed apart from this ledger. `manuscript_draft/` is untouched.
- **Change:**
  - `docs/s4t2_analytic_construction.md`: a Markdown write-up of the closed-form construction and of all its
    proofs, following `main.tex` at `74362a6`, the last commit that changed it. It covers Main Theorem,
    Theorem B, the remarks, Theorem C, Corollary D with Table 11.1, the literature summary, a status section
    and a check-by-check verification map. Numbering follows the compiled manuscript. The manuscript stays
    authoritative.
  - `extra/s4t2_analytic/`: a pinned environment (sympy 1.14.0, mpmath 1.3.0), a runner, and the reference
    outputs of the two manuscript scripts. It also has a CSV export of the products-of-spheres
    classification: 159 products of dimensions 2–12, with almost complexity, status, minimal inputs and a
    witness decomposition.
  - The scripts are **not copied**. Everything calls `manuscript_draft/code/`, so the scripts of §3.3 keep a
    single source.
- **Why:** explicit request by EH: put the analytic construction (scripts, data, environment, and a full
  Markdown write-up) into this repository, following `AGENTS.md` as closely as possible.
  - These paths are outside the scopes of §1.1 and §8.1. They were created on that explicit request (§1.3),
    following the campaign's split between `docs/` and `extra/`.
  - `README.md` was not modified (§1.2). A suggested row for it was given to EH.
- **Verification:**
  - Manuscript scripts, before and after: `verify_construction.py`: `all checks passed`;
    `sphere_products.py`: `Dimension 8 check passed`.
  - LaTeX: **not run** as specified (`latexmk` unavailable). To read the result numbers, a scratch copy
    outside the repository was compiled with `pdflatex`/`bibtex`, with the packages missing locally removed
    from the preamble (bbm, mathdots, mathrsfs, multirow, wrapfig, tikz-cd, comment, appendix, algorithm,
    algpseudocode, blkarray). It gave 34 pages, 0 errors and 0 undefined references. The numbers were taken
    from its `main.aux`. This is not the build of §3.3.
  - `run_all.sh` was run from the repository root and from another directory.
    - The outputs are identical across runs, and identical to a fresh venv built from
      `environment/requirements.txt` (Python 3.10.13).
    - No `__pycache__` is left in `manuscript_draft/`.
  - Fidelity: an independent agent compared the write-up with `main.tex` line by line, covering every
    formula and constant, the result numbers and the mapping to the script's checks.
    - No formula differed.
    - Its 16 points were all fixed: overclaims of what the script checks (positivity in check 1 and
      monotonicity in check 11 are spot-checked and argued by hand), two misplaced check tags, and omitted
      hedges and premises (Remark 11.12 "to our knowledge"; the twistor premise b₁(Z) = 2).
  - `git diff --cached --check`: pass.
    - The CSV is written with `\n` line endings.
    - `data/.gitattributes` exempts the two `.txt` reference outputs from the whitespace check. They are stored
      exactly as `sphere_products.py` prints them, including trailing spaces after empty lists. The script
      itself is not changed.
- **Open points:**
  1. `README.md` (humans only) does not yet list the new files. A suggested row:
     `docs/s4t2_analytic_construction.md`, `extra/s4t2_analytic/`: full Markdown write-up of the analytic
     construction; environment, runner and reference outputs of the manuscript scripts.
  2. The write-up is derived from `main.tex`. When the manuscript changes, it should be updated in a later
     entry. Its header records the base commit.

### 2026-10-11-TSG-1 — Recover the forward map on raw pole caps and locate a radial obstruction
- **Date:** 2026-10-11
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6)
- **Type:** addition
- **Based on:** `338d95c` (the collaborator's pulled branch)
- **Where:** `docs/s4t2_symbolic_programme.md` §11.3 and §11.8;
  `extra/s4t2_explicit/sym/chart_identities.py` C25. No manuscript or construction file was changed.
- **Change:** gave a finite elementary formula for the forward inverse of the inverse-first pole lift
  on both unblended caps, with the explicit projection margin 525/529. Expanded the unchanged south
  raw base for (11.7) through the first radial obstruction and proved that its spinor radius is not C⁶ at the
  pole. This excludes a smooth diffeomorphic lift obtained only by changing phase over that raw base;
  it does not exclude a redesigned base or an elementary inverse of the existing blend. Added exact
  symbolic checks C25a–d and 18 independent numerical compositions C25e.
- **Why:** the stronger finite elementary inverse-pair goal remains open. The cap formula supplies a
  genuine forward piece of it, while the obstruction identifies why a phase-only removal of the
  pole blend cannot work. The formula is derived by reversing two rank-one projectors with a positive
  real scalar. For the obstruction, the exact norm identity for B and the profile/height jets give
  (11.12b). On a real spinor plane, its homogeneous sixth-order term is nonpolynomial: the numerator
  does not vanish at the complex root x=iy of x²+y², with coefficient
  η<1+55/576−48<−46. A C⁶ lift would make this norm C⁶, contradicting its Taylor expansion.
- **Verification:** two agents independently re-derived the south norm expansion and reviewed the
  local/global scope. `sym/chart_identities.py`: `ALL PASS`, including C25a–d and 18 C25e compositions
  with maximum relative error 6.3e−68 (numerical evidence only). All four campaign symbolic suites:
  `ALL PASS`. `verify_construction.py`: `all checks passed`; `sphere_products.py`:
  `Dimension 8 check passed`. `git diff --check`: pass. LaTeX: **not run** (`latexmk` unavailable in the
  active environment). The full legacy numerical safety net was not run because no construction block
  changed. Git fetch/pull and push were left to the collaborator as requested.
- **Open points:** the blend between |z_L|²=ε and 2ε still has no finite elementary forward inverse;
  the cap formula is a global F_cl component only on the (11.17) cap/overlap. Human review of the
  sixth-order obstruction is welcome.

### 2026-10-11-TSG-2 — Repair the source-modulus derivative certificate
- **Date:** 2026-10-11
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6)
- **Type:** statement
- **Based on:** `84d3d62`
- **Where:** `docs/s4t2_symbolic_programme.md` §11.4, steps 2–4;
  `extra/s4t2_explicit/sym/i2_margins.py` I2v. No manuscript or construction file was changed.
- **Change:** corrected a false intermediate derivative claim in the global projector certificate.
  The quadratic root Y₀=e^(2ρ₀) from (11.5) differs in the collar from the physical source modulus
  Y_*=e^(2ρ)=|y|². The source bounds and D_transverse ρ now use Y_*; the derivative of ρ₀ with
  respect to q continues to use Y₀. The core and collar triangular inverses are bounded separately.
  Added the explicit pure-middle bound ||Dt||<2·10¹⁰δ⁻⁴ needed for the displayed ||Dx₄|| bound.
  The global conclusion L=10¹¹⁴ and all fixed counts remain unchanged.
- **Why:** the old inequality |D_transverse ρ|≤1/(4Y₀) is not valid in the collar. For the reduced
  representative, both terms blended in (11.5) lie in [P,0], so q≤Y_*≤1. In the collar,
  D_transverse ρ=p_transverse/(4Y_*r), and q>1/3025 yields
  |D_transverse ρ|≤3025/4<1000. In the core Y_*=Y₀, and the displayed derivative bounds give
  |D_transverse a|<2·10⁵ and a_m≥1/56. In the collar they give
  |D_transverse a|<2·10⁶ and a_m≥1/2000. Both inverse norms remain <10¹⁶.
  The profile inverse and ratio bounds in step 1 give ||Dt||<128·100·10⁶δ⁻⁴<2·10¹⁰δ⁻⁴
  in the pure middle, completing the final inequality.
- **Verification:** two independent agents re-derived the core/collar and pure-middle bounds.
  I2v checks the physical-modulus derivative identities symbolically and all new constant
  comparisons: `ALL PASS`. All four campaign symbolic suites: `ALL PASS`.
  `verify_construction.py`: `all checks passed`; `sphere_products.py`:
  `Dimension 8 check passed`. `git diff --check`: pass. LaTeX: **not run** (`latexmk` unavailable).
  The full legacy numerical safety net was not run because no construction block changed.
  Git fetch/pull and push were left to the collaborator as requested.
- **Open points:** this correction of an intermediate mathematical claim merits human review.
  The finite elementary forward inverse of the pole blend remains open.

### 2026-10-11-TSG-3 — Reduce the blended pole shell to an exact Hopf-orbit equation
- **Date:** 2026-10-11
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6)
- **Type:** addition
- **Based on:** fd04e21
- **Where:** docs/s4t2_symbolic_programme.md §11.3 and §11.8;
  extra/s4t2_explicit/sym/chart_identities.py C26. No manuscript or construction file was changed.
- **Change:** expressed the forward problem for the existing blended pole shell as one finite elementary
  equation in the three real components of the raw source Hopf orbit. Once that orbit is known,
  an explicit positive rank-one projection recovers the spinor and its phase. The equation
  has a unique solution in the certified shell cap for data in its image. Added exact symbolic
  identities C26a–c and 30 sampled two-cap compositions C26d.
- **Why:** the inverse-first pole lift already has an elementary forward formula on the unblended
  caps, but the transition shell is the remaining local inverse-pair problem. Expanding the
  Hopf matrix of the blended spinor and using
  \(P_LQ_S+Q_SP_L=hI+(\Sigma(n_L)+\Sigma(n_S))/2\)
  gives (11.12c), while \(P_Lw\) is a positive real multiple of the unknown spinor,
  giving (11.12d). This removes the circle phase from the unknowns without asserting
  an elementary inverse for the three-variable orbit map.
- **Verification:** two agents independently checked the projector algebra, the norm and phase
  formulas, the domain restriction and the converse. sym/chart_identities.py: ALL PASS,
  including C26a–c and 30 C26d samples with maximum relative residual 4.1e−67
  (numerical evidence only). All four campaign symbolic suites: ALL PASS.
  verify_construction.py: all checks passed; sphere_products.py:
  Dimension 8 check passed. git diff --check: pass. LaTeX: **not run**
  (latexmk unavailable). The full legacy numerical safety net was not run because no
  construction block changed. Git fetch/pull and push were left to the collaborator as requested.
- **Open points:** no finite elementary inverse of the three-variable shell orbit equation is
  supplied; the stronger global finite elementary inverse-pair goal remains open.

### 2026-10-11-TSG-4 — Exclude scalar smoothing of the south raw pole lift
- **Date:** 2026-10-11
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6)
- **Type:** addition
- **Based on:** 4a7fdc0
- **Where:** docs/s4t2_symbolic_programme.md §11.3 and §11.8;
  extra/s4t2_explicit/sym/chart_identities.py C27. No manuscript or construction file was changed.
- **Change:** derived the south raw-base quadratic jet and the resulting nonpolynomial
  cubic term in the canonical raw spinor lift. Proved that multiplying this lift by
  any real scalar tending to one cannot make it C³ at the pole, even when the scalar
  changes the raw orbit radius. Added exact symbolic jet and divisibility checks C27a–c.
- **Why:** the earlier C25 obstruction excludes smooth lifts over the unchanged raw base.
  A scalar taper changes that base and was therefore not covered. With
  \(C_0(p)=p+\Delta_2(p)+O(|p|^3)\), normalized-projector differentiation gives
  \(g_0(z)=z-\Sigma(\Delta_2(\operatorname{Hopf}z))z/(2|z|^2)+O(|z|^5)\).
  On real spinors \(z=(x,y)\), the imaginary first-component cubic has numerator
  \(\kappa xy^2(x^2-y^2)\), nonzero at the complex root \(x=iy\) of \(x^2+y^2\).
  Thus it cannot be a cubic Taylor polynomial. A real scalar tending to one
  leaves this leading imaginary term unchanged.
- **Verification:** two agents independently derived and reviewed the quadratic jet,
  projector cancellation, uniform remainder, and scalar-taper scope. C27a–c in
  sym/chart_identities.py: ALL PASS. An independent SymPy calculation of the B jet
  and real-plane numerator also gave zero residuals and \(-2i\kappa y^5\) at \(x=iy\).
  All four campaign symbolic suites: ALL PASS. verify_construction.py:
  all checks passed; sphere_products.py: Dimension 8 check passed.
  git diff --check: pass. LaTeX: **not run** (latexmk unavailable).
  The full legacy numerical safety net was not run because no construction block changed.
  Git fetch/pull and push were left to the collaborator as requested.
- **Open points:** this excludes scalar rescaling of the canonical raw lift, not
  a change in Hopf direction, a phase-modified lift, or an elementary inverse
  of the existing blended shell.

### 2026-10-11-TSG-5 — Give the pole replacement an elementary inverse in both directions
- **Date:** 2026-10-11
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6)
- **Type:** addition
- **Based on:** 72f2513
- **Where:** docs/s4t2_symbolic_programme.md §11.3, §11.6–11.8;
  SYMBOLIC_CAMPAIGN.md status; extra/s4t2_explicit/cusp.py,
  check_inverse.py and sym/chart_identities.py. No manuscript file was changed.
- **Change:** replaced the pole-shell blend in the elementary inverse-first
  construction by the radial conjugation \(H_j=D^{-1}g_jD\), joined to the
  identity inside \(|z|^2\leq10^{-8}\). Both \(H_j\) and its inverse now have
  finite elementary formulas. Added the corresponding forward formula on the
  source cap, exact symbolic identities, independent sampled compositions,
  shell tensor checks, and precision guards for evaluating the flat seams.
- **Why:** the earlier smooth blend had no elementary forward inverse. The
  increasing radial map \(S=\psi^{-1}T\phi\) compresses the inner boundary
  to zero and is the identity to all orders at its outer boundary; conjugation
  makes the raw lift's nonsmooth pole jet flat while preserving its explicit
  inverse. The spherical-blowup argument in §11.3 proves smoothness at the
  inner seam. The raw defect bound proves that the image contains the required
  cap and that the local formulas agree with the slab on an open overlap.
- **Verification:** two agents independently reviewed the radial identities,
  seam proof, image/domain bounds and forward-cap signs. C28a–d,h–i give exact
  identities and rational margins; C28e–g give 120 two-cap compositions and
  10 inverse-Jacobian samples. The 18-case production seam check has relative
  round-trip error at most 3.0e−59; its precision shortcuts have explicit
  type-(c) bounds in §11.7. `sym/chart_identities.py`, all three other campaign
  symbolic suites, `check_inverse.py --closed`, and
  `check_inverse.py --closed-nijenhuis`: **ALL PASS**. The legacy checks
  `check_tcore.py`, `check_tlift.py`, `check_lift.py`, `check_degree.py`,
  `check_full.py`, `check_inverse.py`, `check_step1.py`, `check_step2.py`,
  `check_glue.py`, `check_hopf.py`, `nijenhuis_round.py 2000`, and
  `check_jround.py` passed; the fine 1121×128×128 degree grid had positive
  determinant at all 18,366,464 nodes and degree 1.000000.
  `verify_construction.py`: `all checks passed`; `sphere_products.py`:
  `Dimension 8 check passed`; `git diff --check`: pass. LaTeX: **not run**
  (`latexmk` unavailable). Git fetch/pull and push were left to the
  collaborator as requested.
- **Open points:** the inverse of the global height transition in (11.3)
  remains implicit. This commit gives an elementary inverse of the pole
  replacement, not yet a finite elementary global inverse pair.

### 2026-10-11-TSG-6 — Complete the finite elementary inverse pair
- **Date:** 2026-10-11
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6)
- **Type:** statement
- **Based on:** 984e4ba
- **Where:** docs/s4t2_symbolic_programme.md §11.2, §§11.4–11.8;
  SYMBOLIC_CAMPAIGN.md status and T15; extra/s4t2_explicit/cusp.py,
  check_inverse.py, sym/chart_identities.py and sym/i2_margins.py.
  No manuscript file was changed.
- **Change:** replaced the inverse-first scalar height by the terminal-cut
  height (11.3)–(11.5), with finite elementary formulas in both directions,
  exact central-axis and collar charts, and an explicit source-side reversal
  F_cl of the existing G_cl. Recomputed the global projector bound as
  L=10^8000 and the associated fixed finite counts. The campaign statement
  is promoted to a finite elementary inverse-pair and tensor certificate,
  subject to human review. Added exact identity checks, implementation of
  the forward map, two-way sampled checks and exact half-cell seam tests.
- **Why:** the earlier global height had no elementary forward inverse. The
  tuned terminal map sends the nonzero inner endpoint to zero, is flat at its
  outer identity branch, and has an explicit quadratic inverse. A radial
  conjugation preserves the exact node model while keeping an elementary
  inverse; at the outer height interval and λ=1 it matches the previous
  exponential and middle collar exactly. The type-(c) derivative estimates
  in §11.4 give positive lower bounds and the revised finite projector
  count. The source cap selector is proved to cover the patched pole image;
  its complement lies in the free slab. Reversing the finite scalar,
  projector-phase and torus blocks gives both compositions F_cl G_cl=id and
  G_cl F_cl=id. The displayed dual-number rules make J_cl a finite
  elementary coefficient expression. Independent agents audited the scalar
  inverses, constants, source domains, signs and torus/deck transitions.
- **Verification:** sym/chart_identities.py (C29a–h exact identities and C29i
  numerical margin), sym/i2_margins.py (I2v), sym/lemmas_1d.py and
  sym/node_remainder.py: **ALL PASS**; symbolic_identities.py: pass.
  check_inverse.py --closed: **ALL PASS**, including 276 two-way height
  cases (maximum error 1.0e−56) and 38 reduced-count global forward cases,
  including four exact half-cell seams (maximum error 1.2e−54).
  check_inverse.py --closed-nijenhuis: **ALL PASS** for its three sampled
  finite-instance tensor points and nonintegrable controls. The legacy
  check_tcore.py, check_tlift.py, check_lift.py, check_degree.py
  1121 128 128, check_full.py, check_inverse.py 3000, four numpy checks,
  nijenhuis_round.py 2000 and check_jround.py 300 all passed; the fine
  degree grid had degree 1.000000 and positive determinant at 18,366,464
  nodes. verify_construction.py: `all checks passed`;
  sphere_products.py: `Dimension 8 check passed`; git diff --check: pass.
  LaTeX: **not run** (latexmk unavailable). Git fetch/pull and push were
  left to the collaborator as requested.
- **Open points:** this is a mathematical statement change requiring human
  review. The certified counts N and K have 8001 and 16004 decimal digits;
  the exact fixed-count products were not numerically executed, and no
  smaller counts are certified. The old fitted evaluator's I3–I5/T14
  flags remain open. A shorter practical formula is a further task.

### 2026-10-11-TSG-7 — Flag the campaign category-(d) mismatch
- **Date:** 2026-10-11
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6)
- **Type:** flag
- **Based on:** cab13eb
- **Where:** SYMBOLIC_CAMPAIGN.md §0, compared with AGENTS.md §8.2 and
  docs/s4t2_symbolic_programme.md §11.8. Only this ledger entry was added.
- **Change:** flagged a pre-existing mismatch: the campaign introduction
  labels a cited classical theorem as category (d), while binding AGENTS.md
  uses category (d) for a labelled numerical observation. The current §11.8
  follows AGENTS.md; no campaign claim was changed by this flag.
- **Why:** a reader might misclassify a numerical check as a proof step, or
  apply the wrong rule to a cited input. AGENTS.md §8.1 limits edits to
  SYMBOLIC_CAMPAIGN.md to its task list and status, so an agent cannot
  silently rewrite this definition in §0.
- **Verification:** read the exact category lists in both files and checked
  their conflict; git diff --check: pass. Mathematical and build checks:
  **not run** (this entry changes no construction or proof).
- **Open points:** a human collaborator should reconcile the §0 category
  list with AGENTS.md §8.2 or explicitly authorize an edit there.

### 2026-10-12-TSG-1 — Direct second-jet checks on the exact pole and middle charts
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** 2e16096
- **Where:** extra/s4t2_explicit/sym/chart_identities.py (C30),
  extra/s4t2_explicit/sym/node_remainder.py (N8),
  extra/s4t2_explicit/symbolic_identities.py (pure-middle check),
  docs/s4t2_symbolic_programme.md status and §11.8,
  SYMBOLIC_CAMPAIGN.md task T16. No manuscript or construction file changed.
- **Change:** added executable exact second-jet checks for the certified
  projection/phase loop steps and cusp torus restoration; differentiated the
  actual six-coordinate inverse map on both inner pole caps, including the
  finite torus height; and computed the actual middle-chart Jacobian, tensor
  square and all 216 Nijenhuis entries. Added an explicit task and scope record
  for the still-open global source-linked verification.
- **Why:** the previous C10 verified naturality for an arbitrary coordinate
  map, but did not differentiate the specific map in (11.18). On the inner
  caps the pole patch is exactly the identity and the four-seed height reduces
  to (A-Y)/(1-AY); N8 constructs and differentiates this actual map. On the
  pure middle, the displayed torus-restoration formula has determinant
  2π(1+|t|²)((|t|²-1)log|t|+2|t|²)>0, since the logarithmic product is
  nonnegative. In both checks, the actual Hessian is symmetric. The executable
  216-component inverse-map identity then gives N=0; a nonsymmetric-frame
  control produces a nonzero component. C30 verifies the exact source loop
  update identities and the one-step second-jet rules, so the fixed enormous
  counts are handled by induction rather than numerical execution.
- **Verification:** chart_identities.py (C30a–n), node_remainder.py (N8),
  symbolic_identities.py (middle c1–c8), lemmas_1d.py and i2_margins.py:
  **ALL PASS**. N8's six production-value comparisons are numerical only
  (maximum errors 1.1e-72 and 1.1e-75). check_inverse.py --closed and
  --closed-nijenhuis: **ALL PASS** before these edits; their construction
  inputs were not changed. verify_construction.py: all checks passed;
  sphere_products.py: Dimension 8 check passed; git diff --check: pass.
  LaTeX: **not run** (latexmk unavailable). Git fetch/pull and push were
  left to the collaborator as requested.
- **Open points:** these checks do not yet give a direct executable audit of
  the exact free-slab round-to-W map, all scalar/height cells, selector and
  deck transitions, or cap/slab and collar seams. T16 records the remaining
  exact branch-circuit work. The previously stated global theorem still rests
  on the written inverse/gluing certificate and generic naturality, and still
  requires human review.

### 2026-10-12-TSG-2 — Extend source-linked local checks of the specific tensor
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** statement
- **Based on:** 8b74839
- **Where:** extra/s4t2_explicit/sym/chart_identities.py (C31–C32),
  extra/s4t2_explicit/sym/node_remainder.py (N9),
  extra/s4t2_explicit/symbolic_identities.py (d–e),
  docs/s4t2_symbolic_programme.md status and §11.8,
  SYMBOLIC_CAMPAIGN.md task T16. No manuscript or construction file changed.
- **Change:** added exact, source-linked projector and second-jet checks on one
  pure-middle and one two-active-pair cusp cell; identified the actual raw
  projector and cusp reference on both cap/slab overlaps and checked equality
  of their torus-restored values and jets; differentiated the finite torus
  correction on one active south collar chamber; and constructed the actual
  six-real-coordinate map and `J=M^{-1}J_0M` on one south raw-pole annulus.
  The latter uses 43 explicitly defined scalar second-jet nodes and six
  stage-bound Hessian circuits. The documents record these as partial regional
  checks and retain T16 as open.
- **Why:** the source raw inverse, round cap and torus restoration compose to
  a concrete elementary map on the stated annulus. The exact CSE circuit
  reconstructs that map and propagates its first and second derivatives;
  the checked chain-rule skew identity and the existing 216-component
  algebraic identity reduce `N_J` to zero for this particular local `J`.
  The round factor is invertible, (11.11) gives the raw singular-value margin
  >1999/2000, and the restoration determinant is
  `-(1-|t|^4) ≠ 0`. N9 uses source formulas to identify the raw cap
  projector with the slab reference before applying the exact phase and cusp
  actions. C32 checks the nonconstant two-pair projector rather than treating
  a meridian projector as an arbitrary smooth function.
- **Verification:** chart_identities.py C31a–d and C32a–d: **ALL PASS**;
  node_remainder.py N9a–i: **ALL PASS**; symbolic_identities.py d1–d5 and
  e1–e7: **ALL PASS**. The e7 raw-Jacobian comparison is numerical only,
  with one witness: measured \(\|Dg-I\|_{\mathrm{op}}=5.4871\cdot10^{-7}\)
  versus proved bound \(3\cdot10^{-4}\), and sampled minimum singular value
  0.9999999270 versus proved uniform lower bound 0.9995. The d5 production
  comparison is also numerical only (error \(3.9233\cdot10^{-61}\)).
  verify_construction.py: `all checks passed`; sphere_products.py:
  `Dimension 8 check passed`; git diff --check: pass. check_inverse.py
  --closed: **ALL PASS** earlier this turn on the unchanged construction
  code; the numerical construction was not modified. LaTeX: **not run**
  (latexmk unavailable). Git fetch/pull and push were left to the
  collaborator as requested.
- **Open points:** this is not a direct executable verification of the
  specific global `J_cl`. The free-slab round-to-W branch and
  its intervening cusp/collar projectors, remaining height and selector
  cells, the north and radial-transition cap chambers, and remaining seam
  jets are open. Human review is required for the newly checked cap/slab
  overlap status and regional claims. The global conclusion still uses the written
  diffeomorphism/gluing certificate and generic naturality, and requires
  human review.

### 2026-10-12-TSG-3 — Check the nonidentity pole radial stage and its inverse
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** 791c775
- **Where:** extra/s4t2_explicit/sym/node_remainder.py (N10) and
  docs/s4t2_symbolic_programme.md status and §11.8. No manuscript or
  construction file changed.
- **Change:** added source-linked exact scalar nodes for the elementary
  radial compression \(S\) and its inverse \(S^{-1}\), checked the primitive
  inverse identities and derivative factors, and differentiated both radial
  maps in four real coordinates through their 16 first and 64 second entries
  apiece. The exact determinant is \(S(r)S'(r)/r>0\) for
  \(10^{-8}<r<1.5\cdot10^{-8}\). The document records this as a stage check,
  not a check of the complete shell tensor.
- **Why:** on this shell the pole patch is
  \(H_j=D^{-1}\circ g_j\circ D\), so the previously checked raw annulus
  cannot be used without both nonidentity radial stages. The code binds
  the stable expm1/log1p source operations to the exact elementary
  definitions, checks the positive quadratic root, and compares direct
  differentiation of each radial map with its stated Hessian formula.
  The positive scalar derivative follows from the displayed product of
  the positive \(\phi'\), \(T'\), and \((\psi^{-1})'\) factors.
- **Verification:** node_remainder.py N10a–d: **ALL PASS**;
  verify_construction.py: all checks passed; sphere_products.py:
  Dimension 8 check passed; git diff --check: pass.
  LaTeX: **not run** (latexmk unavailable). Git fetch/pull and push were
  left to the collaborator as requested.
- **Open points:** N10 does not compose the radial stages with \(g_j\),
  the round cap, or torus restoration, and does not calculate the specific
  six-dimensional \(J_{\mathrm{cl}}\) on the shell. All global T16 gaps
  recorded in the previous entry remain open; the north raw annulus,
  remaining projector branches, and free-slab path are being audited
  separately. Human review remains necessary for the global conclusion.

### 2026-10-12-TSG-4 — Verify the north raw-pole annulus tensor
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** 0cbc16c
- **Where:** extra/s4t2_explicit/symbolic_identities.py (f) and
  docs/s4t2_symbolic_programme.md status and §11.8. No manuscript or
  construction file changed.
- **Change:** added a source-linked north raw-pole annulus circuit for
  the actual six-real map from round coordinates through the north raw
  inverse and finite torus restoration. It defines \(J=M^{-1}J_0M\)
  from a 28-node scalar second-jet program, checks 64 raw/round and
  216 torus/composite Hessian skews, and applies the previously checked
  36/216-component identity to this named tensor. The programme now
  displays the full local formula and states its annulus scope.
- **Why:** the north Hopf rotation, dominant ratio, \(\nu=e^{-2}\)
  normalization, Cramer orientation and torus frame differ from the
  south annulus. The ordered source guards and exact algebraic nodes
  check those signs. The round factor is nonsingular, (11.11) supplies
  \(\sigma_{\min}(Dg_\infty)>1999/2000\), and the north restoration
  determinant is \((1-|t|^4)/\nu^2>0\) on this chamber.
- **Verification:** symbolic_identities.py (f1–f6): **ALL PASS**.
  The separately labelled numerical witness has
  \(\|Dg-I\|_{\mathrm{op}}=3.9588550449509324\cdot10^{-7}\)
  versus its proved pointwise bound \(3\cdot10^{-4}\),
  sampled \(\sigma_{\min}=1.0000000085917364\) versus the proved
  uniform bound \(0.9995\), and production raw-value error
  \(3.5284\cdot10^{-64}\). verify_construction.py: all checks passed;
  sphere_products.py: Dimension 8 check passed; git diff --check: pass.
  LaTeX: **not run** (latexmk unavailable). Git fetch/pull and push were
  left to the collaborator as requested.
- **Open points:** this certificate covers one open north raw annulus,
  not the nonidentity radial shell, selector seams, or the free-slab
  round-to-W circuit. The complete direct verification of the specific
  global \(J_{\mathrm{cl}}\) remains open, and the global statement
  still requires human review.

### 2026-10-12-TSG-5 — Check the raw-reference profile and selector seam
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** 86064f7
- **Where:** extra/s4t2_explicit/sym/node_remainder.py (N11) and
  docs/s4t2_symbolic_programme.md status and §11.8. No manuscript or
  construction file changed.
- **Change:** added an ordered source-path check for both
  dominant-coordinate positive-root profile branches, exact nonempty
  round-meridian witnesses, finite value/first/second real-ratio jets, and
  the removable \(|u_1|=|u_2|\) selector seam of the effective middle
  coordinate \(t\). The seam calculation proves \(D_iD_d=4e^4\), equality
  of the two positive factors and their first two derivatives, and the
  common source phase input \(u_2\) there.
- **Why:** the raw reference is an input to every free-slab projector step.
  Its two ratio branches use different numerical formulas, so their
  smooth agreement at the selector boundary needs an exact check. The
  script compares the ordered source AST with the written formulas and
  checks every radial derivative identity by SymPy.
- **Verification:** node_remainder.py N11a–e: **ALL PASS**;
  verify_construction.py: all checks passed; sphere_products.py:
  Dimension 8 check passed; git diff --check: pass.
  LaTeX: **not run** (latexmk unavailable). Git fetch/pull and push were
  left to the collaborator as requested.
- **Open points:** this checks one scalar generator and its selector
  seam, not the inverse-first cusp height, ten-slot embedding,
  fixed-count free-slab projection/phase composition, or the global
  six-dimensional \(J_{\mathrm{cl}}\). The complete direct T16
  verification and human review of the global certificate remain open.

### 2026-10-12-TSG-6 — Check a shifted cusp projector and log-cut seam
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** 7e963cb
- **Where:** extra/s4t2_explicit/sym/chart_identities.py (C33–C34),
  docs/s4t2_symbolic_programme.md status and §11.8, and
  SYMBOLIC_CAMPAIGN.md T16 status. No manuscript or construction file
  changed.
- **Change:** added an ordered source-path check for the negative-height
  deck-shifted cusp frame, its unit rank-one four-slot projector, and its
  five-variable first and second jets at both cusp signs. Added a
  source-linked active-collar negative-real-log seam certificate: the
  real Cramer solve and two-coordinate cusp reconstruction, south/north
  deck transformations, surviving angular phase, and matching
  value/first/second chart-coordinate jets for the angular entry and
  cusp input. Strict target witnesses show both seam sectors are
  nonempty.
- **Why:** the free-slab meridian projector crosses the negative-height
  cusp frame and principal-log cuts. These branches cannot be inferred
  from the previously checked positive-height frame. On the collar the
  Cramer determinant has manifest lower bounds \(981/5000\) south and
  \(40\) north, so the coordinate reconstruction is unique.
- **Verification:** chart_identities.py C33a–d and C34a–f: **ALL PASS**.
  The separately labelled numerical observations (d), sample size two,
  give south/north determinants 1.2640336891 and 125.5013159820
  against the proved bounds above. verify_construction.py:
  all checks passed; sphere_products.py: Dimension 8 check passed;
  git diff --check: pass. LaTeX: **not run** (latexmk unavailable).
  Git fetch/pull and push were left to the collaborator as requested.
- **Open points:** C34 differentiates the surviving angular entry and
  reconstructed cusp scalar in six local chart variables; it does not
  differentiate the preceding round-meridian map or all ten embedding
  slots. Other projector cells and seams, the full fixed-count
  projection/phase composition, and the direct global \(J_{\mathrm{cl}}\)
  check remain open.

### 2026-10-12-TSG-7 — Verify one south nonidentity radial-shell tensor
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** 7d265d1
- **Where:** extra/s4t2_explicit/symbolic_identities.py (g),
  docs/s4t2_symbolic_programme.md status and §11.8, and
  SYMBOLIC_CAMPAIGN.md T16 status. No manuscript or construction file
  changed.
- **Change:** added a source-linked six-real second-jet circuit on one
  open south radial-transition chamber. It composes the actual round
  cap, nonidentity \(D\), the existing 43-node south raw inverse,
  nonidentity \(D^{-1}\), and finite torus restoration; stores the
  resulting value, \(6\times6\) Jacobian, six Hessian circuits and
  \(J=M^{-1}J_0M\). Three 64-entry four-dimensional and one
  216-entry six-dimensional Hessian-skew checks pass. The previously
  checked (c8) identity then reduces all 36 square and 216 Nijenhuis
  components for this bound circuit.
- **Why:** N10 checked the radial maps separately, but the actual
  shell tensor contains the nonlinear raw inverse between them.
  Exact bounds on \(1.3\cdot10^{-8}<|z_L|^2<1.4\cdot10^{-8}\) make both
  radial branches and the raw source chamber fixed, avoid runtime
  precision shortcuts, and give the uniform raw singular-value
  margin \(>99993/100000\). The written (11.11), N10 and §11.6
  determinant bounds make the local Jacobian invertible.
- **Verification:** symbolic_identities.py (g1–g6): **ALL PASS** on
  an independent full run; verify_construction.py: all checks passed;
  sphere_products.py: Dimension 8 check passed; git diff --check:
  pass. The separately labelled numerical observation (d), sample
  size one, has production radial-patch value error
  \(4.3267\cdot10^{-84}\), raw
  \(\|Dg-I\|_{\mathrm{op}}=1.1986845\cdot10^{-7}\) against its
  pointwise bound \(6.6998742\cdot10^{-5}\), and sampled
  \(\sigma_{\min}=0.9999999867\). LaTeX: **not run**
  (latexmk unavailable). Git fetch/pull and push were left to the
  collaborator as requested.
- **Open points:** the Nijenhuis calculation uses the specific
  source-bound second jets with the generic polynomial identity (c8);
  it does not expand 216 giant component expressions. This checks
  one south transition chamber, not both full shells, the free-slab
  round-to-\(W\) map, or the global \(J_{\mathrm{cl}}\).

### 2026-10-12-TSG-8 — Check one active inverse-first height stage
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** 6f67208
- **Where:** extra/s4t2_explicit/sym/node_remainder.py (N12),
  docs/s4t2_symbolic_programme.md status, and SYMBOLIC_CAMPAIGN.md T16
  status. No manuscript or construction file changed.
- **Change:** added an ordered source-path check and a 30-node finite
  second-jet circuit for one active inverse-first cusp height stage,
  ending at the terminal-forward output \(v\). It gives the value,
  three first, and six independent second partials in
  \((b,q,\lambda)\), with no unevaluated derivative.
- **Why:** the active terminal correction is a nontrivial scalar
  branch inside every free-slab reference. The exact witness
  \(a=b=3/20,\ |t|=2/5,\ \nu=1,\ \lambda=1/2\) and written rational
  bounds certify its source guards and exclude the working-precision
  shortcuts on a neighborhood.
- **Verification:** node_remainder.py N12a–d: **ALL PASS** on an
  independent full run. The separately labelled one-point numerical
  check (d) gives source value error \(<10^{-75}\), \(q\) lower/upper
  margins \(0.0742842,0.0503158\), and terminal \(\eta/\mathrm{gap}\)
  upper margin \(0.2012902\). verify_construction.py: all checks passed;
  sphere_products.py: Dimension 8 check passed; git diff --check:
  pass. LaTeX: **not run** (latexmk unavailable). Git fetch/pull and
  push were left to the collaborator as requested.
- **Open points:** this circuit stops before the subsequent
  \(G^{-1}/R^{-1}/G\) height chain, \(Y\), \(\rho\), the other height
  cells, and the global six-real \(J_{\mathrm{cl}}\) circuit.

### 2026-10-12-TSG-9 — Verify one north nonidentity radial-shell tensor
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** 3c628fc
- **Where:** extra/s4t2_explicit/symbolic_identities.py (h),
  docs/s4t2_symbolic_programme.md status, and
  SYMBOLIC_CAMPAIGN.md T16 status. No manuscript or construction
  file changed.
- **Change:** added a separate source-linked north six-real
  second-jet circuit on one open radial-transition chamber. It
  composes the actual round cap, nonidentity \(D\), north raw
  inverse, nonidentity \(D^{-1}\), and finite torus restoration;
  stores the value, \(6\times6\) Jacobian, six Hessian circuits,
  and \(J=M^{-1}J_0M\). Three 64-entry four-dimensional and one
  216-entry six-dimensional Hessian-skew checks pass. The previously
  checked (c8) identity reduces all 36 square and 216 Nijenhuis
  components for this bound circuit.
- **Why:** the north profile and cusp signs differ from the checked
  south transition. The strict domain bounds fix both radial calls
  to nonidentity branches, exclude working-precision shortcuts,
  and give uniform raw singular-value margin
  \(>99993/100000\). The written N10, (11.11), and §11.6
  bounds make the local Jacobian invertible.
- **Verification:** symbolic_identities.py (h1–h5): **ALL PASS**
  on an independent full run; verify_construction.py: all checks
  passed; sphere_products.py: Dimension 8 check passed;
  git diff --check: pass. The separately labelled numerical
  observation (d), sample size one, has patch-value error
  \(3.6645\cdot10^{-84}\), raw
  \(\|Dg-I\|_{\mathrm{op}}=8.7626736\cdot10^{-8}\) against its
  pointwise bound \(6.6998742\cdot10^{-5}\), and sampled
  \(\sigma_{\min}=1.0000000019\). LaTeX: **not run**
  (latexmk unavailable). Git fetch/pull and push were left to the
  collaborator as requested.
- **Open points:** as in the south check, the 216 Nijenhuis
  reduction is compositional and uses the source-specific Hessian
  symmetries with the generic identity (c8). It covers one north
  chamber, not the whole shell, free slab, or global \(J_{\mathrm{cl}}\).

### 2026-10-12-TSG-10 — Check an active ten-slot collar projector
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** 00f2b8e
- **Where:** extra/s4t2_explicit/sym/chart_identities.py (C35),
  docs/s4t2_symbolic_programme.md status, and
  SYMBOLIC_CAMPAIGN.md T16 status. No manuscript or construction file
  changed.
- **Change:** added an ordered source-path certificate for one strict
  active upper-negative-axis middle-collar projector cell at both
  cusps. The normalized ten-slot source frame has exactly three
  nonzero entries, with south support \((0,1,9)\) and north support
  \((4,5,9)\). Its projector is rank one. The source flat weight,
  cusp frame, normalization, and projector have exact value, first,
  and second-jet circuits in \((\Re t_c,\Im t_c,\rho,\mathrm{one},H)\);
  the height is on an open constant-weight plateau.
- **Why:** prior checks covered pure-middle and pure-cusp projectors
  but no strict cell in which a cusp pair and an angular middle slot
  coexist. The source frame has \(R_f>0\) and the manifest norm bound
  \(D=(1-s)^2+s^2\ge1/2\). Exact south and north round-meridian
  witnesses certify that the selected source branches are nonempty.
- **Verification:** chart_identities.py C35a–c: **ALL PASS** on an
  independent full run. The separately labelled numerical observations
  (d), sample size two, give \(D=0.9773063821\) south and
  \(0.6619583826\) north, measured margins \(0.4773063821\) and
  \(0.1619583826\) to the proved \(1/2\) bound.
  verify_construction.py: all checks passed; sphere_products.py:
  Dimension 8 check passed; git diff --check: pass. LaTeX:
  **not run** (latexmk unavailable). Git fetch/pull and push were
  left to the collaborator as requested.
- **Open points:** the 3-slot Hessian circuit uses exact generic
  Taylor product rules and source-specific scalar primitives; it
  does not expand all 10×10 Hessian expressions. Other projector
  cells, the fixed-count free-slab composition, and the direct
  global \(J_{\mathrm{cl}}\) check remain open.

### 2026-10-12-TSG-11 — Check selected outer collar seam jets
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** da111ab
- **Where:** extra/s4t2_explicit/sym/node_remainder.py (N13),
  docs/s4t2_symbolic_programme.md status, and
  SYMBOLIC_CAMPAIGN.md T16 status. No manuscript or construction file
  changed.
- **Change:** added source-linked exact checks for one lower-angular
  collar/pure-middle seam cell at each cusp, with
  \(\Arg t=-2\pi/3\), south \(|t|=1/2,\ a=1/20\), and north
  \(|t|=2,\ a=9/20\). The positive-root profiles and source Cramer
  solutions select genuine raw round-meridian inputs. The
  \(\lambda=1\) height and flat endpoint give matching first and
  second jets. The ten-slot projector and final phase-corrected
  middle \(W\) coordinate agree modulo deck, even though the raw
  reference tuples differ by a scalar gauge.
- **Why:** direct well-definedness of the free-slab tensor needs
  source-chart seam jets, and equality of raw reference tuples is a
  false test. The exact Cramer determinants exceed \(5/8\) south
  and \(40\) north, and both height \(q\) values lie in
  \((1/2500,1/8)\).
- **Verification:** node_remainder.py N13a–j: **ALL PASS** on an
  independent full run. The separately labelled numerical
  observations (d), sample size two, give Cramer lower-bound
  margins \(0.6498254818,10.39720771\) and lower/upper \(q\)
  margins \(0.0654992845/0.0591007155\) south and
  \(0.0008069875/0.1237930125\) north. The exact source
  phase modulus and all stated first/second seam identities pass.
  verify_construction.py: all checks passed; sphere_products.py:
  Dimension 8 check passed; git diff --check: pass. LaTeX:
  **not run** (latexmk unavailable). Git fetch/pull and push were
  left to the collaborator as requested.
- **Open points:** this is one strict angular and height cell at
  each outer collar seam. The other seam cells, the complete
  free-slab branch circuit, and direct global \(J_{\mathrm{cl}}\)
  verification remain open.

### 2026-10-12-TSG-12 — Carry a source projector through a six-variable jet fold
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** 3833a12
- **Where:** extra/s4t2_explicit/symbolic_identities.py (i),
  docs/s4t2_symbolic_programme.md status, and
  SYMBOLIC_CAMPAIGN.md T16 status. No manuscript or construction
  file changed.
- **Change:** added a six-variable value, first, and symmetric
  second-jet ring with 28 coefficients per scalar. It binds one
  actual C31 pure-middle ten-slot rank-one projector to the source
  AST and carries an arbitrary prior unit phase through the
  normalized projection recurrence. It checks one right-half-plane
  phase increment, local gauge covariance in all 28 coefficients,
  the exact certified \(N,K\) constants at their source call sites,
  and rejection of a skew coordinate Hessian.
- **Why:** the fixed-count projection and phase products are too
  large to execute step by step. The checked recurrence rules admit
  induction at those exact finite counts once every source projector
  and seam supplies valid jets and the written denominator bounds
  apply. The actual C31 cell is one nontrivial instantiation of
  that recurrence.
- **Verification:** symbolic_identities.py (i1–i6): **ALL PASS**
  on an independent full run; its generic product and unary
  Taylor coefficients, source-cell rank-one projector, all 28
  gauge coefficients, phase increment, and negative control pass.
  verify_construction.py: all checks passed; sphere_products.py:
  Dimension 8 check passed; git diff --check: pass. LaTeX:
  **not run** (latexmk unavailable). Git fetch/pull and push were
  left to the collaborator as requested.
- **Open points:** the fixed-count induction is conditional on
  the other projector cells, seams, and denominator premises.
  This does not execute the full free-slab map or directly verify
  the global \(J_{\mathrm{cl}}\).

### 2026-10-12-TSG-13 — Check one pure-middle negative-real log seam
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** da9e9a6
- **Where:** extra/s4t2_explicit/symbolic_identities.py (j),
  docs/s4t2_symbolic_programme.md status, and
  SYMBOLIC_CAMPAIGN.md T16 status. No manuscript or construction
  file changed.
- **Change:** added a source-linked strict pure-middle
  negative-real-axis log-cut check. The positive-root raw profile
  admits an exact round-meridian witness. The canonical Cramer
  coordinate shifts by one when the local log increases by
  \(2\pi i\), while the ten-slot projector is \(E_{99}\) on both
  sides. For a common nonzero endpoint vector, the phase-corrected
  middle \(W\) representatives differ by the explicit
  \(e_1+e_4\) deck vector, so their local two-jets match.
- **Why:** well-definedness of a directly differentiated tensor
  requires checking argument-cut transitions in source coordinates.
  This pure-middle seam is distinct from the previously checked
  active-collar cut.
- **Verification:** symbolic_identities.py (j1–j5): **ALL PASS**
  on an independent full run. The exact Cramer determinant is
  \(>3/2\) on the strict wedge; at \(t=-3/4\) its excess above
  \(225/128\) is \(175\log(4/3)/256>0\). The separately labelled
  one-point numerical value is \(0.196657666715085\).
  verify_construction.py: all checks passed; sphere_products.py:
  Dimension 8 check passed; git diff --check: pass. LaTeX:
  **not run** (latexmk unavailable). Git fetch/pull and push were
  left to the collaborator as requested.
- **Open points:** the endpoint comparison presumes the finite
  projection denominators and a common endpoint vector. Other
  source paths and log seams, the full free-slab circuit, and
  direct global \(J_{\mathrm{cl}}\) verification remain open.

### 2026-10-12-TSG-14 — Complete one active height log-radius circuit
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** 84ec2ae
- **Where:** extra/s4t2_explicit/sym/node_remainder.py (N14),
  docs/s4t2_symbolic_programme.md status, and
  SYMBOLIC_CAMPAIGN.md T16 status. No manuscript or construction
  file changed.
- **Change:** extended the strict N12 active south height witness
  through the production \(G^{-1}\), \(R^{-1}\), \(G\), positive
  quadratic \(Y\), and log-radius \(\rho\) branches. The full
  40-node circuit has a value, three first, and six independent
  second partials in \((b,q,\lambda)\), without unevaluated
  derivatives. Exact outward rational intervals certify the
  branch guards and positive denominators.
- **Why:** the previous source-bound circuit stopped at the
  terminal-forward value. A direct free-slab tensor needs the
  height derivative through \(\rho\), including the inverse
  \(R\) branch.
- **Verification:** node_remainder.py N14a–d: **ALL PASS** on
  an independent full run. The type-(c) rational interval
  certificate uses a \(10^{-30}\) grid, exponential Taylor
  degree 50 with an explicit geometric tail, logarithm odd
  degree 101 with an explicit tail, and exact integer square-root
  floors. It proves \(G^{-1}(v)\in(1/10,19/100)\),
  \(R^{-1}G^{-1}(v)\in(1/10,1/5)\),
  \(Y\in(1/10,1/4)\), and \(\rho\in(-1,-1/2)\).
  The separately labelled one-point numerical comparison (d)
  has \(|\rho_{\rm DAG}-\rho_{\rm source}|=3.1109\cdot10^{-91}\);
  measured \(G^{-1}\) lower/upper margins
  \(0.06316123157/0.02683876843\), \(R^{-1}<1/5\)
  margin \(0.03670151219\), and \(Y\) lower/upper margins
  \(0.09136272386/0.05863727614\).
  verify_construction.py: all checks passed; sphere_products.py:
  Dimension 8 check passed; git diff --check: pass. LaTeX:
  **not run** (latexmk unavailable). Git fetch/pull and push were
  left to the collaborator as requested.
- **Open points:** this is one strict active south height cell.
  Other heights, full six-real free-slab composition, and global
  direct \(J_{\mathrm{cl}}\) verification remain open.

### 2026-10-12-TSG-15 — Factor two local free-slab Jacobians
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** a9d8092
- **Where:** extra/s4t2_explicit/sym/chart_identities.py (C36),
  docs/s4t2_symbolic_programme.md status, and
  SYMBOLIC_CAMPAIGN.md T16 status. No manuscript or construction
  file changed.
- **Change:** factored the six-real Jacobian of the source-linked
  round-to-middle-\(W\) map on strict neighborhoods of both C35
  witnesses. The proof follows the raw quotient, rational profile,
  active height, Cramer, finite-projection fibre phase, and torus
  restoration branches. In the fixed target order
  \((\Re t,\Im t,\Re z,\Im z,\Re w,\Im w)\), both determinants have
  the explicit negative factorizations recorded in the programme;
  every remaining factor is manifestly positive.
- **Why:** direct verification of \(J=M^{-1}J_0M\) requires a
  nonzero actual source Jacobian, including the free slab. This
  gives the first source-linked six-real determinant cells there.
- **Verification:** chart_identities.py C36a–c: **ALL PASS**
  on an independent full run. SymPy checks the raw signed Jacobian,
  both ratio/Cramer signs, the positive rational-profile
  derivatives, height derivative, source-path AST, fibre
  equivariance, torus-restoration determinant, and the two
  factorizations. The exact \(1-s^2\) factor exceeds \(15/16\)
  on \(9/100<s<1/4\); C35 provides the strict source witnesses.
  verify_construction.py: all checks passed; sphere_products.py:
  Dimension 8 check passed; git diff --check: pass. LaTeX:
  **not run** (latexmk unavailable). Git fetch/pull and push were
  left to the collaborator as requested.
- **Open points:** this uses the prior finite-projection
  nonvanishing bound and covers only two free-slab neighborhoods.
  Their full second-jet \(J\) circuits, other cells and seams,
  and global direct \(J_{\mathrm{cl}}\) verification remain open.

### 2026-10-12-TSG-16 — Factor the exact finite projection phase
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** a940814
- **Where:** extra/s4t2_explicit/symbolic_identities.py (k),
  docs/s4t2_symbolic_programme.md status, and
  SYMBOLIC_CAMPAIGN.md T16 status. No manuscript or construction
  file changed.
- **Change:** derived the production \(N+1\) normalized
  rank-one projection result as the final unit reference times
  an initial normalized overlap and \(N\) consecutive normalized
  scalar overlaps. The executable check binds the source loop
  and exact certified \(N\), verifies the scalar update and all
  28 coefficients of its six-variable second jet, telescopes
  a variable reference gauge, and checks the local argument/deck
  convention. A same-cell C31 reference gives one positive
  source-specific adjacent overlap.
- **Why:** the certified projection count cannot be expanded
  literally. The scalar formula gives a finite exact expression
  with an inductive jet rule that can be used by a direct
  six-real \(h\) circuit when its path cells and denominators
  are certified.
- **Verification:** symbolic_identities.py (k1–k4): **ALL PASS**
  on an independent full run. The algebraic one-step identity
  and six-variable Taylor coefficients are exact; the C31
  same-cell overlap is manifestly positive. The full-count
  induction is conditional on nonzero overlaps and source
  reference jets at every sample.
  verify_construction.py: all checks passed; sphere_products.py:
  Dimension 8 check passed; git diff --check: pass. LaTeX:
  **not run** (latexmk unavailable). Git fetch/pull and push were
  left to the collaborator as requested.
- **Open points:** this result alone does not instantiate all
  meridian reference cells, certify every argument seam, or
  construct the full free-slab \(h/M/H/J\) circuit. Global
  direct \(J_{\mathrm{cl}}\) verification remains open.

### 2026-10-12-TSG-17 — Certify one all-mesh free-slab cusp cell
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** f035a49
- **Where:** extra/s4t2_explicit/sym/chart_identities.py (C37),
  docs/s4t2_symbolic_programme.md status, and
  SYMBOLIC_CAMPAIGN.md T16 status. No manuscript or construction
  file changed.
- **Change:** certified a strict near-south round chamber in which
  every meridian reference of the source-defined fixed-count
  projection takes the same positive-\(x\) south-cusp branch.
  The exact source path, open witness, uniform \(>3/5\) initial
  and adjacent overlap bounds, four-seed reduction
  \(H_f=a/(1-q)\), and positive factored six-real Jacobian are
  recorded with their coordinate orders. The checker uses
  induction over the certified count; the production evaluator
  requires max_work=None to name this astronomical finite
  expression.
- **Why:** a direct source-specific tensor circuit needs both
  valid differentiable branch formulas along the entire fixed
  projection path and an invertible six-real Jacobian.
- **Verification:** chart_identities.py C37a–d: **ALL PASS**
  on an independent full run; a second independent reviewer
  also ran the full script to **ALL PASS** and checked the
  determinant signs and every uniform branch inequality.
  The exact projection denominator bound is \(>3/5\);
  the type-(d) one-point initial-overlap sample is
  \(0.9999999999999999999512\), measured excess
  \(0.3999999999999999999512\). The separate \(N=1\)
  endpoint sample is not a certified-mesh minimum.
  verify_construction.py: all checks passed; sphere_products.py:
  Dimension 8 check passed. LaTeX: **not run**
  (latexmk unavailable). Git fetch/pull and push were left
  to the collaborator as requested.
- **Open points:** C37 supplies no six Hessians or named \(J\).
  Other free-slab cells, seams, and the global direct
  \(J_{\mathrm{cl}}\) verification remain open.

### 2026-10-12-TSG-18 — Verify one source-bound all-mesh free-slab tensor
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** d0f1822
- **Where:** extra/s4t2_explicit/symbolic_identities.py (l),
  docs/s4t2_symbolic_programme.md status and certificate, and
  SYMBOLIC_CAMPAIGN.md T16 status. No manuscript or construction
  file changed.
- **Change:** composed the actual source meridian reference and
  fixed \(N=4\cdot10^{8000}+1\) normalized projection product on
  one strict south free-slab chamber. Forty source scalar nodes
  have six first and 21 symmetric second-jet assignments each.
  The initial and later overlap schemas and their exact product
  induction give the endpoint gauge; the source four-seed height,
  Cramer correction and torus restoration then give a six-real
  \(h\), its \(6\times6\) Jacobian \(M\), six Hessians, and the
  named unevaluated concrete \(J=M^{-1}J_0M\). The 36 square
  and 216 Nijenhuis reductions use these actual Hessians and
  the checked c8 algebraic identity. The previously checked
  C37d determinant proves \(M\) invertible on this subchamber.
- **Why:** this is a direct executable symbolic check of the
  specific tensor on an open free-slab cell, including the
  source's astronomical fixed projection count.
- **Verification:** symbolic_identities.py (l1–l7): **ALL PASS**
  on a full run. The independently checked C37a–d source
  chamber and determinant give the nonzero Jacobian. The
  separate type-(d) \(N=1\) six-output source comparison has
  error \(5.6882\cdot10^{-81}\); it is a transcription check
  of a different finite map, not a certified-count execution.
  verify_construction.py: all checks passed; sphere_products.py:
  Dimension 8 check passed; git diff --check: pass. LaTeX:
  **not run** (latexmk unavailable). Git fetch/pull and push
  were left to the collaborator as requested.
- **Open points:** only one strict south free-slab chamber is
  certified here. Other free-slab cells, chart seams, and the
  global direct \(J_{\mathrm{cl}}\) verification remain open.

### 2026-10-12-TSG-19 — Check the pure-cusp argument-selector overlap
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** a18964d
- **Where:** extra/s4t2_explicit/sym/chart_identities.py (C38),
  docs/s4t2_symbolic_programme.md status and certificate, and
  SYMBOLIC_CAMPAIGN.md T16 status. No manuscript or construction
  file changed.
- **Change:** proved an open source argument-selector overlap
  in both pure cusp charts. Exact round witnesses and strict
  perturbations meet both selector sides. On the lower wedge
  the source heights differ by one deck unit, but the resulting
  ten-slot reference and projector agree with all first and
  second jets. The four-seed height shifts by one; after torus
  restoration, the deck transition is
  \((z_1,z_2,z_3)\mapsto(1/(\nu z_2),\nu z_1z_2^2,
  z_3+2\sigma z_1z_2)\), with \(\sigma=1\) south and
  \(-1\) north and complex determinant \(1\). The upper
  wedge has identical source argument values. The projector
  equality propagates through any fixed-count meridian step
  with nonzero projection denominator.
- **Why:** equality on an open overlap checks the normal as
  well as tangential jets of this artificial selector, which
  a seam-value comparison alone would not establish.
- **Verification:** chart_identities.py C38a–e: **ALL PASS**
  on the author's full run, an independent full run, and a
  separate referee run. The source AST, exact round inverse,
  deck relation, height shift, frame/projector two-jets and
  holomorphic transition were independently checked. The
  proved \(|t_c|<1/10\) guard has separately labelled
  type-(d) lower-wedge samples \(0.0039838456\) south and
  \(0.0294368583\) north, with measured gaps \(0.0960161544\)
  and \(0.0705631417\). verify_construction.py: all checks
  passed; sphere_products.py: Dimension 8 check passed;
  git diff --check: pass. LaTeX: **not run** (latexmk
  unavailable). Git fetch/pull and push were left to the
  collaborator as requested.
- **Open points:** fold propagation uses the established
  nonzero-overlap premise for other meridian cells. Other
  selectors and the global direct \(J_{\mathrm{cl}}\)
  verification remain open.

### 2026-10-12-TSG-20 — Expose exact-formula branches for symbolic auditing
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** a38fbb9
- **Where:** extra/s4t2_explicit/cusp.py,
  extra/s4t2_explicit/symbolic_identities.py (m1–m3),
  docs/s4t2_symbolic_programme.md status and §11.7–11.8, and
  SYMBOLIC_CAMPAIGN.md T16 status. No manuscript file changed.
- **Change:** added an opt-in context that follows the six
  elementary source formula branches formerly replaced by
  working-precision identity returns. The default numerical
  evaluator is unchanged. Source AST checks identify all six
  guards and test nested and exception-safe restoration. Two
  type-(d) strict-interior samples exhibit formula/default value
  differences: \(9.0129315017\cdot10^{-99}\) for the pole radial
  map and \(3.23193414504\cdot10^{-92}\) for the height radial
  map at 240 digits; the symbolic programme separately proves
  the strict radial inequalities on their open active regions.
- **Why:** source comparisons in direct symbolic checks must
  follow the finite elementary formulas, including where a
  finite-precision shortcut would return the identity.
- **Verification:** symbolic_identities.py (m1–m3): **ALL PASS**
  on the full run; an independent referee checked that unwrapping
  the six opt-in gates recovers the original source AST, exercised
  nested/exception restoration and endpoint cases, and reproduced
  both samples. chart_identities.py: **ALL PASS** before the
  separate C39 edit; node_remainder.py: **ALL PASS**.
  check_tcore.py, check_tlift.py, check_lift.py, check_full.py,
  check_inverse.py 3000, check_inverse.py --closed,
  check_inverse.py --closed-nijenhuis, check_degree.py 1121 128 128,
  nijenhuis_round.py 2000, check_step1.py, check_step2.py,
  check_glue.py, check_hopf.py, and check_middle.py: pass.
  check_jround.py 300: **ALL PASS** (including 300 uniform
  points and 100 points at each chart switch).
  verify_construction.py:
  all checks passed; sphere_products.py: Dimension 8 check passed;
  git diff --check: pass. LaTeX: **not run** (latexmk
  unavailable). Git fetch/pull and push were left to the
  collaborator as requested.
- **Open points:** this mode uses finite-precision mpmath
  arithmetic and does not execute the astronomical fixed counts.
  The global direct \(J_{\mathrm{cl}}\) second-jet and overlap
  certificate remains open.

### 2026-10-12-TSG-21 — Certify a dominant-coordinate selector overlap
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** 2b503c6
- **Where:** extra/s4t2_explicit/sym/chart_identities.py (C39a–e),
  docs/s4t2_symbolic_programme.md status, §11.8 and certificate,
  and SYMBOLIC_CAMPAIGN.md T16 status. No manuscript or
  construction file changed.
- **Change:** proved that both source dominant-ratio positive
  roots give the same effective \(t\) and real two-jets on a
  strict open band across \(|u_1|=|u_2|\). Explicit constants
  keep this band in the active south collar, with nonzero
  Cramer and torus determinants and reference norm. Both
  branches therefore have the same actual three-slot
  normalized reference, projector and restored output through
  second order. An exact round seam witness and perturbations
  on both selector sides are checked. A separate one-point
  type-(d) production sample records the measured margins.
- **Why:** the local \(J\) circuits must agree across an
  artificial source selector before a global direct certificate
  can be claimed.
- **Verification:** chart_identities.py C39a–e: **ALL PASS**
  on my full run and independently on the author's and
  referee's full runs. The referee also checked the source
  reference and embedding at the algebraic seam witness and
  strict side witnesses. The proved margins are
  \(q>1/100,\ q<1/8\), Cramer \(>981/5000\), torus
  \(>15/16\), and norm square \(\ge1/2\). The one-point
  observed gaps are respectively \(0.0623622466\),
  \(0.0526377534\), \(1.0927862054\),
  \(0.0441843611\), and \(0.1888800593\).
  verify_construction.py: all checks passed;
  sphere_products.py: Dimension 8 check passed;
  git diff --check: pass. LaTeX: **not run** (latexmk
  unavailable). Git fetch/pull and push were left to the
  collaborator as requested.
- **Open points:** equality propagates through each finite
  projector step under its nonzero-denominator premise; this
  certificate does not execute the global fixed-count loop.
  Other source seams and global \(J_{\mathrm{cl}}\) coverage
  remain open.

### 2026-10-12-TSG-22 — Check a north flat infinity selector join
- **Date:** 2026-10-12
- **Collaborator:** Tancredi Schettini Gherardini
- **Agent:** Codex (GPT-6), with subagents
- **Type:** addition
- **Based on:** e38ba92
- **Where:** extra/s4t2_explicit/sym/chart_identities.py (C40a–f),
  docs/s4t2_symbolic_programme.md status, §11.8 and certificate,
  and SYMBOLIC_CAMPAIGN.md T16 status. No manuscript or
  construction file changed.
- **Change:** checked one north \(|\mathrm{small}|=1/2\)
  source selector chamber with exact round witnesses at
  \(s=499/1000,1/2,501/1000\). The infinity-side
  height and sparse ten-slot projector have matching
  second jets with the middle side at this flat seam.
  The terminal-height identity at \(\lambda=1\),
  source Cramer signs, flat radial and embedding weights,
  mixed projector jets, local phase cancellation and
  strict four-seed torus domain are checked. The
  final phase/restoration join is conditional on matching
  finite-fold endpoint jets, a nonzero final phase,
  and nonzero projection denominators. A two-point
  type-(d) production sample reports margins.
- **Why:** a global direct tensor certificate needs exact
  second-jet joins at source selectors, including the
  normal derivative at flat boundaries.
- **Verification:** chart_identities.py C40a–f, including
  c2 and c3: **ALL PASS** on the author's and my
  independent full runs.
  An independent referee re-derived the Cramer, height,
  flat cutoff, projector, phase and torus algebra and
  audited the conditional scope. Proven bounds are
  \(q>1/2025\), Cramer \(>30\), torus \(>9/10\),
  and norm square \(>1/4\). The two-point type-(d)
  minimum observed excesses are respectively
  \(0.000709386291679574016\),
  \(19.9925817884944262\), \(0.036998497999\)
  and \(0.783279071736023484\); the left strict
  interior height difference is
  \(-1.37624625388097242\cdot10^{-87}\) at
  190 digits. verify_construction.py: all checks passed;
  sphere_products.py: Dimension 8 check passed;
  git diff --check: pass. LaTeX: **not run**
  (latexmk unavailable). Git fetch/pull and push were
  left to the collaborator as requested.
- **Open points:** intermediate meridian projector joins
  and their simultaneous intersections are not covered by
  C40. The fixed-count global \(J_{\mathrm{cl}}\) direct
  second-jet certificate remains open. The next unchecked
  targets recorded in the T16 task list are the cusp/middle
  \(|t_c|=3/10\) join at both signs, the pure-middle
  \(|t|=1\) join and its argument-cut corner, and the
  all-branch guarded six-real circuit.

### 2026-10-12-EH-1 — Certify the cusp/middle |t_c|=3/10 source join at both signs
- **Date:** 2026-10-12
- **Collaborator:** Edward Hirst
- **Agent:** Claude Code (Claude Opus 5.5), with one independent reviewer subagent
- **Type:** addition
- **Based on:** ea14d2a
- **Where:** extra/s4t2_explicit/sym/chart_identities.py (new section C41a–f
  at the end), docs/s4t2_symbolic_programme.md (status line, §11.8 table
  row and certificate paragraph), SYMBOLIC_CAMPAIGN.md (T16 status and
  next tasks). No manuscript or construction file changed.
- **Change:** added C41, the first T16 target listed by TSG-22: the source
  selector \(|t_c|=3/10\) between the c0/ci core reference
  (`elementary_cusp_reference`) and the middle data built from
  `elementary_cusp_logradius` and the Cramer solve, at both cusp signs,
  including the restored holomorphic chart transition. Let \(M_0\) be the
  right-hand formulas with \(\lambda\), \(1-b_0\) (or \(1-b_\infty\)) and
  \(b_m\) set to zero. C41 proves the following along the whole seam,
  for each argument chart of \(a\) and each fixed branch of \(\log t\):
  (b) the core reference is \(g\) times the \(M_0\) reference
  \((t_c/e^\rho,e^\rho)\), with \(g=t_c/|t_c|\) for \(a>0\) and \(g=1\) for
  \(a\le0\), in all three cells, so the ten-slot projectors are
  identical; (c) the right-hand source differs from \(M_0\) only by terms
  with a factor flat at \(|t_c|=3/10\); (d) for every \(\phi\) and a common
  final unit vector, the c0/ci torus restoration equals the paper's
  holomorphic chart change \(\Phi\) (south) or \(\Theta\) (north) of the
  middle restoration, via \(D-cc=H_fc\). The complex determinants are
  \(i\) and \(-i/t^4\), and \(\Phi\), \(\Theta\) are bound to
  `tfull.c0_to_mid` and `ci_to_mid`. C41a binds the source, C41e gives
  the band and the exact witnesses, and C41f is a type-(d) production
  comparison.
- **Why:** the direct certificate of the specific \(J_{\mathrm{cl}}\)
  needs every source selector to join with matching jets. This one
  changes chart (c0/ci versus middle), so equality of raw tuples is not
  the right test: the references differ by a unit gauge, and the outputs
  differ by a holomorphic chart change. Because the core reference
  equals the \(\lambda=0\) middle continuation identically, the join
  holds along the whole seam, not in one chamber as in C39–C40.
- **Verification:** `python3 extra/s4t2_explicit/sym/chart_identities.py`:
  **ALL PASS** (C1–C41; sympy 1.14.0, mpmath 1.3.0, conda env
  `nijenhuis`, Python 3.10.19). Negative controls, run separately in a
  scratch script, fail as they should: a mutated source line, the north
  correction \(D\) used at the south, a non-flat factor, and a wrong cell
  formula. Explicit constants (type (c)): on \(2/25<|t_c|^2<1/10\),
  \(1/25<q<1/10\) south and \(1/1375<q<1/10\) north (both in
  \((1/2500,1/8)\)), Cramer determinant \(>2r^2(1+r^2)>0\), torus
  determinants \(1-|t_c|^4>0\) and \(|t_c|^4-1<0\), and \(|a|<3/5\) (from
  \(\cos(\pi/5)<9/10\)); the witnesses have \(|X|<3/5\). C41f, numerical
  observation (d), sample size 18 at 120 digits. Per sign it uses four
  pairs (a, arg t_c) covering a<=0 (with the deck frame), a=0, a>0 in
  cell 0, and a>0 in cell 1 (half-centred chart), at the left and seam
  points, plus one right point. Maximum projector difference
  \(2.4\cdot10^{-121}\); maximum restored-output difference
  \(3.9\cdot10^{-120}\); right-side flat weights \(3.8\cdot10^{-87}\);
  \(|\rho-\rho(\lambda=0)|\le1.1\cdot10^{-87}\). At the exact seam,
  rounding sends the production raw reference to either branch, so the
  core formula is evaluated directly there. An independent reviewer
  subagent audited C41 adversarially. It re-derived (b) and (d) by hand
  and ran its own checks: 1200 production points over all cells, ten
  values of arg t (including the cut) and both signs, with worst
  identity errors \(6\cdot10^{-79}\); 360 right-side points showing the
  flat order; a separate symbolic proof of the restore link; and a
  reduced-count whole-chain comparison. Verdict: sound, with caveats.
  All its actionable points are fixed in this commit:
  - the lattice solve (`ez`, `ew`, the determinant helper, `den`) is now
    bound, and the helper is checked to be the cofactor determinant;
  - \(p_m=g\,p_c\) is derived from the conjugate-linear source inner
    product, with a linear-slot negative control, instead of being
    assumed;
  - the restore loop's `yy` is linked symbolically to the final c0/ci
    \(W\)-point at both signs, instead of assuming a common \(H_f\);
  - three trivial checks are replaced by real ones (composition with
    the `tfull` chart maps, \(1+\Re u_1=1-k/Q>0\), and the identity
    \(\det-2r^2(1+r^2)=(1+r^2)(r^2-1)\log r\ge0\) for every \(r>0\));
  - the exposition now states that (d) is exact for \(E\) in the image
    of \(P_{M_0}\) and holds up to flat terms on the right, and that
    only second order is checked;
  - C41f is widened to all cells and both argument charts.
  The full script was rerun after these fixes: ALL PASS (127 s).
  Baseline before the edit: `lemmas_1d.py`, `chart_identities.py`,
  `node_remainder.py`, `i2_margins.py` and `symbolic_identities.py`:
  **ALL PASS**. `verify_construction.py`: `all checks passed`;
  `sphere_products.py`: `Dimension 8 check passed` (before and after).
  `git diff --check`: pass. LaTeX: **not run** (no manuscript change).
  The legacy numerical safety net was not run, because no construction
  block (`cusp.py`, `tcore.py`, `tlift.py`, …) changed.
- **Open points:** the reviewer found a numerical edge (type (d)) in the
  production evaluator, not in the mathematics. The selector compares
  with `mp.mpf('.3')` at dps+50, while the embedding's guard
  \(|xy|\le .3\) uses dps+40, where `.3` can round lower. A point within
  about \(10^{-(\mathrm{dps}+40)}\) of the seam could therefore select
  c0 and then raise in the embedding. It was triggered only in a direct
  call with mixed precision, never through the production chain (80
  exact-seam trials). `cusp.py` is left unchanged; fixing this needs a
  construction edit together with the safety net. Propagation through
  the fold is still conditional on
  the nonzero projection denominators and on the other source seams. C38
  certifies the argument selector only for \(|t_c|<1/10\); its
  continuation to \(1/10\le|t_c|<1/2\), which meets this seam at a corner,
  is not covered, and is now listed as a T16 task. The pure-middle
  \(|t|=1\) selector and its argument cuts remain open. The global
  \(J_{\mathrm{cl}}\) statement still requires human review.

### 2026-10-12-EH-2 — Certify all pure-middle source seams by one deck-invariance lemma
- **Date:** 2026-10-12
- **Collaborator:** Edward Hirst
- **Agent:** Claude Code (Claude Opus 5.5), with one independent reviewer subagent
- **Type:** addition
- **Based on:** 36c66ee
- **Where:** extra/s4t2_explicit/sym/chart_identities.py (new section C42a–f
  at the end), docs/s4t2_symbolic_programme.md (status line, §11.8 table
  row and certificate paragraph), SYMBOLIC_CAMPAIGN.md (T16 status and
  next tasks). No manuscript or construction file changed.
- **Change:** added C42, the next T16 target listed by TSG-22: the
  pure-middle \(|t|=1\) selector and its argument cuts. C42 proves one
  lemma that holds on every open subset of \(1/2<|t|<2\). Let \(G\)
  be generated by \(x_1\mapsto x_1+1\), \(x_4\mapsto x_4+1\) and the
  log-branch change \((x_1,L)\mapsto(x_1-x_4,L+2\pi i)\). Modulo \(G\),
  the source middle data equal data built from continuous branches of
  \(\arg u_1\), \(\arg u_2\), \(\arg t\) and \(\log t\). Hence the
  ten-slot projector is that of the continuous data. For a final fold
  vector \(E\) that is smooth there, \(\Psi(W)\) changes by the lattice
  vectors \((m-kn)e_1+ne_4(L)\), so the restored output changes by
  holomorphic deck translations (complex determinant \(1\)). The
  smoothness of \(E\) is a premise; it rests on the projectors along
  the whole meridian, hence on the other seams and on nonzero overlaps. The lemma covers the \(|t|=1\) selector, where the two
  continuous \(x_4\)-formulas differ by an integer, and the principal
  cuts of \(u_1\), \(u_2\) and \(t\), with their corners. Section (j) of
  `symbolic_identities.py` is the special case of the \(t\)-cut on
  \(1/2<|t|^2<3/4\).
- **Why:** the direct certificate of the specific \(J_{\mathrm{cl}}\)
  needs every artificial selector to glue. In the pure middle all the
  selectors are integer ambiguities of lattice coordinates, so a
  group-invariance argument treats them all at once, instead of one
  chamber at a time. Integer ambiguities in \(x_1\) or \(x_4\) are
  harmless only because the embedding is invariant under them and the
  restoration commutes with lattice translations; C42b and C42d check
  both facts.
- **Verification:** `python3 extra/s4t2_explicit/sym/chart_identities.py`:
  **ALL PASS** (C1–C42). C42a proves that the pure middle is the
  inverse-ratio branch: in the direct branch \(|\mathrm{small}|\le e^{-1}\)
  (type (b): \(\mathrm{root}\ge2\sqrt{B\,rr}\)). C42b–d are identities
  (type (a)). C42e gives exact witnesses at \(|t|=999/1000,\,1,\,1001/1000\)
  (\(\theta=\pi/2\)) and at the \(t\)-cut corner (\(\theta=\pi\)), with
  \(|X|<3/5\). C42f, numerical observation (d), sample size 32: across eight
  seams (\(|t|=1\); the \(t\)-cut at \(|t|=9/10,\,11/10,\,7/5\); the corner;
  the \(u_1\)-cut at \(11/10\) and \(7/5\); the \(u_2\)-cut), the jump of
  \(\Psi\) between production points at distance \(\varepsilon\) is an
  integer lattice vector plus a residual of order \(\varepsilon\): at most
  \(2.9\cdot10^{-12}\) at \(\varepsilon=10^{-12}\), and \(10^{-12}\) times
  smaller at \(\varepsilon=10^{-24}\). The projector jump is at most
  \(3.0\cdot10^{-123}\). An independent reviewer subagent audited C42.
  It re-derived every step against `cusp.py`, found no missed selector,
  and ran its own production tests: 300 random pure-middle points; 300
  random branch choices; 50 seam crossings at \(\varepsilon=10^{-15}\) and
  \(10^{-30}\), including every pairwise and the triple corner of the
  cuts with \(|t|=1\) and points in the \(\lambda\) band; and the full
  production \(G\) with a reduced fold count across ten seams. Verdict:
  sound, with caveats. Fixed in this commit:
  - the premise that \(E\) is smooth is now stated in every place where
    the restored-output conclusion appears;
  - the scope is stated as open subsets of \(1/2<|t|<2\);
  - vacuous checks are replaced: principal-\(\Arg\) values of
    \(\beta_9\) at exact points, the weight supports, and the
    flat-step arguments as factored identities;
  - C42f now includes the \(\lambda\) band.
  The full script was rerun after these fixes: ALL PASS.
  `verify_construction.py`: `all checks passed`;
  `sphere_products.py`: `Dimension 8 check passed`. `git diff --check`:
  pass. LaTeX: **not run** (no manuscript change). The legacy numerical
  safety net was not run (no construction block changed).
- **Open points:** the reviewer notes that `elementary_tensor_numeric`
  (`cusp.py`) differentiates \(G\) directly in the chart. Within
  `mp.diff`'s step of a pure-middle seam it would differentiate across a
  lattice jump, so \(J\) must be read through the quotient (that is,
  through the continuous data). This is not a C42 defect. The joins of
  the pure middle with the collars at
  \(|t|=1/2\) and \(|t|=2\) are covered only in the chambers of N13 and C40.
  The C38 argument selector for \(1/10\le|t_c|<1/2\) and the
  all-branch six-real circuit with global determinant and overlap bounds
  remain open. The global \(J_{\mathrm{cl}}\) statement still requires
  human review.

### 2026-10-12-EH-3 — Extend the argument-selector overlap to every |t_c|<1/2
- **Date:** 2026-10-12
- **Collaborator:** Edward Hirst
- **Agent:** Claude Code (Claude Opus 5.5), with one independent reviewer subagent
- **Type:** addition
- **Based on:** d88746a
- **Where:** extra/s4t2_explicit/sym/chart_identities.py (new section C43a–f
  at the end), docs/s4t2_symbolic_programme.md (status line, §11.8 table
  row and certificate paragraph), SYMBOLIC_CAMPAIGN.md (T16 status and
  next tasks). No manuscript or construction file changed.
- **Change:** added C43, which extends C38 (certified only for
  \(|t_c|<1/10\)) to every \(|t_c|<1/2\), at both cusp signs. On the lower
  overlap \(-\pi<\Arg p<-5\pi/6\), where the half-centred chart has
  \(a_1=a_0+1\), the two charts give the same ten-slot projector. For a
  smooth common final fold vector (a stated premise), the restored
  outputs are related by the holomorphic cusp deck \(T\) of C38d in the
  core, and by the lattice vector \(+e_4\) (south) or \(-e_4\) (north) in
  the collar. The key identity is an exact deck shift
  \(H_f(\hat\psi p)=H_f(p)+1\) of the four-seed height on \(1/4<Y_0<2\),
  for every \(q<1/8\), obtained from \(K(u)+K(1/u)=1\) and the linear
  seeds; it replaces C38c's seed-branch computation, which needed
  \(q<1/100\). The other ingredients are:
  - \(Y_0\in(1/4,2)\) on both charts, in the core and the collar;
  - \((x_1,y_1)=\hat\psi(x_0,y_0)\) in the core, for every moment \(f\);
  - in the collar, \(\rho_1=\rho_0+P\) and the lattice identity
    \((-iP,0)=\pm(e_4-c\cdot e)\mp(\theta/2\pi)e_1\): the ten-vector
    changes by \(e^{-i\theta_c}\), the corrected \(x_1\) is unchanged
    modulo \(1\), and the correction \(D\) is invariant.
- **Why:** this selector meets the C41 seam at a corner and lies on
  every meridian that crosses the lower wedge of the core or the collar.
  C38's height computation depended on a small-\(q\) seed regime; the
  telescoping identity of §10.2 makes the shift exact for every \(q\)
  that occurs.
- **Verification:** `python3 extra/s4t2_explicit/sym/chart_identities.py`:
  **ALL PASS** (C1–C43). The explicit-constant bounds of C43b are
  type (c): core \(|f|>e^{-1/6}-e^{1/6}/8>1/4\); collar
  \(|f|\ge e^{-k/12}(1-e^{-5k/6})>e^{-2/3}(1-e^{-5/3})>1/4\) with
  \(2<k=-\log q<8\) (from \(e^2<8\) and \(e^8>2500\)), where the collar
  bounds \(1/1375<q<1/8\) come from the unimodality of \(q\) in
  \(|t_c|^2\in(9/100,1/4)\); in both,
  \(Y_0<1+\sqrt{1/8}<2\). C43e gives exact witnesses at \(|t_c|=1/5\)
  and \(2/5\), on both sides of the selector and at both signs, with a
  rational parametrization of the selector angle and \(|X|<3/5\). C43f,
  numerical observation (d), sample size 8 at 120 digits: projector
  difference \(\le8.6\cdot10^{-122}\), deck residual
  \(\le1.7\cdot10^{-120}\), and measured \(\Delta x_4=+1\) south, \(-1\)
  north. Measured margins (type (d)) beside the proven
  \(|f|>1/4\) and \(Y_0\in(1/4,2)\): C43f gives \(\min|f|=0.6335\) and
  \(Y_0\in[0.6356,0.8663]\); the reviewer's scan gives \(|f|\ge0.578\) and
  \(Y_0\in[0.580,0.996]\). An independent reviewer subagent audited C43.
  It re-derived every claim against `cusp.py` (including the deck
  directions, consistent with C41's chart maps at the corner) and ran:
  - 1568 forced-chart production cases (both cusps; 16 radii in
    \([0.1,0.495]\) including the corner; seven values of \(\Arg t_c\));
  - seam crossings at \(\varepsilon=10^{-20}\) and \(10^{-40}\);
  - the reduced-count production fold across the seam at six points;
  - a bounds scan;
  - 17 mutations of the script, all of which fail.
  Worst residuals were \(3.5\cdot10^{-80}\) (deck) and \(1.1\cdot10^{-81}\)
  (projector). Verdict: sound, with minor caveats. Fixed in this commit:
  - the collar \(q\) range is now proved rather than assumed;
  - tautological sub-checks and a dead variable are removed, and the
    lower-wedge check of the witnesses is added;
  - the seed-cancellation wording, the \(|p|=1\) premise and the C43b
    label are corrected;
  - measured margins are recorded, and the C41 paragraph now points to
    C43.
  The full script was rerun after these fixes: ALL PASS.
  `verify_construction.py`: `all checks
  passed`; `sphere_products.py`: `Dimension 8 check passed`.
  `git diff --check`: pass. LaTeX: **not run** (no manuscript change).
  The legacy numerical safety net was not run (no construction block
  changed).
- **Open points:** the restored-output conclusion rests on a smooth
  final fold vector, as for C40–C42. The remaining T16 work is the
  collar/pure-middle joins beyond the N13 and C40 chambers, then the
  all-branch six-real circuit with global determinant and overlap
  bounds. The global \(J_{\mathrm{cl}}\) statement still requires human
  review.
