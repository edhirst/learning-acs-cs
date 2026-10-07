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
