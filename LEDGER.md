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
