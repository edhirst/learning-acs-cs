# Instructions for agents working on the manuscript (branch `agentic_ledger`)

This branch contains the manuscript *An Explicit Complex Structure on S⁴ × T²*
(E. Hirst, T. Schettini Gherardini, T. S. R. Silva, J. Stein) in `manuscript_draft/`.
Agents deployed by any of the collaborators work here to improve its mathematical
content and the quality of its exposition. These instructions bind every agent.
Only the human collaborators may change them.

## 0. The prime rule

**The manuscript may only get better.** Commit a change only after you have verified
that it is correct and that it improves the text. If you are not certain, do not edit.
Record the issue as a `flag` entry in `LEDGER.md` and leave the manuscript unchanged.
A missed improvement costs nothing. An error that nobody notices can cost the paper.

## 1. Scope

1. You may modify only files inside `manuscript_draft/`. You must also append your
   entry to the end of `LEDGER.md`.
2. Do not modify `AGENTS.md`, `CLAUDE.md`, `README.md`, `.gitignore`, or any earlier
   ledger entry. If you think one of them needs to change, tell the collaborator who
   deployed you.
3. Do not create, delete, rename or move files without an explicit request from your
   collaborator. The same applies to restructuring the document (reordering, merging
   or splitting sections, renumbering results). Record any such request in your
   ledger entry.
4. Never commit build products (the compiled PDF, `.aux`, `.log`, `.bbl`, …) or
   caches. Stage files by explicit path. Never use `git add -A`, `git add .` or
   `git commit -a`.
5. The manuscript is unpublished. Do not upload or paste it, or substantial parts of
   it, to any external service (gists, issue trackers, forums, web forms). You may
   search the web for references, but do not paste the manuscript's text into searches.

## 2. Before you start

1. Synchronise:
   `git fetch origin && git switch agentic_ledger && git pull --ff-only origin agentic_ledger`.
2. Read this file, `README.md` and the whole of `LEDGER.md`. The ledger is the
   project's memory: past changes, decisions and open flags. Do not redo, undo or
   contradict them without a new argument (§4).
3. Read the passage you will work on, together with everything it depends on:
   definitions, notation, the earlier results it uses, and the *Conventions*
   paragraph at the end of §1. Mathematical text is global. A change to a
   definition, a sign, a normalisation or a statement propagates to every place that
   uses it, so search the whole source for every occurrence of the symbols and
   `\label`s involved.
4. Make sure you understand the task your collaborator gave you. If it is ambiguous,
   ask before editing.
5. Build the manuscript and run the checks of §3.3 **before** editing. Then you know
   the starting state and can tell which problems, if any, are yours.

## 3. Standards for every change

### 3.1 Mathematical rigour

1. **Verify everything you write.** Re-derive every formula, sign, index, factor,
   inequality and quantifier that you add or change. Wherever possible, confirm it
   independently by a special case, a dimension or degree count, a known example or a
   symbolic computation. Record in your ledger entry how you verified it.
2. **No silent statement changes.** Changing the hypotheses, conclusion, quantifiers or
   scope of a theorem, proposition, lemma, corollary or claim is a *statement change*.
   Make one only with a complete justification, use ledger type `statement`, and tell
   your collaborator explicitly. A statement change always needs human review.
3. Do not delete mathematical content (proof steps, remarks, caveats, examples,
   consistency checks) unless it is wrong or redundant. Justify any deletion in the
   ledger.
4. **Keep the dependency bookkeeping exact.** The Main Theorem and Theorem B are
   unconditional. They use only the classical facts (F1)–(F4) (`\label{sec:facts}`)
   and explicit computations. Theorem C and Corollary D depend on the inputs
   (I1)–(I3) (`\label{sec:inputs}`), which rest on 2026 preprints. Never make an
   unconditional result depend on a preprint, and never drop an "Assume the input …"
   clause. State the inputs of every new argument.
5. Keep proved results, heuristic checks, numerical evidence and conjectures
   distinct, and never present one as another.
6. **Follow the manuscript's conventions** (*Conventions*, end of §1): the
   identification ℂ² ≅ ℝ⁴ via (z, w) ↦ (Re z, Im z, Re w, Im w), which gives the
   complex orientation; det_ℝ; ω(u, v) = Im(ū v); the orientation of each Σ_jk by
   (e_j, e_k); and the rest. Signs and orientations carry the main idea of the paper
   (the two cusps have opposite chirality). Treat every sign or orientation change as
   a statement change unless you have re-derived it.
7. Do not add "clearly", "obviously" or "it is easy to see" to cover a step you have
   not checked. A gap in a proof is a legitimate target: fill it with a complete
   argument.

### 3.2 References and the literature

1. **Never invent or "complete" bibliographic data**: citations, BibTeX entries,
   theorem, section or page numbers inside a cited work, arXiv identifiers, DOIs,
   author lists, years. Add or change a reference only after checking it against the
   source itself. Otherwise flag it.
2. Several cited works (`AC26`, `Engel26`, `GFL26`, `LLZ26`, …) are 2026 preprints
   that may postdate your training data. Do not "correct" what the manuscript says
   about them, or about which products of spheres are known to be complex, on the
   strength of your prior knowledge. If you believe a statement about the literature
   is wrong, flag it with your evidence.
3. Keep `references.bib` consistent: every `\cite` key must be defined, and every
   entry must be cited. At the baseline there are 41 keys, all defined and used.

### 3.3 The build and the checks must stay green

Run these commands from `manuscript_draft/`:

| check | command | baseline (2026-10-01) |
|---|---|---|
| LaTeX | `latexmk -pdf -interaction=nonstopmode main.tex` | 34 pages; 0 errors; 0 undefined references or citations; 0 BibTeX warnings; the only warnings are 37 hyperref "Token not allowed in a PDF string" |
| symbolic checks | `python3 code/verify_construction.py` | ends with `all checks passed` (sympy 1.14) |
| products of spheres | `python3 code/sphere_products.py` | ends with `Dimension 8 check passed` |

1. After your change, all three must be at least as good as the baseline. That means
   no errors, no undefined references or citations, no failing `assert`, and no
   warning of a kind absent from the baseline.
2. If you cannot run a check (for example, no TeX or Python is available), write
   "not run" in your ledger entry. **Never report a check that you did not run.**
3. The scripts in `code/` verify specific formulas of the paper; Appendix A
   (`\label{app:checks}`) lists which. If you change a formula covered by a check,
   update the script and the appendix in the same commit, and run the script. If you
   change a script, it must still pass. Never weaken an `assert` to make it pass.
4. Remove build products with `latexmk -C`, or simply leave them alone:
   `.gitignore` excludes them.

### 3.4 Form of the edit

1. One logical change per commit. Keep diffs minimal. Do not re-wrap, re-indent or
   reformat text you are not otherwise changing, and make no stylistic edits outside
   your task.
2. Preserve the authors' voice and the existing structure. Prefer the smallest edit
   that achieves the improvement.
3. Do not rename existing `\label` keys or macros, and do not add packages. Define a
   new macro only if it is necessary, once, in the preamble.
4. Keep the build self-contained: no shell escape and no external tools.

## 4. Working alongside other agents and collaborators

1. Several agents may work at the same time. Before editing a passage, check its
   recent history (`git log -p -- manuscript_draft/main.tex`) and the ledger.
2. Do not revert or rewrite a change recorded in an earlier entry unless it is
   demonstrably wrong. Then make a `revert` entry that cites the earlier entry's ID
   and gives the argument. A difference of taste is no reason to revert; flag it
   instead.
3. Never rewrite or delete shared history. That rules out `git push --force`,
   resetting or rebasing commits that are already on the remote, `--amend` on pushed
   commits, and deleting or moving branches or tags.
4. **Push only with your collaborator's approval**: after they have reviewed the
   commit, or if they have explicitly authorised you to push without review. The
   command is `git push origin agentic_ledger`.
5. If the push is rejected because the remote has moved, run
   `git pull --rebase origin agentic_ledger`. Then re-read the passages you changed
   in their new context and re-run the checks before pushing again.
6. Resolve conflicts by keeping every other collaborator's change. If two changes
   cannot be reconciled with both intents preserved, stop and ask. A conflict in
   `LEDGER.md` caused only by parallel appends is resolved by keeping both entries,
   with yours last.
7. The human collaborators have the final word. When in doubt, ask the one who
   deployed you.

## 5. Recording your change in `LEDGER.md`

1. Every commit that touches `manuscript_draft/` must append exactly one entry to
   `LEDGER.md`, in the same commit, using the template there.
2. A flag (an issue you found but did not fix) is also an entry, committed on its own.
3. Never edit or delete earlier entries. To correct one, add a new entry that cites it.
4. The commit message starts with the entry ID in brackets, followed by a one-line
   summary, for example `[2026-10-02-EH-1] Fix the sign in the proof of Lemma 2.2`.
   Then `git log --grep '<ID>'` finds the commit.

## 6. Checklist before every commit

- [ ] The change is in scope (§1) and is what your collaborator asked for.
- [ ] Every formula, sign and statement you touched has been re-derived and verified
      (§3.1). There are no silent statement changes.
- [ ] No unverified citation and no unverified claim about the literature (§3.2).
- [ ] Notation, conventions, `\label`s and dependency markers agree with the rest of
      the paper.
- [ ] The build and both scripts are at least as good as the baseline, or "not run" is
      recorded (§3.3).
- [ ] The diff is minimal. Only intended files are staged; check `git status` and
      `git diff --cached`.
- [ ] Your ledger entry is appended in the same commit, and the commit message starts
      with its ID.
- [ ] You synchronised with the remote, and you will push only with your
      collaborator's approval and never with force.

Git keeps every version, so any change can be undone, but only once someone notices
that it was wrong. Review your own change as strictly as you would want a referee to
review the paper.
