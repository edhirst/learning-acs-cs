# `agentic_ledger`: agent-assisted revision of the manuscript

This branch holds only the manuscript *An Explicit Complex Structure on S⁴ × T²* and the
infrastructure for improving it with AI agents deployed by any of the collaborators.
Since 2026-10-07 it also hosts the **symbolic campaign** (see `SYMBOLIC_CAMPAIGN.md`), a shared workspace for the collaborators and their agents. Apart from the two vendored modules listed below, it shares no files with the code branches of this repository.

| path | contents |
|---|---|
| `manuscript_draft/` | the manuscript: `main.tex`, `references.bib`, `utphys.bst`, `figures/`, and `code/` (the verification scripts of Appendix A) |
| `AGENTS.md` | binding instructions for every agent. Read automatically by tools that work in the repository and support `AGENTS.md` (e.g. OpenAI Codex, also when launched from ChatGPT, and the GitHub Copilot coding agent). A chat assistant such as ChatGPT must be given the file explicitly (see below) |
| `CLAUDE.md` | makes Claude Code load `AGENTS.md` automatically |
| `LEDGER.md` | append-only record of every change: when, by whom, why, and how it was verified |
| `.gitignore` | keeps LaTeX build products and local files out of the repository |
| `SYMBOLIC_CAMPAIGN.md` | **entry point of the symbolic campaign** (since 2026-10-07): the goal (a rigorous, constructive, symbolic-elementary proof that the round S⁴ × T² carries the paper's complex structure, clear and easy to review), the state of the art, the task list, and the rules for contributing |
| `docs/s4t2_explicit_diffeo.md` | the construction log of the explicit diffeomorphism F: Y → S⁴ × T² (§12 = the verified numerical construction; §13 = the agreed proof scheme) |
| `docs/s4t2_symbolic_programme.md` | the symbolic programme: block-by-block specification of F and G = F⁻¹, lemma list with status, certificates, open conditions |
| `extra/s4t2_explicit/` | the code of the construction (numpy reference, torch implementation, numerical checks, fits) and `sym/`, the sympy certificates |
| `src/geometry/ambient.py`, `src/loss/ambient_loss.py` | two modules vendored verbatim from the code branch `no_overlaps` (commit `dcf2bab`): the ambient Nijenhuis tensor used by the numerical checks. Frozen on this branch |

## Deploying an agent

1. In any clone of the repository:
   `git fetch origin && git switch agentic_ledger && git pull --ff-only origin agentic_ledger`.
2. Start the agent in the repository root and give it a precise task, for example:
   > Read AGENTS.md and follow it strictly. Task: check the proof of Lemma~\ref{lem:slabs}
   > line by line. Fix only errors you can prove; flag everything else in LEDGER.md.

   If your tool does not read `AGENTS.md` by itself, tell the agent to read it first.
3. Before pushing, review the agent's commit: the diff (`git show`), its new
   `LEDGER.md` entry, and the build and check results it reports. Then push with
   `git push origin agentic_ledger`, or authorise the agent to push. `AGENTS.md`
   forbids agents to push without your approval.

If you use **Codex** (from ChatGPT, the Codex CLI or an IDE), it works directly in the
repository and reads `AGENTS.md`. Proceed as above, on the branch `agentic_ledger`, and
review its diff before anything is pushed or merged.

## Using ChatGPT, or another chat assistant without access to the repository

In a chat window the assistant cannot run git or the checks, so you apply its proposal,
run the checks and make the commit. `AGENTS.md` §7 tells the assistant what to deliver.

1. **Before uploading the unpublished manuscript**, check your data settings. In ChatGPT
   you can, for example, turn off training on your conversations (Settings → Data
   controls), or use a temporary chat.
2. Synchronise: `git fetch origin && git switch agentic_ledger && git pull --ff-only origin agentic_ledger`.
3. In a new conversation, upload `AGENTS.md`, `LEDGER.md` and the files the task needs,
   usually `manuscript_draft/main.tex` (add `references.bib` or the scripts in `code/`
   if the task involves them). Then write, for example:
   > Read AGENTS.md and follow it strictly. You have no access to the repository, so
   > follow its §7. Task: check the proof of Lemma~\ref{lem:slabs} line by line.
   > Propose only changes you can prove, as a unified diff or as exact OLD/NEW pairs,
   > and flag everything else.
4. Apply the proposal.
   - **Diff:** save it as `change.patch`, then run
     `git apply --check change.patch && git apply change.patch`. Chat assistants often
     get the line counts in hunk headers wrong. If git calls the patch corrupt, add
     `--recount` to both commands. If it still fails, ask for OLD/NEW pairs instead.
     Delete `change.patch` afterwards.
   - **OLD → NEW pairs:** make each replacement in your editor. Each OLD must occur
     exactly once in the file.
5. Check with `git diff` that it contains exactly the proposed changes and nothing else.
6. From `manuscript_draft/`, run the build and both scripts (`AGENTS.md` §3.3).
7. Complete the ledger entry drafted by the assistant. Fill in **Based on** with
   `git rev-parse --short HEAD` (run before committing), and fill in **Verification**
   with what you actually ran and its results. Never record a check that nobody ran.
8. Append the entry to `LEDGER.md`. Commit `LEDGER.md` and the changed files, staged by
   explicit path, with the message `[<ID>] <summary>`. Then push.

## Reviewing and undoing changes

- History of the manuscript: `git log --oneline -- manuscript_draft/`. One change:
  `git show <hash>`. Who last changed a line: `git blame manuscript_draft/main.tex`.
- Everything since the baseline: `git diff <baseline> -- manuscript_draft/`. The
  baseline is the commit of ledger entry `2026-10-01-TSG-1`; find it with
  `git log --grep 2026-10-01-TSG-1`. To give it a name:
  `git tag ms-baseline <hash> && git push origin ms-baseline`.
- To undo a change, run `git revert <hash>` and add a `revert` entry to the ledger.
  Never rewrite history that has already been pushed.

## The symbolic campaign

`SYMBOLIC_CAMPAIGN.md` is the entry point: read it first, then `docs/s4t2_symbolic_programme.md`.
Agents working on the campaign obey `AGENTS.md` §8 (which extends the scope to the campaign
files) and append a ledger entry per commit exactly as for the manuscript. Every mathematical
claim added to the campaign documents must be backed by a script in `extra/s4t2_explicit/sym/`
or by an elementary argument with explicit constants written out in full.

## Rules in brief

- Agents edit only `manuscript_draft/` (manuscript tasks) or the campaign files listed in `AGENTS.md` §8 (campaign tasks), and append a ledger entry in the same commit.
- Every change must be verified and must be an improvement. When unsure, an agent
  records a flag in the ledger instead of editing.
- Mathematical statements are never changed silently, and citations are never invented.
- The LaTeX build and the two verification scripts must stay at least as good as the
  baseline recorded in the ledger.
- No force-pushes, and no rewriting of history.
- Only the human collaborators change `AGENTS.md`, `CLAUDE.md`, `README.md` and
  `.gitignore`, in a commit with an `infra` ledger entry.
