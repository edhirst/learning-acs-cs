# `agentic_ledger`: agent-assisted revision of the manuscript

This branch holds only the manuscript *An Explicit Complex Structure on S⁴ × T²* and the
infrastructure for improving it with AI agents deployed by any of the collaborators.
It shares no files with the code branches of this repository.

| path | contents |
|---|---|
| `manuscript_draft/` | the manuscript: `main.tex`, `references.bib`, `utphys.bst`, `figures/`, and `code/` (the verification scripts of Appendix A) |
| `AGENTS.md` | binding instructions for every agent. Read automatically by tools that support `AGENTS.md` (e.g. Codex and the GitHub Copilot coding agent) |
| `CLAUDE.md` | makes Claude Code load `AGENTS.md` automatically |
| `LEDGER.md` | append-only record of every change: when, by whom, why, and how it was verified |
| `.gitignore` | keeps LaTeX build products and local files out of the repository |

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

## Reviewing and undoing changes

- History of the manuscript: `git log --oneline -- manuscript_draft/`. One change:
  `git show <hash>`. Who last changed a line: `git blame manuscript_draft/main.tex`.
- Everything since the baseline: `git diff <baseline> -- manuscript_draft/`. The
  baseline is the commit of ledger entry `2026-10-01-TSG-1`; find it with
  `git log --grep 2026-10-01-TSG-1`. To give it a name:
  `git tag ms-baseline <hash> && git push origin ms-baseline`.
- To undo a change, run `git revert <hash>` and add a `revert` entry to the ledger.
  Never rewrite history that has already been pushed.

## Rules in brief

- Agents edit only `manuscript_draft/`, and append a ledger entry in the same commit.
- Every change must be verified and must be an improvement. When unsure, an agent
  records a flag in the ledger instead of editing.
- Mathematical statements are never changed silently, and citations are never invented.
- The LaTeX build and the two verification scripts must stay at least as good as the
  baseline recorded in the ledger.
- No force-pushes, and no rewriting of history.
- Only the human collaborators change `AGENTS.md`, `CLAUDE.md`, `README.md` and
  `.gitignore`, in a commit with an `infra` ledger entry.
