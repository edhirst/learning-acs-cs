# Checked CS continuation from the two-exponential pilot

The user's first pilots completed 1,000 steps. One exponential reached Nijenhuis
RMS 1.70889 and adversarial maximum found 2.79320 in 3.93 minutes; two reached
1.52273 and 2.53700 in 5.42 minutes. Both retained ACS errors below 5e-15.
These are exploratory measurements on different samples, not CS certificates
or proof of an architectural advantage. The next bounded experiment continues
the two-exponential candidate, without resetting its learned deformation.

## What the new runner does

Use `extra/cs_training.py` with `configs/ambient/s4xs1xs1_cs/continue.yaml`.
Unlike the pilot-only wrapper, it checks the source CS architecture and its
recorded ACS base hash, loads its weights unchanged and records the source hash.
The original full-audited ACS and its differentiable polar correction remain
frozen. The ACS source and audit are checked again before publishing a result.

- 10,000 **additional** Adam updates, two exponentials, 2x128 MLP, batch 32,
  uniform sampling, float64, one thread, fresh sampling seed 107.
- New Adam and cosine schedule: 1e-4 to 1e-5. No L-BFGS or tail sampling.
- 1,024 fixed validation points, evaluated every 250 steps; save the lowest
  validation Nijenhuis RMS among candidates with ACS error <=1e-10 and ambient
  conjugation condition number <=1e6. These are numerical guardrails, not
  global bounds. The initial source candidate participates in selection.
- Monitor Nijenhuis mean/RMS/maximum, floor-free ACS error, condition number
  of G = exp(B2)exp(B1), and Frobenius norm of J.
- Atomic restart state every 100 steps, including optimizer, scheduler,
  training generator, fixed validation points and selected best model.
  SIGTERM/SIGINT stop cooperatively after the current work unit. Per-allocation
  application budget is 180 minutes; PBS is four hours. Initialization and final
  reports are not interrupted internally, so retain this scheduler margin.
- `state.pt` is the last training state; `best.pt` is the selected candidate.
  `final.pt` is published only after all requested updates and final reporting.
  An incomplete run exits 75 and can be resumed with identical numerical config.
- Final diagnostics use 2,048 fresh points and 32 adversarial starts/40 steps.
  They are exploratory, not a full independent CS audit. No chart report is
  generated here (`n_chart`/`n_initial` are unused compatibility fields inherited
  from the pilot configuration schema).

`numerically_validated_cs` remains false. Small algebraic ACS error is not
integrability. Conjugation also stays within the reference's homotopy class.
Neither finite searches nor this one lineage resolve CS existence.

## Publish and submit

Locally, commit and push these tracked files on `s4t2`:

```bash
git add extra/cs_training.py extra/test_cs_training.py \
  configs/ambient/s4xs1xs1_cs/continue.yaml docs/s4xt2_cs_continuation.md \
  docs/s4xt2_reproduce_and_cs.md
git commit -m "Add checked resumable CS continuation"
git push origin s4t2
scp scratch/submit_S4xS1xS1_CS_continue.sh \
  ehirst@lovelace:/home/lovelace/proj/proj925/ehirst/nijenhuis/
```

On Lovelace, pull `s4t2`, then submit from the parent `nijenhuis` directory:

```bash
(
set -euo pipefail
cd /home/lovelace/proj/proj925/ehirst/nijenhuis/github
git pull --ff-only origin s4t2
git cat-file -e HEAD:extra/cs_training.py
cd ..
REPO="$PWD/github"
BASE="$PWD/outputs/s4xs1xs1_acs/polish3_s307/polish_trilin_uniform/s307"
SOURCE="$PWD/outputs/s4xs1xs1_cs/pilots_audited_s7/stages2_s7/final.pt"
ROOT="$PWD/outputs/s4xs1xs1_cs/continue_stages2_s107"
OUT="$ROOT/run"
test -f "$SOURCE"
HASH=$(sha256sum "$SOURCE" | awk '{print $1}')
[[ "$HASH" =~ ^[0-9a-f]{64}$ ]]
mkdir "$ROOT"
JOB=$(qsub -v "REPO_DIR=$REPO,BASE_DIR=$BASE,BASE_AUDIT=$BASE/final.audit.full.seed900001.json,INIT_FROM=$SOURCE,SOURCE_SHA256=$HASH,OUT=$OUT" submit_S4xS1xS1_CS_continue.sh)
printf 'CS job: %s\n' "$JOB"
[[ "$JOB" =~ ^[0-9]+([.][A-Za-z0-9._-]+)?$ ]]
printf 'role\tjob_id\toutput\ncs_continue\t%s\t%s\n' "$JOB" "$OUT" > "$ROOT/submission.tsv"
SUMMARY=$(qsub -W "depend=afterany:$JOB" -v "REPO_DIR=$REPO,SUMMARY_ROOT=$ROOT" submit_S4xS1xS1_ACS_summary.sh)
printf 'summary\t%s\t%s\n' "$SUMMARY" "$ROOT" >> "$ROOT/submission.tsv"
printf 'Summary: %s/summary.%s.txt\n' "$ROOT" "$SUMMARY"
)
```

Do not rerun this block if only the summary submission fails; the CS job may
already be queued. Inspect `submission.tsv` and `qstat` first. After completion,
send `run/result.json` and `run/progress.jsonl` from the new root.

Only after an interruption (not after completion), resubmit the same output:

```bash
cd /home/lovelace/proj/proj925/ehirst/nijenhuis
BASE="$PWD/outputs/s4xs1xs1_acs/polish3_s307/polish_trilin_uniform/s307"
OUT="$PWD/outputs/s4xs1xs1_cs/continue_stages2_s107/run"
qsub -v "REPO_DIR=$PWD/github,BASE_DIR=$BASE,BASE_AUDIT=$BASE/final.audit.full.seed900001.json,OUT=$OUT,RESUME=1" submit_S4xS1xS1_CS_continue.sh
```

The summary attached to the first allocation does not wait for this replacement
job; submit another summary if resuming. Keep the original CS pilot and ACS
artifact unchanged. Local tests use analytic ACS fixtures, not the inaccessible
remote learned weights; full HPC continuation has not been run locally.
