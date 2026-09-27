# Read-only diagnostic after the first 10k CS continuation

The reported `continue_stages2_s107/run/final.pt` reduced fixed-validation
Nijenhuis RMS from 1.54069 to 0.46033. Fresh RMS was 0.47060, but adversarial
search found 2.11663 (versus fixed-validation maximum 0.92992). At sampled
points, the maximum norm of J rose from 3.42 to 10.70 and the maximum
conjugation condition number from 5.23 to 12.16. ACS error remained near 1e-14.
This motivates stronger diagnostics before another training allocation.

## What runs

`extra/diagnose_cs.py` and `configs/ambient/s4xs1xs1_cs/diagnostic.yaml`:

- Load the learned CS and the same full-audited, polar-corrected ACS base.
  Verify checkpoint identities and input hashes; freeze parameters, retain
  coordinate derivatives. No optimization of model weights occurs.
- Seed 920107: 16,384 fresh uniform samples; 5,120 structured points (S4 poles
  crossed with both circle grids, and nearby perturbed points).
- Three independent searches, seeds 921107/922107/923107, each selecting 128
  starts from 4,096 points and taking 100 ascent iterations. All endpoints are
  retained, not just the single maximum. Values/input gradients use chunks of
  16 to limit peak memory. Every endpoint's residual is nondecreasing through
  the line search.
- At the *same points*, record Nijenhuis norm, floor-free ACS error, tangency,
  conjugation condition number and Frobenius norm of J. Report the worst 16
  Nijenhuis witnesses with coordinates and their associated diagnostics.
- Save all uniform/structured/adversarial coordinates and per-point measurements
  in `diagnostic_points.pt`, including group labels and checkpoint fingerprint.
  If later used for replay, they become training data and cannot serve as final
  independent validation points.
- `report.json` is an exploratory diagnostic, **not a certificate**. The maximum
  found is a lower bound on the global supremum. No integrability acceptance flag
  is granted, even if numerical values happen to be small.

The default application guard is three hours; PBS requests four hours/10 GB.
There is no training dependency: the CS checkpoint already exists. Runtime has
not been established on the remote learned checkpoint. SIGTERM/SIGINT/application
timeout write an `incomplete` report with completed sample groups and exit 75.
A hard kill can leave no report. This diagnostic does not resume; preserve a
partial report and use a new output directory if another allocation is needed.

## Publish

From the local repository:

```bash
git add extra/diagnose_cs.py extra/test_diagnose_cs.py \
  configs/ambient/s4xs1xs1_cs/diagnostic.yaml docs/s4xt2_cs_diagnostic.md
git commit -m "Add read-only CS witness and conditioning diagnostics"
git push origin s4t2
scp scratch/submit_S4xS1xS1_CS_diagnostic.sh \
  ehirst@lovelace:/home/lovelace/proj/proj925/ehirst/nijenhuis/
```

## Submit on Lovelace

```bash
(
set -euo pipefail
cd /home/lovelace/proj/proj925/ehirst/nijenhuis/github
git pull --ff-only origin s4t2
git cat-file -e HEAD:extra/diagnose_cs.py
cd ..
BASE="$PWD/outputs/s4xs1xs1_acs/polish3_s307/polish_trilin_uniform/s307"
CS="$PWD/outputs/s4xs1xs1_cs/continue_stages2_s107/run/final.pt"
OUT="$PWD/outputs/s4xs1xs1_cs/diagnostic_stages2_s107"
test -f "$CS"
test ! -e "$OUT"
HASH=$(sha256sum "$CS" | awk '{print $1}')
[[ "$HASH" =~ ^[0-9a-f]{64}$ ]]
mkdir "$OUT"
qsub -v "REPO_DIR=$PWD/github,BASE_DIR=$BASE,BASE_AUDIT=$BASE/final.audit.full.seed900001.json,CHECKPOINT=$CS,SOURCE_SHA256=$HASH,OUT=$OUT" \
  submit_S4xS1xS1_CS_diagnostic.sh
)
```

Read the report once finished:

```bash
cat /home/lovelace/proj/proj925/ehirst/nijenhuis/outputs/s4xs1xs1_cs/diagnostic_stages2_s107/report.json
```

Send that JSON for review; keep the binary `diagnostic_points.pt` on the HPC.
The job log is `logs/S4T2_cs_diagnostic_JOBID.log`. The existing campaign-summary
script reads training results, not this diagnostic format; no summary job is
needed. No additional training is submitted by this script. Matched uniform
versus replay continuations are the next decision after reviewing the witnesses.

Local tests use synthetic ACS/CS fixtures, including actual Nijenhuis input
gradients, monotone sphere-constrained ascent, paired metrics at saved witnesses,
immutable inputs, hash/base rejection, interruption and PBS exit-code handling.
They do not constitute a local rerun of the remote learned-field diagnostic.
