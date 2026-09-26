# Reproduce the learned ACS; begin bounded CS experiments

After reviewing the completed one-/two-exponential pilots, use the
[checked continuation handoff](s4xt2_cs_continuation.md) for longer CS runs.

## Validated ACS artifact

The user-supplied full audit of
`polish3_s307/polish_trilin_uniform/s307/final.pt` passed all configured numerical
checks: fresh raw RMS 9.13759e-5, maximum found 2.12031e-4, corrected identity
maximum 2.69093e-15. It used 100,000 uniform points, 5,120 structured points and
three 1,024-start/200-step adversarial searches. This is numerical validation,
not a rigorous uniform bound or an integrability result. The remote weights
have not been copied here or independently rerun locally.

SHA256: `687b3164cf8aa91a202333ba8b5c4aeb0f5aa936f19d629c846990df6e2a0e19`.

## Production ACS recipe

Use `extra/acs_recipe.py` with
`configs/ambient/s4xs1xs1_acs/reproduce.yaml`. Historical pilot configurations
and the generic trainer's defaults are intentionally unchanged. The new recipe
pins every current ACS Config field and reproduces the successful sequence:

| Stage | Seed for root seed 7 | Steps | Batch | Adam schedule |
|---|---:|---|---:|---|
| Initial polynomial | 7 | 5,000 | 192 | 5e-3 to 1e-5 |
| Initial full refinement | 7 | 5,000 | 192 | 1e-3 to 1e-5 |
| Polish 1 | 107 | 10,000 | 768 | 1e-5 to 1e-6 |
| Polish 2 | 207 | 10,000 | 768 | 1e-5 to 1e-6 |
| Polish 3 | 307 | 10,000 | 768 | 1e-5 to 1e-6 |

All stages use float64, one thread, uniform samples, full ambient 9x9 outputs,
`all2+trilin` features and a 4x256 MLP. The initial stage retains its gated,
rollback-protected 200-iteration L-BFGS budget: the observed seed-7 run rejected
its first five iterations. Polishing disables L-BFGS. Each polish loads the
previous **best** `final.pt` and starts a new optimizer/scheduler and sampling
stream. No oracle weights, labels, or analytic reference enter training.

This is not equivalent to one 35,000-step refinement. Historical best-checkpoint
selection also discarded some late updates (s107 selected step 9,000; s207 step
8,500). The recipe repeats that selection logic rather than forcing last states.

Run from the repository:

```bash
python extra/acs_recipe.py --config configs/ambient/s4xs1xs1_acs/reproduce.yaml \
  --seed 7 --out /absolute/new/acs_recipe_s7
```

Stages are saved as `00_s7`, `01_s107`, `02_s207`, `03_s307`. The final candidate
path is recorded in `recipe_result.json`. A fresh seed 8 produces 8/108/208/308.
It tests repeatability; convergence is not guaranteed. Even seed 7 need not be
bitwise identical across PyTorch/BLAS versions. The known run used PyTorch
2.9.1 on Lovelace; local tests use 2.5.1. To obtain the *same field*, preserve
and load the validated checkpoint rather than retraining.

The whole invocation has a three-hour application guard by default. An
interrupted stage saves its existing Runner state and exits 75. Resubmit the
same config, root seed and OUT with `--resume`; completed stages are checked
against recorded hashes and skipped. An independent full audit is still
required for each new candidate. Training completion does not grant acceptance.

## CS: use the validated candidate now, not a new ACS run

`extra/cs_pilot.py` checks a passing **full** ACS audit and its checkpoint hash
before starting the existing ambient CS trainer. It records the base identity
in both manifest and final CS checkpoint. It reconstructs the learned field K,
freezes its parameters (not its coordinate derivatives), and uses the same
differentiable polar correction as the audit to obtain J_base.

The network receives the nine ambient coordinates and predicts one or two
9x9 matrix fields. For B = P net(x) P, it forms
`J = exp(B) J_base exp(-B)` (or two successive conjugations). The output layer
is initialized to zero, so initially J equals J_base. The loss is the existing
round-metric mean squared Nijenhuis norm. Neither orthogonality nor a Kähler
condition is imposed. Every smooth exponential deformation stays homotopic to
its base; this is a search from this ACS, not an exhaustive search over all
topological classes. No existence or convergence guarantee for CS follows from
the ACS audit.

Start with `configs/ambient/s4xs1xs1_cs/smoke.yaml`: 20 steps, one exponential,
2x128 MLP, batch 16, Adam 1e-4 to 1e-5. Review finite losses, ACS residuals,
runtime and Nijenhuis measurements before submitting the optional 1,000-step
pilot (batch 32). Compare one and two exponentials with `--stages 1/2`. These
are proposed pilot hyperparameters, not empirically validated CS settings.
No adaptive/tail sampling or L-BFGS is enabled in this initial comparison.

The pilot wrapper requires a fresh output directory. It does not silently use
the legacy trainer's automatic resume. A wall-stopped run is marked incomplete,
its candidate is named `partial.pt`, and the wrapper exits 75. The legacy wall
guard covers training, not all initial/final diagnostics; the PBS allocation is
the hard limit. This is a bounded smoke/pilot path, not yet a production long-run
CS workflow. Do not extend it to 80,000 steps until profiling and review.

Review `result.json`, `progress.jsonl`, and `cs_pilot_manifest.json`:

- Initialization at the corrected ACS and finite training/validation losses.
- Actual ACS identity residual `dense.acs_floorfree_sup` remains small
  (target <=1e-10, not just an
  `exact_acs_by_construction` label). Conjugation can amplify floating-point
  errors when ill-conditioned. Legacy `acs_sup` retains its 1e-10 reporting
  floor; the new diagnostic removes it without changing training losses.
- Dense and adversarial Nijenhuis errors, not only the training mean.
- No wall stop or nonfinite batches. A completed pilot is **not** a CS success.

Later long runs need independent fixed-set comparisons, condition monitoring,
CS-specific held-out validation, and checked restart/best-candidate handling.
ACS reproducibility on seeds 8 and 9 remains useful but need not prevent this
bounded CS smoke from using the already validated seed-7 lineage.

## Lovelace submission

Publish the new tracked files and copy the two new `scratch/` scripts separately
(scratch remains gitignored):

```bash
git add extra/acs_recipe.py extra/cs_pilot.py extra/campaign_unit.py \
  extra/test_acs_recipe_cs_pilot.py extra/test_acs_pbs.py extra/test_acs_training.py \
  configs/ambient/s4xs1xs1_acs/reproduce.yaml configs/ambient/s4xs1xs1_cs \
  docs/s4xt2_reproduce_and_cs.md docs/s4xt2_acs_hpc.md
git commit -m "Preserve successful ACS recipe and add audited CS pilots"
git push origin s4t2
scp scratch/submit_S4xS1xS1_ACS_reproduce.sh scratch/submit_S4xS1xS1_CS_v2.sh \
  ehirst@lovelace:/home/lovelace/proj/proj925/ehirst/nijenhuis/
```

On Lovelace:

```bash
cd /home/lovelace/proj/proj925/ehirst/nijenhuis/github
git pull --ff-only origin s4t2
git cat-file -e HEAD:extra/cs_pilot.py
git cat-file -e HEAD:configs/ambient/s4xs1xs1_acs/reproduce.yaml
cd ..
BASE="$PWD/outputs/s4xs1xs1_acs/polish3_s307/polish_trilin_uniform/s307"
qsub -v "REPO_DIR=$PWD/github,BASE_DIR=$BASE,BASE_AUDIT=$BASE/final.audit.full.seed900001.json,RECIPE=smoke" \
  submit_S4xS1xS1_CS_v2.sh
```

No dependency is needed: the base already exists and has passed its audit.
Do not use the old `submit_S4xS1xS1_CS.sh` (80k steps, unguarded audit handoff).
The new job prints its output path and logs under `logs/S4T2_cs_JOBID.log`.
Send its three JSON/JSONL files for review before the pilot jobs.

To reproduce ACS independently (optional, separate from the CS smoke):

```bash
qsub -v SEED=7 submit_S4xS1xS1_ACS_reproduce.sh
# Later independent starts: SEED=8 and SEED=9.
# Only after an interruption, using the EXACT original output root:
# qsub -v SEED=7,RESUME=1,OUT=/absolute/original/root submit_S4xS1xS1_ACS_reproduce.sh
```

## Local verification

The focused recipe, ACS, CS handoff, PBS and product-geometry tests pass,
including actual bounded CS optimization, gradients through polar correction
and one/two exponentials, the full-width ACS base, exact interrupted recipe
replay, checksum rejection, and incomplete-run handling. Analytic fixtures are
used only in tests; no local long training or remote learned-weight validation
is claimed. The older `extra/test_all.py` still has the ten previously recorded
float32/float64 failures (nine passes). PyTorch 2.5.1 emits a matrix-exponential
batching performance warning, so HPC profiling is part of the CS smoke review.
