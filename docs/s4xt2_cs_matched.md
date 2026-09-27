# Matched CS continuation after the s107 diagnostic

The reviewed checkpoint has fresh uniform Nijenhuis RMS 0.462289, an adversarial
maximum found of 2.255012, ACS error 1.54e-14, and tangency error 6.79e-15.
This is an almost-complex field with substantial remaining integrability error,
not a numerically validated complex structure. Numerical optimization is not
guaranteed to find an integrable structure.

## Experiment

`extra/cs_matched.py` extends the resumable continuation without changing the
ambient geometry, polar-corrected ACS base, two-exponential conjugation, network
width, Nijenhuis definition, or existing uniform trainer defaults.

Both arms use `configs/ambient/s4xs1xs1_cs/matched.yaml`: seed 207, 2,000 additional
Adam updates, batch 32, float64, one CPU thread, learning rate 1e-4 to 1e-5.
The optimizer is restarted in both arms. These are a controlled pair, not two
independent random-start confirmations and not matched wall-clock experiments.

- **Uniform:** 32 independent uniform points per update.
- **Replay:** the same first 24 uniform points, plus 8 replay points: 4 exact
  witnesses and 4 tangent-perturbed, sphere-retracted points. Initial witnesses
  are adversarial endpoints from `diagnostic_stages2_s107/diagnostic_points.pt`.
  Residual-ranked chord-distance filtering prevents duplicates dominating the
  buffer. Every 250 updates, use a fresh pool and input-gradient ascent to
  refresh up to 256 witnesses, using only training-side samples.
- **Common selection:** 4,096 fresh fixed uniform points plus 64 independent
  source-based adversarial endpoints, generated before any training with a
  separate seed. Select the smallest maximum of (uniform RMS/source RMS,
  challenge maximum/source challenge maximum), subject to numerical ACS,
  tangency and conditioning gates. The starting checkpoint participates, so the
  selected step can be zero if no checkpoint improves both metrics. These
  repeated selection sets are validation data, not an independent final audit.
- **Monitoring:** retain Nijenhuis quantiles, ACS/tangency errors, conjugation
  condition number and J Frobenius norm. No orthogonality constraint or norm
  penalty has been added. The conditioning gate is a loose numerical guard,
  not a compactness guarantee; growth must be reviewed after the run.
- **Fresh diagnostics:** each PBS job automatically evaluates its selected
  checkpoint using `diagnostic_matched.yaml`: seed 940207, 16,384 uniform,
  5,120 structured, and three 128-start/100-step searches. The same seeds are
  used for both arms. The adversarial endpoints can differ between arms.
  These diagnostics have not supplied gradients or selected checkpoints.

The old diagnostic witnesses are now training data. Do not reuse that old
report as independent validation. None of these finite tests gives a global
supremum bound or an integrability certificate.

## Publication and submission

From the **local checkout**, commit and push the implementation (not scratch):

```bash
cd /Users/ed/Documents/Projects/DG/Nijenhuis/github
git add extra/cs_pilot.py extra/cs_training.py extra/cs_matched.py extra/test_cs_matched.py configs/ambient/s4xs1xs1_cs/matched.yaml configs/ambient/s4xs1xs1_cs/diagnostic_matched.yaml docs/s4xt2_cs_matched.md
git commit -m "Add matched uniform and targeted CS continuations"
git push origin s4t2
scp scratch/launch_S4xS1xS1_CS_matched.sh scratch/submit_S4xS1xS1_CS_matched.sh scratch/submit_S4xS1xS1_CS_matched_summary.sh ehirst@lovelace:/home/lovelace/proj/proj925/ehirst/nijenhuis/
```

On **Lovelace**, paste commands only, without terminal prompt characters:

```bash
cd /home/lovelace/proj/proj925/ehirst/nijenhuis/github
git switch s4t2
git pull --ff-only origin s4t2
cd /home/lovelace/proj/proj925/ehirst/nijenhuis
bash launch_S4xS1xS1_CS_matched.sh
qstat -u "$USER"
```

Run the launcher with **bash**, not qsub. It validates nonempty inputs/scripts and
the two reviewed hashes before submitting anything, creates a unique directory
under `outputs/s4xs1xs1_cs/matched_s207_*`, submits two jobs, then submits an
`afterany` summary dependent on both. No manual hold is needed. It prints all
three job IDs, records them in `submission.tsv`, and prints the exact command to
read `comparison.txt`. `comparison.json` retains both training and diagnostic
reports. Missing/failed arms are explicit and cannot be declared a completed
comparable pair. A completed pair still never sets `numerically_validated_cs`.

The scripts remain under ignored `scratch/` locally and must be copied with SCP.
They request serial/10 GB/four hours per training-and-diagnostic job. Each phase
has a separate 90-minute application guard, leaving PBS headroom. Actual runtime
of the new training arm on the learned checkpoint has not been measured locally.

## Failure and restart

- If submission fails, do not blindly rerun the launcher: inspect `qstat` and
  the printed directory's `submission.tsv`. Already submitted jobs are retained;
  a new invocation creates a new experiment, not a resume.
- A cooperative training stop exits 75, saves `state.pt`, and does not run the
  diagnostic. Re-submit **that arm's** PBS script with the same environment and
  `RESUME=1`; the model, optimizer, scheduler, buffer, fixed challenge and both
  training RNG streams are restored. Input hashes and config must still match.
- Completed training cannot be resumed into the same directory. If only its
  diagnostic failed, preserve the training result and rerun the diagnostic alone
  into a new empty directory, using the matched diagnostic config. Do not repeat
  training or overwrite the old partial diagnostic.
- Stop checks happen between updates and inside input-gradient searches;
  validation/final evaluation complete their current operation before stopping.
  A hard kill may lose updates since the last atomic checkpoint (every 100).
- The summary is a read-only comparison of model artifacts. It writes only
  `comparison.json`/`comparison.txt`; it does not choose a winner or launch more
  optimization. Once the final diagnostics inform subsequent training choices,
  they are development evidence, not an untouched final certification set.

## Local verification

Tests exercise actual small-network CS differentiation and optimization, common
validation and uniform RNG streams, replay/jitter geometry, input hash checks,
diversity, both-arm exact interrupted/resumed trajectories, transactional search
interrupts, tail-aware selection, completion summaries, and mocked PBS submission
and dependency failures. They use synthetic fixture checkpoints, not the remote
learned checkpoint, and do not establish convergence of the production search.
