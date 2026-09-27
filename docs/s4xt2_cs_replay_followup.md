# Lower-rate replay-mixture comparison from the selected s207 replay field

The first targeted run ended at its selected step 2,000 with fixed-validation
RMS 0.443039 and challenge maximum 1.084134. A fresh diagnostic found RMS
0.443675 and maximum 1.569371. The last 500 updates still improved RMS, but
challenge error was not monotone. This motivates another bounded experiment,
not a claim of integrability or guaranteed convergence.

## Exact inputs and experiment

Both arms start at:

```
/home/lovelace/proj/proj925/ehirst/nijenhuis/outputs/s4xs1xs1_cs/matched_s207_hs8qhz51/replay/final.pt
```

Witnesses come from that run's `diagnostic/diagnostic_points.pt`, **not** the older
`diagnostic_stages2_s107` directory. `extra/cs_followup_inputs.py` uses only the
Python standard library on the login node. Before submitting, it checks the
completed training status, seed and selected step; the recorded parent hash;
the completed diagnostic's checkpoint and points-file hashes; and the shared
audited ACS base and audit hashes. The current checkpoint hash is read and
verified on the HPC, not guessed from the parent hash in `result.json`.

The two YAMLs differ only in `sampling.target_count`:

| Arm | Uniform / targeted per batch | Configuration |
|---|---|---|
| replay25 | 24 / 8 | replay_followup25.yaml |
| replay12p5 | 28 / 4 | replay_followup12p5.yaml |

Common settings: seed 307, 5,000 additional updates, fresh Adam, cosine LR
2e-5 to the existing 1e-5 floor, batch 32, float64, two exponential stages,
hidden widths [128,128], same fixed validation/challenge protocol and replay
refresh every 250 steps. This is a new continuation, not an exact optimizer
resume of the completed run. Both arms consume the same uniform stream; replay
streams diverge because their sample counts differ. Architecture, geometry,
Nijenhuis loss and checkpoint selection rules are unchanged.

Each job automatically runs a post-selection diagnostic using seed 941307,
16,384 uniform points, 5,120 structured points and three adversarial searches
with 128 starts and 100 steps. This uses new random seeds, not the new training
witnesses. Structured grids include deterministic points as before. Each phase
has a guard (120 minutes training, 90 minutes diagnostic) within a four-hour,
10 GB serial allocation. The summary job waits for both, even if a job fails.

The comparison checks that the only intended configuration difference is 8
versus 4 targeted points. It requires matching source/base/initial validation,
and completed diagnostics whose hashes match the selected checkpoints. Missing
or failed arms cannot form a completed comparable pair. No automatic winner is
chosen and no finite diagnostic is labelled an integrability certificate.

## Publish locally

```bash
cd /Users/ed/Documents/Projects/DG/Nijenhuis/github
git add extra/cs_matched.py extra/cs_followup_inputs.py extra/test_cs_replay_followup.py configs/ambient/s4xs1xs1_cs/replay_followup25.yaml configs/ambient/s4xs1xs1_cs/replay_followup12p5.yaml configs/ambient/s4xs1xs1_cs/diagnostic_replay_followup.yaml docs/s4xt2_cs_replay_followup.md
git commit -m "Add lower-rate CS replay-mixture follow-up"
git push origin s4t2
```

After the push succeeds, copy **all three** scripts; the two existing wrappers
have changed to accept the new configurations and summary experiment:

```bash
scp scratch/launch_S4xS1xS1_CS_replay_followup.sh scratch/submit_S4xS1xS1_CS_matched.sh scratch/submit_S4xS1xS1_CS_matched_summary.sh ehirst@lovelace:/home/lovelace/proj/proj925/ehirst/nijenhuis/
```

## Run on Lovelace

```bash
cd /home/lovelace/proj/proj925/ehirst/nijenhuis/github
git switch s4t2 && git pull --ff-only origin s4t2
```

After that succeeds:

```bash
cd /home/lovelace/proj/proj925/ehirst/nijenhuis
bash launch_S4xS1xS1_CS_replay_followup.sh
qstat -u "$USER"
```

Run the launcher with **bash**, not qsub. It prints three job IDs and the exact
command to read `comparison.txt`. Outputs are in a fresh directory:

```
outputs/s4xs1xs1_cs/replay_followup_s307_XXXXXXXX/
  submission.tsv
  replay25/
  replay12p5/
  comparison.txt
  comparison.json
```

Do not blindly re-run the launcher after a partial submission failure: it creates
a new experiment each time. Inspect the printed `submission.tsv` and `qstat`
first. Training exit 75 leaves a resumable state but skips the diagnostic.
For a training resume, use the same arm/config/output environment with RESUME=1.
If training completed and only the diagnostic failed, preserve it and rerun
only the diagnostic into a new empty directory; do not retrain the completed arm.

Send `comparison.txt` after completion. Review average and tail errors together
with conditioning and J norms before deciding on any further allocation.
Scripts remain gitignored and are transferred by SCP, not Git.
