# Long, monitored CS continuation from the selected s307 12.5% replay field

This experiment allows **50,000 additional Adam updates**, ten times the last
5,000-update run. It continues the selected step-4,750 field at:

```
/home/lovelace/proj/proj925/ehirst/nijenhuis/outputs/s4xs1xs1_cs/replay_followup_s307_YIFBe02B/replay12p5/final.pt
```

The latest diagnostic witnesses in that run's `diagnostic/` become training
replay. The previously learned ACS base is unchanged. `cs_followup_inputs.py
--stage s307` verifies the parent lineage, selected step, 12.5% arm, completed
diagnostic hashes and full ACS audit before any submission.

## Training and monitoring

- Seed 407, batch 32, 28 uniform and 4 targeted points; half the targeted points
  are perturbed and retracted. Existing periodic training replay refresh remains.
- Constant learning rate **1e-5**. Adam restarts once from the selected source;
  model/optimizer/scheduler/RNG/buffer state is resumed exactly across allocations.
- Same ambient field, two-exponential architecture, float64 and one CPU thread.
- Fixed uniform validation has 4,096 points. Initial selection challenge uses
  an independent pool of 8,192 points, 128 starts and 80 ascent steps.
- Every 2,500 updates, a separate validation-only search uses 2,048 fresh pool
  points, 64 starts and 60 ascent steps. Its endpoints are added to a growing
  monitoring set, never to the training buffer. Evaluation occurs every 500
  updates. Retained candidates are re-evaluated whenever that set expands.
- Separate candidates are retained: lowest uniform RMS subject to a tail gate,
  and lowest monitored tail subject to an RMS gate. Each gate allows at most
  5% regression relative to the source's initial measured error. Both candidates
  are compared on the same validation/challenge/monitoring points. The monitored
  tail is the largest residual found across those sets, not a global bound.
- Checkpoints are atomic every 100 updates and at allocation boundaries.
  Pending validation searches are recorded: after an interruption, the same
  observation finishes before another training update.

Operational early stops request review (not success or a mathematical obstruction):

- Maximum J norm doubles, or conjugation condition triples, relative to the
  source **on the same fixed uniform/challenge sets**.
- ACS/tangency/numerical gates fail, or the retained candidates cannot satisfy
  the average/tail selection gates on an expanded monitoring set.
- Neither retained RMS nor monitored tail improves by 0.5% over a 10,000-update
  window. Because monitoring grows, this is a conservative compute-budget
  heuristic, not proof of a stationary point or nonexistence.

After a review stop, retained candidates are exported for diagnosis when available;
the remaining training slots become no-ops. The reports clearly say
`stopped_for_review`, not completed training or certified CS. Nonfinite optimization
failures instead leave the last atomic checkpoint and fail the allocation.

## PBS allocation chain

The launcher submits five sequential training jobs, each capped at 10,000 updates,
followed by a diagnostic job and a summary job. Every training/diagnostic job uses
the known serial/10 GB/four-hour request; there is no unverified long queue limit.
Each training allocation has a 210-minute application guard, leaving PBS headroom.
Based on the previous 43-minute/5,000-update run, allow several hours of actual
compute plus queue time; stronger monitoring adds overhead. Runtime is not guaranteed.

Dependencies use `afterany` so failures do not leave the reporting chain blocked.
However, each allocation must verify the previous allocation's successful receipt
and exact state hash before resuming. Missing/mismatched receipts fail promptly;
the chain never silently resumes after an unsafe failure. The diagnostic checks
that the last slot completed safely and candidates are ready.

A productive wall-limited allocation saves state and permits the next slot to
continue. If the five slots are insufficient, the run remains explicitly
`incomplete`, diagnostics do not run, and the summary reports that more runtime
is needed. Signals/zero-progress stops exit 75 rather than consuming the chain.

At completion or a review stop, `final_average.pt` and `final_tail.pt` are diagnosed
with seed 942407: 16,384 uniform points, 5,120 structured points, three searches
of 128 starts/100 steps. Each diagnostic has a 90-minute guard. Both reports must
match their checkpoint hashes. As before, deterministic structured grids can
repeat; uniform/search seeds are new. No report claims integrability certification.

## Publish locally

```bash
cd /Users/ed/Documents/Projects/DG/Nijenhuis/github
git add extra/cs_long.py extra/cs_followup_inputs.py extra/test_cs_long.py configs/ambient/s4xs1xs1_cs/long_replay.yaml configs/ambient/s4xs1xs1_cs/diagnostic_long_replay.yaml docs/s4xt2_cs_long.md
git commit -m "Add monitored long CS continuation with chained restarts"
git push origin s4t2
```

After the push succeeds, copy the two new ignored scripts:

```bash
scp scratch/launch_S4xS1xS1_CS_long.sh scratch/submit_S4xS1xS1_CS_long.sh ehirst@lovelace:/home/lovelace/proj/proj925/ehirst/nijenhuis/
```

## Launch on Lovelace

```bash
cd /home/lovelace/proj/proj925/ehirst/nijenhuis/github
git switch s4t2 && git pull --ff-only origin s4t2
```

After that succeeds:

```bash
cd /home/lovelace/proj/proj925/ehirst/nijenhuis
bash launch_S4xS1xS1_CS_long.sh
qstat -u "$USER"
```

Use **bash**, not qsub, for the launcher. It creates a unique campaign under
`outputs/s4xs1xs1_cs/long_replay_s407_XXXXXXXX/`, records seven job IDs in
`submission.tsv`, and prints the exact command to read `run/comparison.txt`.
The summary includes both candidates and diagnostics without choosing a winner.

Do not rerun the launcher after a partial submission failure: inspect `qstat`
and `submission.tsv`. Old results and source files are not overwritten. Never
fabricate a success receipt. If another allocation is required, submit the
training PBS script with the same environment and the next allocation number;
update ALLOCATIONS for diagnosis accordingly. After a failed allocation, inspect
its failure before resubmitting that allocation number. If only a diagnostic
fails, retain both candidate checkpoints and rerun that diagnostic separately
into a new empty directory; do not restart the entire training campaign.

## Verification scope

Local tests use small synthetic ACS/CS fixtures but actual Nijenhuis derivatives
and Adam updates. They exercise dual selection, constant LR, monitor expansion,
interruption during a pending observation, exact resume including optimizer and
all RNG streams, plateau/growth gates, hash-bound allocation receipts, no-op
terminal slots, source preflight, summary hashes and mocked PBS job chains.
They do not reproduce the remote learned field or establish convergence to a CS.
