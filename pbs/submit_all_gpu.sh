#!/bin/bash
# Submit all 7 GPU sweep jobs at once (phase 1 baselines + phase 2 hypothesis tests).
# Each job runs on miggpu24h (24h walltime, half-A100, 20 GB VRAM, 30 GB RAM).

set -e
cd "$(dirname "$0")/.."

echo "Submitting 7 GPU jobs to miggpu24h..."
echo
for pbs in \
    pbs/train_gpu_B1_s3s3_cs.pbs \
    pbs/train_gpu_B2_s5s1_cs.pbs \
    pbs/train_gpu_B3_s4s2_cs.pbs \
    pbs/train_gpu_B4_s4s2_acs.pbs \
    pbs/train_gpu_H1_s3s3_light.pbs \
    pbs/train_gpu_H2_s4s2_small.pbs \
    pbs/train_gpu_H3_s3s3_lbfgs.pbs
do
    jid=$(qsub "$pbs")
    printf "  %-40s %s\n" "$pbs" "$jid"
done
echo
echo "Monitor with: qstat -u $USER"
echo "Or: tail -f logs/gpu_*.log"
