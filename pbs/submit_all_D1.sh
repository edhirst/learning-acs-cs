#!/bin/bash
# Submit all 5 Round-1 derivative-program jobs (miggpu24h).
set -e
cd "$(dirname "$0")/.."

echo -n "pbs/train_gpu_D1_s3s3_cs.pbs: "; qsub pbs/train_gpu_D1_s3s3_cs.pbs
echo -n "pbs/train_gpu_D1_s5s1_cs.pbs: "; qsub pbs/train_gpu_D1_s5s1_cs.pbs
echo -n "pbs/train_gpu_D1_s4s2_acs.pbs: "; qsub pbs/train_gpu_D1_s4s2_acs.pbs
# echo -n "pbs/train_gpu_D1_s3s1_cs.pbs: "; qsub pbs/train_gpu_D1_s3s1_cs.pbs
# echo -n "pbs/train_gpu_D1_s2s2_cs.pbs: "; qsub pbs/train_gpu_D1_s2s2_cs.pbs

echo; echo "Monitor with: qstat -u $USER"
