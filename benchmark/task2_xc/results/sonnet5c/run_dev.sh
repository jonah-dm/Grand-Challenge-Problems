#!/bin/bash
# usage: run_dev.sh TAG BASIS  (env S5_* passed through)
cd /home/user/Grand-Challenge-Problems/benchmark/task2_xc
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NTHREADS=1 python score.py dev_set kit/dev_reference.json file:results/sonnet5c/submission.py --basis $2 --out results/sonnet5c/runs/$1.json > results/sonnet5c/runs/$1.log 2>&1
