cd /home/user/Grand-Challenge-Problems/benchmark/task2_xc/results/opus6t/work
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
printf '%s\n' "0.25 0.25" "0.5 0.5" "0.625 0.625" "0.75 0.75" "0.875 0.875" "1.0 1.0" "0.75 0.25" "1.0 0.5" | xargs -P 3 -L 1 sh -c 'python lieb1d.py $0 $1 softH >> v1.out 2>> v1.err'
