#!/bin/bash

# benchmark 'tall skinny' kernels on a system with 2 CPU sockets and 2 GPUs with GHOST

set -o verbose
 mpirun_rrze -np 4 ./phist_bench 2> /dev/null
 mpirun_rrze -np 1 env GHOST_TYPE=cuda ./Dmvec_benchmarks 10000000 40 4 20 2> /dev/null |grep BENCH
 mpirun_rrze -np 1 env GHOST_TYPE=work ./Dmvec_benchmarks 10000000 40 4 20 2> /dev/null |grep BENCH
 mpirun_rrze -np 2 env GHOST_TYPE=work ./Dmvec_benchmarks 10000000 40 4 20 2> /dev/null |grep BENCH
 mpirun_rrze -np 4 ./Dmvec_benchmarks 10000000 40 4 20 2> /dev/null |grep BENCH

