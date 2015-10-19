#!/usr/bin/env bash

set -e
export OMPI_MCA_hwloc_base_binding_policy=none


# accuracy tests
for acc in 0 1 2 3 4 5; do
  for nb in 1 2 4; do
    for dev1 in 1. 1.e-3 1.e-6 1.e-9 1.e-12; do
      for dev2 in 1. 1.e-3 1.e-6 1.e-9 1.e-12; do
        logname="orthog_fusedRR_small_nv20_nb${nb}_acc${acc}_dev_${dev1}_${dev2}.log"
        echo "Run ${logname}"
        OMP_NUM_THREADS=1 PHIST_NUM_THREADS=1 ./Dorthog_fusedRR 1000 $nb 20 10 $acc $dev1 $dev2 1.e-12 &> ${logname}
      done
    done
  done
done

# performance tests
for acc in 2 3; do # 2 3 4 5
  for nb in 1 2 4; do
    for nv in 10 20 40; do
      dev=1.e-6

      omp_threads=48
      logname="orthog_fusedRR_large_nv${nv}_nb${nb}_acc${acc}_dev${dev}_nt${omp_threads}.log"
      echo "Run ${logname}"
      OMP_NUM_THREADS=${omp_threads} PHIST_NUM_THREADS=${omp_threads} ./Dorthog_fusedRR 4000000 $nb $nv 300 $acc $dev $dev 1.e-12 &> ${logname}
    done
  done
done
