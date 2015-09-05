#!/usr/bin/env bash


# accuracy tests
for acc in 0 1 2 3 5 7; do
  for nb in 1 2 4 8; do
    for dev1 in 1. 1.e-3 1.e-6 1.e-9 1.e-12; do
      for dev2 in 1. 1.e-3 1.e-6 1.e-9 1.e-12; do
        echo "Run orthog_fusedRR_nb${nb}_acc${acc}_dev_${dev1}_${dev2}.log"
        OMP_SCHEDULE=guided,1000 PHIST_NUM_THREADS=12 mpirun -np 2 ./Dorthog_fusedRR_times 1000 $nb 20 10 $acc $dev1 $dev2 &> orthog_fusedRR_nb${nb}_acc${acc}_dev_${dev1}_${dev2}.log
      done
    done
  done
done
