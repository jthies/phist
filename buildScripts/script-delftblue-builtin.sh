#!/bin/bash

# configure phist on DelftBlue,
# demonstrates MPIEXEC settings such that
# you can run "make test" on the login node
# and the tests are executed via srun on compute nodes.

module load 2022r2
module load intel-mkl
module load openmpi
module load trilinos

export CC=`which mpicc`
export FC=`which mpif90`
export CXX=`which mpicxx`

cmake -DMPIEXEC=`which srun` \
      -DMPIEXEC_NUMPROC_FLAG="-n" \
      -DMPIEXEC_MAX_NUMPROCS="4" \
      -DMPIEXEC_PREFLAGS="--time=00:05:00;--mem-per-cpu=1GB" \
      \
      -DPHIST_KERNEL_LIB=builtin \
      ..
