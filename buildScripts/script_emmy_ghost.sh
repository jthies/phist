#!/bin/bash

# build PHIST with kernels provided by GHOST,
# which should be installed in the same directory
# set as "CMAKE_INSTALL_PREFIX" so that it is found.
# 
# Use Trilinos to provide TSQR for the block JaDa (subspacejada)
# method.

export LIKWID_HOME=/apps/likwid/3.1.1
export TRILINOS_HOME=${HOME}/Trilinos/11.12.1

export CC="mpicc -mt_mpi"
export CXX="mpicxx -mt_mpi"
export FC="mpif90 -mt_mpi"
SRC_DIR=${HOME}/essex/phist

cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=/home/hpc/essex/hpc023/work/essex/installs \
      -DPHIST_KERNEL_LIB=ghost \
      -DPHIST_MVECS_ROW_MAJOR=1 \
      -DPHIST_SDMATS_ROW_MAJOR=0 \
      -DPHIST_USE_SELL=ON \
      -DPHIST_SELL_C=1 \
      -DPHIST_SELL_SIGMA=1 \
	${SRC_DIR}
