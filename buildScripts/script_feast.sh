#!/bin/bash
# This is a script that should be used for the development of the 
# the FEAST/CARP-CG prototype. In this first step, we provide all
# required core operations (CARP etc.) in the fortran kernel implementation.

# compile without Trilinios and override env var for kernel lib
unset PHIST_KERNEL_LIB
unset TRILINOS_HOME
unset Trilinos_HOME

# use ParMETIS for global (re-)partitioning
PARMETIS_DIR=/tools/modulesystem/tools/parmetis/parmetis-4.0.3/install/sled11.x86_64.gcc-4.3.4.release.int64
# use ColPack for node-local graph coloring
COLPACK_DIR=/tools/modulesystem/tools/ColPack/ColPack-1.0.8/install/sled11.x86_64.gcc-4.3.4.release
# use MPI compilers
CC=mpicc
CXX=mpicxx
FC=mpif90
SRC_DIR=${HOME}/essex/phist
cmake -DCMAKE_BUILD_TYPE=Release \
      -DPHIST_KERNEL_LIB=fortran \
      -DPARMETIS_DIR=${PARMETIS_DIR} \
      -DCOLPACK_DIR=${COLPACK_DIR} \
      -DCMAKE_C_COMPILER=mpicc \
      -DCMAKE_CXX_COMPILER=mpicxx \
      -DCMAKE_Fortran_COMPILER=mpif90 \
	${SRC_DIR}
