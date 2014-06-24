#!/bin/bash
# This is a script that should be used for the development of the 
# the FEAST/CARP-CG prototype. In this first step, we provide all
# required core operations (CARP etc.) in the fortran kernel implementation.

# use Fortran kernels for experimental CARP implementation
PHIST_KERNEL_LIB=fortran
# compile without Trilinios
unset TRILINOS_HOME
unset Trilinos_HOME

# use ParMETIS for global (re-)partitioning
PARMETIS_DIR=/tools/modulesystem/tools/parmetis/parmetis-4.0.3/install/sled11.x86_64.gcc-4.8.2.release.int64
# use ColPack for node-local graph coloring
COLPACK_DIR=/unsecured/thie_jo/essex/installs
# use MPI compilers
CC=mpicc
CXX=mpicxx
FC=mpif90
SRC_DIR=${HOME}/essex/phist
cmake -DCMAKE_BUILD_TYPE=Release \
      -DPARMETIS_DIR=${PARMETIS_DIR} \
      -DCOLPACK_DIR=${COLPACK_DIR} \
      ${SRC_DIR}
