#!/bin/bash
# This is a script that should be used for the development of the 
# the FEAST/CARP-CG prototype. In this first step, we provide all
# required core operations (CARP etc.) in the fortran kernel implementation.

# use Fortran kernels for experimental CARP implementation
PHIST_KERNEL_LIB=fortran
# compile without Trilinios
unset TRILINOS_HOME
unset Trilinos_HOME
# use MPI compilers
CC=mpicc
CXX=mpicxx
FC=mpif90
SRC_DIR=${HOME}/essex/phist
cmake -DCMAKE_BUILD_TYPE=Debug ${SRC_DIR}
