#!/bin/bash

# build PHIST with the builtin (fortran) kernels,
# enable repartitioning by ParMetis
# enable graph coloring for using the MC_CARP-CG method

unset TRILINOS_HOME

CC="mpicc -mt_mpi"
CXX="mpicxx -mt_mpi"
FC="mpif90 -mt_mpi"
LIKWID_HOME=/apps/likwid/3.1.1/

cmake \
      -DPHIST_KERNEL_LIB=builtin \
      -DCMAKE_BUILD_TYPE=Release \
      -DCOLPACK_DIR=/home/hpc/essex/hpc023//local/ \
      -DPARMETIS_DIR=/home/hpc/essex/hpc023/local/ \
                    /home/hpc/essex/hpc023/essex/phist/

#      -DLIKWID_PERFMON=ON \
#      -DCMAKE_INSTALL_PREFIX=~/work/essex/installs/ \
