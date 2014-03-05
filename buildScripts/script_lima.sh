#!/bin/bash
export CC=/apps/intel/mpi/4.1.1.036/rrze-bin-intel/mpicc
export FC=/apps/intel/mpi/4.1.1.036/rrze-bin-intel/mpif90
export CXX=/apps/intel/mpi/4.1.1.036/rrze-bin-intel/mpicxx
export LIKWID_HOME=/apps/likwid/likwid-3.1beta
cmake -DCMAKE_CXX_FLAGS="-DMPICH_SKIP_MPICXX -openmp -Wno-unused-variable -mkl=sequential" \
      -DPHIST_ENABLE_MPI=on \
      -DCMAKE_C_FLAGS="-openmp -Wno-unused-variable" \
      -DCMAKE_C_COMPILER="/apps/intel/mpi/4.1.1.036/rrze-bin-intel/mpicc" \
      -DCMAKE_CXX_COMPILER="/apps/intel/mpi/4.1.1.036/rrze-bin-intel/mpicxx" \
      -DLAPACK_LIBS="mkl_intel_lp64" \
      -DCMAKE_BUILD_TYPE=Debug \
      -DLIKWID_PERFMON=Off \
	..
