#!/bin/bash

# I included this example script for configuring phist with Tpetra as a kernel library, where Tpetra supports CUDA,
# because this combination is rather tricky to build. I use Trilinos 12.13.0 (the present develop branch of Trilinos)
#
# - I disable CCache, it causes trouble with CUDA in general.
#
# - I explicitly disable support for float and complex data types because I couldn't install Trilinos with both CUDA and 
# those types enabled. In principle phist can cope with this situation, but in my tpetra/cuda experiments I got strange 
# compiler errors concerning 'lost' type casts when building the tests.


# to configure and build the libs, use nvcc_wrapper provided by Kokkos as C++ compiler,
# and add the compiler flag --expt-extended-lambdas, which Kokkos needs
export CXX=mpicxx
export CC=mpicc
export OMPI_CXX=nvcc_wrapper
export OMPI_CC=gcc
export OMPI_MPICXX_CXXLFAGS="-std=c++11 --expt-extended-lambda -fopenmp"

cmake -DPHIST_KERNEL_LIB=tpetra -DCMAKE_BUILD_TYPE=Release -DPHIST_USE_CCACHE=OFF .. || exit 1
#-DPHIST_ENABLE_SP=OFF -DPHIST_ENABLE_COMPLEX=OFF .. || exit 1

make libs -j || make libs || make libs -j 1 || exit 2

# now restore the default OpenMPI behavior
unset OMPI_CXX
unset OMPI_CC
unset OMPI_MPICXX_CXXFLAGS

# and build the drivers and tests
make -j || make || make -j 1 || exit 3




