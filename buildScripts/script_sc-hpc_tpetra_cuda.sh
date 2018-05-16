#!/bin/bash

export PHIST_BASE=${HOME}/essex/phist/
export TRILINOS_HOME=/unsecured/thie_jo/local/

module purge
spack load cmake@3.11.0 %gcc@5.4.0
spack load intel-mkl@2018.1.163

source <(spack module loads --dependencies openmpi+cuda %gcc@5.4.0)

export CMAKE_PREFIX_PATH=${TRILINOS_HOME}:${CMAKE_PREFIX_PATH}
export CC=mpicc
export CXX=mpicxx
export OMPI_CXX=/unsecured/thie_jo/local/bin/nvcc_wrapper
export OMPI_CXXFLAGS="-O2 -fopenmp -expt-extended-lambda -std=c++11"

cmake -DPHIST_KERNEL_LIB=tpetra \
      -DBLA_VENDOR="Intel10_64lp_seq" \
      -DPHIST_PERFCHECK=ON \
      -DCMAKE_BUILD_TYPE=Release \
      -DPHIST_ENABLE_COMPLEX=OFF \
      -DPHIST_ENABLE_SP=OFF \
      -DPHIST_BENCH_LARGE_N=1000000000 \
      ${PHIST_BASE}

make -j 
