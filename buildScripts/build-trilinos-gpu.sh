#!/bin/bash
module purge
source <(spack module loads --dependencies openmpi+cuda %gcc@5.4.0)
spack load cmake
spack load intel-mkl@2018.1.163
spack load doxygen

module list

BLASLIB=$MKLROOT/lib/intel64/libmkl_intel_ilp64.so
PREFIX=/unsecured/thie_jo/local

export CC=mpicc
export FC=mpif90
export CXX=mpicxx
export OMPI_CC=gcc
export OMPI_FC=gfortran
#export OMPI_CXX=/home/thie_jo/trilinos/packages/kokkos/bin/nvcc_wrapper
export OMPI_CXX=/home/thie_jo/bin/nvcc_wrapper
export NVCC_WRAPPER_DEFAULT_COMPILER=`which g++`

cmake -DCMAKE_INSTALL_PREFIX:PATH=${PREFIX} \
      -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_SHARED_LIBS:BOOL=ON \
      -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=OFF \
      -DTrilinos_ENABLE_Epetra:BOOL=ON \
      -DTrilinos_ENABLE_EpetraExt:BOOL=ON \
      -DTrilinos_ENABLE_Ifpack:BOOL=ON \
      -DTrilinos_ENABLE_Ifpack2:BOOL=ON \
      -DTrilinos_ENABLE_Amesos:BOOL=ON \
      -DTrilinos_ENABLE_Amesos2:BOOL=ON \
      -DTrilinos_ENABLE_AztecOO:BOOL=ON \
      -DTrilinos_ENABLE_Anasazi:BOOL=ON \
      -DTrilinos_ENABLE_Belos:BOOL=ON \
      -DTrilinos_ENABLE_Galeri:BOOL=ON \
      -DTrilinos_ENABLE_Tpetra:BOOL=ON \
      -DTrilinos_ENABLE_Kokkos:BOOL=ON \
      -DTrilinos_ENABLE_NOX:BOOL=ON \
      -DTrilinos_ENABLE_ML:BOOL=ON \
      -DTrilinos_ENABLE_Zoltan:BOOL=ON \
      -DTrilinos_ENABLE_Isorropia:BOOL=ON \
      -DTPL_ENABLE_MPI:BOOL=ON \
      -DTrilinos_ENABLE_OpenMP:BOOL=ON \
      -DTPL_ENABLE_MKL:BOOL=OFF \
      -DTPL_ENABLE_CUDA:BOOL=ON \
      -DTPL_ENABLE_BLAS:BOOL=ON \
      -DTPL_BLAS_LIBRARIES="${BLASLIB}" \
      -DTPL_ENABLE_LAPACK:BOOL=ON \
      -DTPL_LAPACK_LIBRARIES="${BLASLIB}" \
      \
      -DTrilinos_ENABLE_FLOAT=ON \
      -DTrilinos_ENABLE_COMPLEX=ON \
      -DTrilinos_ENABLE_COMPLEX_FLOAT=ON \
      -DTrilinos_ENABLE_COMPLEX_DOUBLE=ON \
      -DKokkosKernels_INST_FLOAT=OFF \
      -DKokkosKernels_INST_DOUBLE=ON \
      -DKokkosKernels_INST_COMPLEX_FLOAT=OFF \
      -DKokkosKernels_INST_COMPLEX_DOUBLE=ON \
      -DTpetra_INST_FLOAT=OFF \
      -DTpetra_INST_DOUBLE=ON \
      -DTpetra_INST_COMPLEX_FLOAT=OFF \
      -DTpetra_INST_COMPLEX_DOUBLE=ON \
      \
      -DKokkos_ENABLE_Cuda_Lambda:BOOL=ON \
      -DTrilinos_CXX11_FLAGS="-std=c++11 --expt-extended-lambda" \
      -DKokkos_ENABLE_Cuda_UVM:BOOL=ON \
      -DTpetra_INST_CUDA:BOOL=ON \
      -DTpetra_INST_SERIAL:BOOL=ON \
      -DTpetra_INST_OPENMP:BOOL=ON \
      -DTpetra_INST_THREADS:BOOL=ON \
      -DNVCC_WRAPPER=/home/thie_jo/bin/nvcc_wrapper \
      -DCUDA_NVCC_FLAGS="-gencode arch=compute_70,code=sm_70" \
      .. || exit $LINENO

make -j || make -j 16 || make -j 1 || exit $LINENO
make install || exit $LINENO
