#!/bin/bash
unset TRILINOS_HOME
module unload intel64
module add intel64/15.0up01
cmake -DPHIST_KERNEL_LIB=builtin \
      -DPHIST_KERNEL_LIB_BUILTIN_AVX2_PREC_KERNELS:BOOL=ON \
      -DMPACK_DIR=${HOME}/local \
      ..
