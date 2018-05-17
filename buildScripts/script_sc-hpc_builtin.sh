#!/bin/bash
PHIST_BASE=${HOME}/essex/phist/

module purge
CMP=gcc@7.3.0

spack load ${CMP}
spack  load cmake %${CMP}
#spack load intel-mkl@2018.1.163

source <(spack module loads --dependencies openmpi +thread_multiple +pmi schedulers=slurm %${CMP})
source <(spack module loads --dependencies trilinos %${CMP})


cmake -DPHIST_KERNEL_LIB=builtin \
      -DBLA_VENDOR="Intel10_64lp_seq" \
      -DPHIST_PERFCHECK=ON \
      -DCMAKE_BUILD_TYPE=Release \
      -DPHIST_BENCH_LARGE_N=1000000000 \
      ${PHIST_BASE}
