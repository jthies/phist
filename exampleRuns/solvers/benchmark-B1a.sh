#!/bin/bash

# compute the right-most eigenpairs of a nonsymmetric PDE problem using Anasazi block Krylov-Schur

DRIVER=./phist_Danasazi_krylov_schur
ARGS="BENCH3D-128-B1 I 0 10 LR 1e-6 300 2 40"

MPIEXEC=orterun
MPIFLAGS="-np 1 --bind-to none"

# with GHOST
#ENV="env OMP_SCHEDULE=static"

#with Tpetra (CPU)
#ENV="env OMP_PROC_BIND=false"

#with Tpetra+CUDA
ENV="env OMP_PROC_BIND=false env CUDA_MANAGED_FORCE_DEVICE_ALLOC=1 env CUDA_LAUNCH_BLOCKING=1"


# run the example
${MPIEXEC} ${MPIFLAGS} ${ENV} ${DrIVER} ${ARGS}
