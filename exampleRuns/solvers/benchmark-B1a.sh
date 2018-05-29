#!/bin/bash

# compute the right-most eigenpairs of a nonsymmetric PDE problem using Anasazi block Krylov-Schur

DRIVER=./phist_Danasazi_krylov_schur
ARGS="BENCH3D-128-B1 I 0 10 LR 1e-6 300 2 40"

# in this file you should make settings for the kernel library and machine (MPIEXEC, MPIFLAGS and ENV)
source benchmark-env.sh

# run the example
${MPIEXEC} ${MPIFLAGS} ${ENV} ${DRIVER} ${ARGS}
