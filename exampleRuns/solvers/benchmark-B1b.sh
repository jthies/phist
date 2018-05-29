#!/bin/bash

# compute the left-most eigenpairs of a nonsymmetric PDE problem using PHIST subspacejada (aka bjdqr)

DRIVER=./phist_Dsubspacejada
ARGS="BENCH3D-128-B1 I jadaOpts-B1b.txt"

# in this file you should make settings for the kernel library and machine (MPIEXEC, MPIFLAGS and ENV)
source benchmark-env.sh

# run the example
${MPIEXEC} ${MPIFLAGS} ${ENV} ${DRIVER} ${ARGS}
