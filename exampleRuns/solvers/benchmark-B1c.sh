#!/bin/bash

# compute the left-most eigenpairs of a nonsymmetric PDE problem using PHIST subspacejada (aka bjdqr)
# and the AMG preconditioner ML. This example only works with the epetra kernels.

DRIVER=./phist_Dsubspacejada
ARGS="BENCH3D-128-B1 I jadaOpts-B1c.txt"

MPIEXEC=orterun
MPIFLAGS="-np 56 --bind-to core"

# with Epetra
ENV="env OMP_NUM_THREADS=1"



# run the example
${MPIEXEC} ${MPIFLAGS} ${ENV} ${DrIVER} ${ARGS}
