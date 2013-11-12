#!/bin/bash -x
SCALAR=Z
DRIVER=${SCALAR}jdqr
MATFILE=${SCALAR}jadaTestMat.mm
# we look for 10 eigenpairs
NEIG=10
# at the right-most side of the spectrum
WHICH="SR"
# require approx. 100*eps in accuracy
TOL=2.2e-14
# allow at most 250 iterations
MAXIT=250
# iterate between 10 and 25 vectors
MINBAS=10
MAXBAS=25

./${DRIVER} ${MATFILE} ${NEIG} ${WHICH} ${TOL} ${MAXIT} ${MINBAS} ${MAXBAS}
