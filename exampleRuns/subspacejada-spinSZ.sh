#!/bin/bash -x

# the phist driver 'subspacejada' can be run either with a parameter file or with a 
# long range of input parameters on the command line, this example shows the second
# option, omputing some of the lowest eigenmodes of the Heisenberg spin chain model.

NP=2
DRIVER=phist_Dsubspacejada
MATFILE=spinSZ20
MASSMAT="I"
# symmatric problem
SYM=1
# we look for 20 eigenpairs
NEIG=20
BLOCKSIZE=4
# at the right-most side of the spectrum
WHICH="SR"
# require approx. 100*eps in accuracy
TOL=1e-8
# allow at most 250 iterations
MAXIT=250
INNER_ITERS=10
# iterate between 10 and 25 vectors
MINBAS=28
MAXBAS=60

mpirun -np ${NP} ./${DRIVER} ${MATFILE} ${MASSMAT} ${SYM} ${NEIG} ${WHICH} \
       ${TOL} ${MAXIT} ${BLOCKSIZE} ${MINBAS} ${MAXBAS} ${BLOCKSIZE} ${INNER_ITERS} 0 0 1 1
