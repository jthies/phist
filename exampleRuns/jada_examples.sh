#!/bin/bash

# compute the largest few eigenvalues of a 3D laplacian
mpirun -np 2 ./phist_Dsubspacejada BENCH3D-32-A0 I opts-herm-LM.txt
