#!/bin/bash

# compute the smallest few eigenvalues of a 3D laplacian
mpirun -np 2 ./phist_Dsubspacejada BENCH3D-32-A0 I opts-herm-SM.txt

# same problem with IFPACK preconditioner (will only work with Epetra or Tpetra).
# ifpack_options.xml can be found in exampleRuns/precon/
mpirun -np 2 ./phist_Dsubspacejada BENCH3D-32-A0 I opts-herm-SM-ifpack.txt

# solve a generalized EVP
mpirun -np 2 ./phist_Dsubspacejada elast_A.mm elast_B.mm opts-elast.txt
