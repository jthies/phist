This directory contains some examples how one can compute eigenvalues with two PHIST drivers,
the Trilinos/Anasazi block Krylov Schur method and the PHIST subspacejada (BJDQR) method.

The environment variables and MPI options give an indication on how to get good performance from PHIST for different 
kernel libraries. E.g., for builtin, ghost and tpetra we try to pin threads from within PHIST, so it is crucial to 
prevent pinning by OpenMPI or OpenMP. The variant using Tpetra on GPUs requires a set of environment variables, which 
are displayed by Tpetra when running the program without setting them.

Note that the PHIST drivers, when run without arguments, print some usage information.
