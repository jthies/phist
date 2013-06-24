0) Introduction

the DLR contribution to the ESSEX (Equipping Sparse Solvers for Exa-Scale)
has the working title "PHIST" (Pipelined Hybrid-parallel Iterative Solver Toolkit).
The svn repository can be checked out using the command (on DLR systems):

  svn co https://svn-dmz.sistec.dlr.de/svn/dfg-essex/trunk dfg-essex

1) Installation

make a build directory:

  mkdir build

you need cmake, on the DLR systems, use

  module add application/cmake-2.8.8

you need to select a "kernel library" that provides the basic operations.
We currently support two kernel libs from the Trilinos project, 
epetra (only real double precision, MPI only) and tpetra (single, double, real, complex,
MPI and node-level parallelism). Soon we will also support ghost. On the DLR systems:

  module add trilinos/trilinos-11.2.4
  export PHIST_KERNEL_LIB=epetra

Now go to the build directory, configure and compile:

  cd build/
  cmake ..
  make
  
To run the few serial unit tests that we already have, type
  
  make test
  
or 

  make check

(to get verbose error messages)

The code coverage of the tests can be shown using
  
  make coverage

3) An example of using kernel functions can be found in

  examples/kernels/main-kernel-usage.c

