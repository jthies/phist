0) Introduction

the DLR contribution to the ESSEX (Equipping Sparse Solvers for Exa-Scale)
has the working title "PHIST" (Pipelined Hybrid-parallel Iterative Solver Toolkit).
The git repository can be checked out using the command

  git clone git@bitbucket.org:essex/phist

1a) Installing extra ESSEX libraries

We support several different implementations of the basic linear algebra operations
in PHIST, which are set by the variable PHIST_KERNEL_LIB (see below). If you choose
ghost as kernel library, a number of example drivers are activated which also require 
the ESSEX Physics libraries to be installed:

  https://bitbucket.org/essex/physics

Follow the instructions in order to build and install the libraries.
Obviously you also need to build and install GHOST on your system:

  https://bitbucket.org/essex/ghost

1b) PHIST Installation

make a build directory:

  mkdir build

you need cmake, on the DLR systems, use

  module add cmake

you need to select a "kernel library" that provides the basic operations.
We currently support two kernel libs from the Trilinos project, 
epetra (only real double precision, MPI only) and tpetra (single, double, real, complex,
MPI and node-level parallelism), and ghost (from the essex project). On the DLR systems:

  module add trilinos
  export PHIST_KERNEL_LIB=tpetra

On the RRZE systems proceed as follows to build PHIST with GHOST (substitute $PREFIX and $TRILINOS_HOME with according paths).

  module load intel64
  module load cmake
  cd build/
  PHIST_KERNEL_LIB=ghost CC=icc CXX=icpc 
  cmake .. -DTrilinos_HOME=$TRILINOS_HOME -DGHOST_DIR=$PREFIX/lib/ghost -DESSEX_INSTALL_DIR=$PREFIX/ -DMPIEXEC=mpirun_rrze
  make 

It is sufficient to set GHOST_HOME and TRILINOS_HOME as environement variables (or Ghost_HOME, Trilinos_HOME)
If you want to use likwid, also set LIKWID_HOME and pass -DLIKWID_PERFMON to cmake.

On LiMa at RRZE you can use the script provided in buildScripts/script_lima.sh to get 
started.

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

