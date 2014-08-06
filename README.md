0) What is PHIST?
=================

The DLR contribution to the ESSEX (Equipping Sparse Solvers for Exa-Scale)
has the working title "PHIST" (Pipelined Hybrid-parallel Iterative Solver Toolkit).
The git repository can be checked out using the command

  git clone git@bitbucket.org:essex/phist

PHIST contains 

* an abstract `kernel' interface layer defining the basic operations commonly used
  in iterative linear algebra solvers
  
* implementations of the interface using
  - own optimized Fortran 90 kernels (MPI+OpenMP)
  - Trilinos (epetra or tpetra)
  - GHOST (developed in ESSEX, MPI+X)

* some algorithms implemented using the interface layer:
  - Jacobi-Davidson eigenvalue solvers for nonsymmetric matrices,
    suitable for finding a few exterior eigenvalues 
    (standard eigenvalue problems up to now)
    + a simple single-vector JDQR with GMRES preconditioner
    + reduced communication block JDQR with (pipelined) GMRES or MinRes preconditioning
  
  - MinRes and GMRES linear solvers
  - Hybrid parallel CARP-CG (Gordon & Gordon, 2010), a row projection method suitable as inner linear solver in FEAST (interior eigenvalue computations)

Depending on the kernel library, real and complex, single and double precision versions are 
available. The interface functions and all algorithms provided have a simple C interface
callable from basically any high level programming language. Data structures are only passed 
around as opaque void pointers to increase portability. 

The kernel library to be used has to be chosen before configuring/building the project by 
setting the environment variable PHIST_KERNEL_LIB to either

* fortran
* epetra
* tpetra
* ghost

Stand-alone version
===================

The PHIST project can be compiled and used without any additional dependencies if the
`fortran' kernel lib is used. For better performance, one should use the optional 
third-party libraries (TPLs) ParMETIS (for repartitioning the matrix) and ColPack (for 
enabling intra-node parallelism in CARP-CG), or switch to ghost+TPLs. Note that we currently 
do not support the use of GPUs/Xeon PHI even with ghost as kernel lib.

Dependencies and optional packages
==================================

* fortran
** optional:
*** ParMETIS: http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview
*** ColPack: http://cscapes.cs.purdue.edu/coloringpage/software.htm
*** essex/physics: https://bitbucket.org/essex/physics
* ghost
** required:
*** essex/ghost: https://bitbucket.org/essex/ghost
*** essex/physics: https://bitbucket.org/essex/physics
** optional: depends on ghost installation
* epetra/tpetra
** required:
*** Trilinos: http://trilinos.sandia.gov/
** optional: depends on Trilinos installation

PHIST Installation
==================

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
  cmake .. -DTrilinos_HOME=$TRILINOS_HOME -DGHOST_DIR=$PREFIX/lib/ghost -DESSEX-PHYSICS_DIR=$PREFIX/lib/essex-physics -DMPIEXEC=mpirun_rrze
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

