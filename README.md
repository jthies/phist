------------------
# What is PHIST? #
------------------

The DLR contribution to the ESSEX (Equipping Sparse Solvers for Exa-Scale)
has the working title "PHIST" (Pipelined Hybrid-parallel Iterative Solver Toolkit).
The git repository can be checked out using the command

  git clone git@bitbucket.org:essex/phist

PHIST contains 

* an abstract `kernel' interface layer defining the basic operations commonly used in 
  iterative linear algebra solvers
* implementations of the interface using
    * builtin sample implementation (MPI+OpenMP)
    * GHOST (developed in ESSEX, MPI+X)
    * Trilinos (epetra or tpetra)
* some algorithms implemented using the interface layer:
    * Jacobi-Davidson eigenvalue solvers for nonsymmetric matrices, suitable for finding a 
    few exterior eigenvalues (standard eigenvalue problems up to now)
        * a simple single-vector JDQR with GMRES preconditioner
        * reduced communication block JDQR with (pipelined) GMRES or MinRes preconditioning
    * MinRes and GMRES linear solvers
    * Hybrid parallel CARP-CG (Gordon & Gordon, 2010), a row projection method suitable as 
    inner linear solver in FEAST (interior eigenvalue computations)
* an interface to the Trilinos library Belos for ghost, which allows using some more 
iterative linear solvers.

Depending on the kernel library, real and complex, single and double precision versions are 
available. The interface functions and all algorithms provided have a simple C interface
callable from basically any high level programming language. Data structures are only passed 
around as opaque void pointers to increase portability. 

The kernel library to be used has to be chosen before configuring/building the project by 
either

* passing -DPHIST_KERNEL_LIB=<choice> to cmake.
* using ccmake and modifying the corresponding entry.

Choices supported right now are:

* builtin (supplied with phist, chosen by default)
* ghost
* epetra
* tpetra

-----------------------
# Stand-alone version #
-----------------------

The PHIST project can be compiled and used without any additional dependencies if the
`builtin' kernel lib is used. For better performance, one should use the optional 
third-party libraries (TPLs) ParMETIS (for repartitioning the matrix) and ColPack (for 
enabling intra-node parallelism in CARP-CG), or switch to ghost+TPLs. Note that we currently 
do not support the use of GPUs/Xeon PHI even with ghost as kernel lib.

--------------------------------------
# Dependencies and optional packages #
--------------------------------------

* builtin
    * optional:
        * ParMETIS:       http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview
        * ColPack:        http://cscapes.cs.purdue.edu/coloringpage/software.htm
        * essex/physics:  https://bitbucket.org/essex/physics (if not available, some
          example physics problems like the spin chain and graphene have been copied into
          the phist repo and will be used, but they might be outdated).
* ghost
    * required:
        * essex/ghost:    https://bitbucket.org/essex/ghost
        * essex/physics:  https://bitbucket.org/essex/physics
    * optional: depends on ghost installation
* epetra/tpetra
    * required:
        * Trilinos: http://trilinos.sandia.gov/
  - * optional: depends on Trilinos installation

Typically CMake will automatically find the TPL if you pass the
variable (TPL_NAME)_DIR to the cmake command, for instance:

  cmake   -DPHIST_KERNEL_LIB=ghost \
          -DGHOST_DIR=<path to ghost lib dir> \
          <path to phist dir>

--------------------------------
# PHIST Installation (details) #
--------------------------------

Some sample build scripts are available in phist/buildScripts, they may provide
a good starting point, although they might be somewhat outdated and paths may have to be 
adjusted.

make a build directory:

  mkdir build

you need cmake and MPI, on the DLR systems, use

  module add cmake openmpi

As a start, use the builtin kernel lib without TPLs:

  cd build/
  cmake -DPHIST_KERNEL_LIB=builtin ..

Available targets for building just a certain component (libraries, 
drivers, or tests) are

  make libs
  make drivers
  make tests

Required for installation via

  make install

is only the libs target.

------------------------------
## Testing the installation ##
------------------------------

To run some unit tests type
  
  make test
  
or 

  make check

(to get verbose error messages)

The code coverage of the tests can be shown using
  
  make coverage

----------------------------------------------------
## Installation instructions for the RRZE systems ##
----------------------------------------------------

NOTE: This section is a bit outdated as Trilinos is now no longer
a dependency of phist.

 Proceed as follows to build PHIST with GHOST (substitute 
$PREFIX and $TRILINOS_HOME with according paths).

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

-----------------------
# Directory structure #
-----------------------

* *phist*                main project directory
    * *src*              primary source directory
        * *kernels*      abstract interface definition
            * *<subdir>* subdirectories with concrete implementations
        * *core*         common intermediate level routines like orthogonalization
        * ... (TODO: finish documentation)

-------------------
# Example Drivers #
-------------------

By default, only double precision versions of libraries and drivers are compiled. There 
are four types of drivers

* *test drivers* that follow the naming convention phist-X.Y.Z-<name>-test. These 
executables are run by the 'make test' command, but they can also be run as stand-alone
programs for debugging purposes. Their sources are in phist/test/ and are based on the 
googletest framework.

* *benchmark drivers* like DcrsMat_times_mvec_times which are written to benchmark certain 
kernel operations

* *algorithm drivers* that read a matrix from a file and solve e.g. an eigenvalue problem. 
Examples are phist_Djdqr, phist_Zsubspacejada etc. These drivers take a range of command 
line 
parameters, a list of which is printed when the executable is called without input args.

* *application drivers* that do not read a matrix from a file but build it from a function 
defined in the essex/physics repository. Examples of such drivers are phist_spinChain (which 
uses Djdqr) or phist_graphene_carp (which solves linear systems with graphene matrices and 
complex shifts using CARP-CG).
