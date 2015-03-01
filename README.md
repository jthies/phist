PHIST documentation
===================

Disclaimer
==========

This is a pre-alpha release, so **please do not expect anything to work!** Just 
kidding... But some things may, and probably will, be broken.
Nevertheless, please report any bugs, issues and feature requests to the [issue 
tracker](https://bitbucket.org/essex/phist/issues).


What is PHIST?
==============

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
    * Jacobi-Davidson eigenvalue solvers for symmetric and nonsymmetric matrices, suitable for finding a 
      few exterior eigenvalues (standard eigenvalue problems up to now)
        * a simple single-vector JDQR with GMRES preconditioner
        * reduced communication block JDQR with blocked GMRES or MinRes preconditioning
    * blocked MinRes and GMRES linear solvers
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

* passing `-DPHIST_KERNEL_LIB=<choice>' to cmake.
* using ccmake and modifying the corresponding entry.

Choices supported right now are:

* builtin (supplied with phist, chosen by default)
* ghost
* epetra
* tpetra


Stand-alone version
===================

The PHIST project can be compiled and used without any additional dependencies if the
`builtin' kernel lib is used. For better performance, one should use the optional 
third-party libraries (TPLs) ParMETIS (for repartitioning the matrix) and ColPack (for 
enabling intra-node parallelism in CARP-CG), or switch to ghost+TPLs. Note that we currently 
do not support the use of GPUs/Xeon Phi even with ghost as kernel lib.


Dependencies and optional packages
==================================

* builtin
    * optional:
        * ParMETIS:       http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview
          (configured with IDXTYPEWIDTH=64)
        * ColPack:        http://cscapes.cs.purdue.edu/coloringpage/software.htm
        * essex/physics:  https://bitbucket.org/essex/physics (if not available, some
          example physics problems like the spin chain and graphene have been copied into
          the phist repo and will be used, but they might be outdated).
* ghost
    * required:
        * essex/ghost:    https://bitbucket.org/essex/ghost
          (configured with suitable block sizes!)
    * optional:
        * essex/physics:  https://bitbucket.org/essex/physics
        * additional libraries depend on the ghost installation
* epetra/tpetra
    * required:
        * Trilinos: http://trilinos.sandia.gov/
  - * optional: depends on the Trilinos installation

Typically CMake will automatically find the TPL if you pass the
variable `(TPL_NAME)_DIR' to the cmake command, for instance:

    > cmake   -DPHIST_KERNEL_LIB=ghost \
              -DGHOST_DIR=<path to ghost lib dir> \
              <path to phist dir>


PHIST Installation (details)
============================

Some sample build scripts are available in phist/buildScripts, they may provide
a good starting point, although they might be somewhat outdated and paths may have to be 
adjusted.

make a build directory:

    > mkdir build

you need cmake and MPI, on the DLR systems, use

    > module add cmake openmpi

As a start, use the builtin kernel lib without TPLs:

    > cd build/
    > cmake -DPHIST_KERNEL_LIB=builtin ..

Available targets for building just a certain component (libraries, 
drivers, or tests) are

    > make libs
    > make drivers
    > make tests

Required for installation via

    > make install

are only the targets libs and drivers.

On typical clusters it is recommended to set the C/C++/Fortran compilers to
the MPI wrappers to be used. This is in contrast to the way CMake usually handles MPI
(by trying to extract compiler flags and library info from the wrapper and then using
the compilers directly). We find the method used in PHIST more portable and to the intention
of cluster administrators who provide the wrappers. For example, use

    > cmake [...] \
            -DCMAKE_C_COMPILER="[path/]mpicc" \
            -DCMAKE_CXX_COMPILER="[path/]mpicxx" \
            -DCMAKE_Fortran_COMPILER="[path/]mpif90" \
            [...]


Testing the installation
------------------------

To run some unit tests type
  
    > make test
  
or 

    > make check

(to get verbose error messages)

The code coverage of the tests can be shown using
  
    > make coverage


Installation instructions for the RRZE systems
----------------------------------------------

 Proceed as follows to build PHIST with GHOST (substitute 
$PREFIX with the path where ghost and essex-physics are installed).
This is the easiest way to make PHIST find GHOST. It is also possible
to pass the `GHOST_DIR' variable to cmake, insterad. Note that essex-physics
is an optional package for phist, as is Trilinos.

    > module load intel64
    > module load cmake
    > cd build/
    > export TRILINOS_HOME=<...> # optional, to enable using TSQR and Teuchos timers
    > export CC="mpicc -mt_mpi"               \
             CXX="mpicxx -mt_mpi"             \
             FC="mpif90 -mt_mpi"
    > cmake  -DPHIST_KERNEL_LIB=ghost       \
             -DGHOST_DIR=$PREFIX/lib/ghost  \
             -DMPIEXEC=mpirun_rrze          \
           ..
    > make 

If you want to use likwid, also set `LIKWID_HOME' and pass `-DLIKWID_PERFMON' to cmake.


Directory structure
===================

* *phist*                main project directory
    * *src*              primary source directory
        * *tools*        definitions and helper functions for output, timing,
                         instrumentalization, ...
        * *kernels*      abstract interface definition
            * *<subdir>* subdirectories with concrete implementations
        * *core*         common intermediate level routines like orthogonalization
        * *krylov*       blocked Krylov subspace methods (mostly linear solvers)
        * *jada*         blocked Jacobi-Davidson eigenvalue solvers
        * *feast*        FEAST eigenvalue solver (unfinished)
    * *drivers*          source for executables with a command-line interface
                         to the main algorithms implemented in phist
    * *test*             source & data of our extensive test suite for kernels
                         and algorithms
      * *matrices*       test data
      * *<subdir>*       subdirectories with tests for corresponding directories
                         in *src*
    * *python*           python interface to phist (incomplete, untested)
    * *examples*         source for additional executables that did not make it
                         into the *drivers* directory
    * *cmake*            additional cmake modules, mostly for finding TPLs
    * *doc*              Doxygen documentation (built with make doc if doxygen
                         is available)
    * *exampleRuns*      output & scripts of calculations done by us, probably
                         outdated
    * *lib*              third-party libraries included (currently only googletest)


Example Drivers
===============

By default, only double precision versions of libraries and drivers are compiled. There 
are four types of drivers

* *algorithm drivers* that read a matrix from a file (or use a predefined matrix generator)
  and solve e.g. an eigenvalue problem. 
  Examples are `phist_Djdqr', `phist_Zsubspacejada' etc. These drivers take a range of command 
  line 
  parameters, a list of which is printed when the executable is called without input args.

* *benchmark drivers* like `phist_DcrsMat_times_mvec_times' which are written to benchmark certain 
  kernel operations

* *test drivers* that follow the naming convention phist-X.Y.Z-<name>-test. These 
  executables are run by the `make test' command, but they can also be run as stand-alone
  programs for debugging purposes. Their sources are in phist/test/ and are based on the 
  googletest framework.
