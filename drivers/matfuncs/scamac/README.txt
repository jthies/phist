Disclaimer
----------

This ScaMaC version (0.8) is a development version, and not ready for distribution.
No guarantee is given, no responsibility assumed.
This version has been tested on 4 different systems (3 Linux + 1 MacOs) without problems.
No functionality is really broken and the interface is nearly converged,
but the TODO list remains long, the documentation is a joke, and most matrix examples are still missing.

Availability
------------

All ScaMaC files are found in the ESSEX repository at Bitbucket.
(https://bitbucket.org/essex/matrixcollection)
For now, the repository is not open to the public,
and everybody is asked to *not* distribute the files outside of the ESSEX project.

Dependencies
------------

  General build dependencies:
  - CMake (version >= 3.0 [probably])
  - C compiler
  
  ScaMaC toolkit ("scamact"):
  - BLAS & LAPACK [for Lanczos & spectrum]
  - libpng [for sparsity pattern plots]

  SLEPc minimal working example:
  - SLEPc and, therefore, PETSc.
    Tested versions: PETSc 3.7.6
                     SLEPc 3.7.4
    but other versions should work just fine (requirements are minimal)

Installation
------------

1. Put all files in a directory "FOO", e.g., clone the repository.

2. Create a build directory "GOO". We prefer out-of-source builds.

3. In the directory "GOO", call CMake from the terminal as "cmake FOO".
   The cmake script supports various options, which can be set during the call to CMake,
   or after CMake has finished through "ccmake .", which starts an interactive session.
   Options [with default value ON/OFF] include:
   SCAMAC_USE_64               [ON ]   -- build library with int64 index type
   SCAMAC_USE_BLAS             [OFF]   -- use external BLAS    (required for spectrum/Lanczos)
   SCAMAC_USE_LAPACK           [OFF]   -- use external LAPACK  (required for spectrum/Lanczos)
   SCAMAC_USE_OPENMP           [ON]    -- build with OpenMP    (required for min. work ex.)
   SCAMAC_USE_PNG              [OFF]   -- build with libpng    (required for pattern plots)

4. Build the library, the ScaMaC toolkit ("scamact"), and all minimal working examples with "make".
   Try "make test" afterwards.
   "make install" puts all libraries & executables into the following default locations:
    libraries "libscamac" and "libscamactools" -> GOO/lib
    headers -> GOO/include
    toolkit "scamact" -> GOO/toolbox
    minimal working examples -> GOO/mwe
    minimal working examples that require external libraries -> GOO/mwe_external

5. To build the minimal working examples in GOO/mwe_external:
   - for mwe_slepc:
     set environment variables PETSC_DIR, PETSC_ARCH, SLEPC_DIR for your PETSc/SLEPc installation
     "make mwe_slepc" in GOO/mwe_external" builds the SLEPc example.
     [We here rely on the PETSc/SLEPc makefiles.]
   - for mwe_mpi:
     "make mwe_mpi_count" in GOO/mwe_external" builds the MPI example with mpicc.


Usage
-----

For matrix generation, only libscamac.a is required.
With minimal requirements, compilation is possible in the form of
$(CC) -I$(SCAMAC_INCDIR) -o mwe mwe.c -L$(SCAMAC_LIBDIR) -lscamac -lm
for example
gcc -IGOO/include -o mwe_serial_query mwe_serial_query.c -LGOO/lib -lscamac -lm
will work.

Basic use of the ScaMaC library functions is illustrated through the minimal working examples:

 - mwe_serial_query.c      : Basic example - obtain the dimension of a matrix, without generating any rows.
 - mwe_serial_scale.c      : Basic example - scale up the dimension of a matrix by using different matrix parameters.
                             [This example is expected to terminate with an "error OVERFLOW" message.]
 - mwe_serial_count.c      : Basic example - generate a matrix row by row, and count the non-zeros.
 - mwe_openmp_count.c      : Basic example with OpenMP - generate the matrix rows in parallel, and count the non-zeros.
 - mwe_openmp_statistics.c : Assemble statistics of a ScaMaC matrix.
 - mwe_plot_serial.c       : Plot the sparsity pattern of a ScaMaC matrix.
 - mwe_plot_openmp.c       : ~~~

 - mwe_mpi_count.c         : Basic example with MPI - generate the matrix rows in parallel, and count the non-zeros.
 - mwe_slepc.c             : Integration with SLEPc eigensolvers -- compute an eigenvalue.

With the exception of mwe_serial_scale, all MWEs expect an argument string on the command line, in the form given below, say,
mwe_* Anderson,Lx=100,Ly=50,Lz=20,t=2.0,ranpot=4.0,seed=0x3245,boundary_conditions=periodic
which sets some parameters while the others keep their default values.

The ScaMaC toolkit "scamact" (beware of incorrect hyphenation!) provides several options for exploration of the ScaMaC matrices.
For example, try the following:

$ scamact
  This prints a usage info.

$ scamact list
  This lists the available matrices
  In the output list, you recognize, e.g., the example "Hubbard".

$ scamact info Hubbard
  This prints the description of the example "Hubbard", including a list of the parameters of this example.
  Among the parameters, we recognize, e.g., the parameters "n_sites" and "n_fermions", which take integer values.

$ scamact query Hubbard
  This prints information about the example "Hubbard" to the respective parameter values, including the matrix dimension.
  Every parameter has a default value. To set a parameter, we use the standard syntax of ScaMaC argument strings, namely:

$ scamact query Hubbard,n_sites=12
  Notice how the matrix dimension has increased. Now, try
$ scamact query Hubbard,n_sites=2*n,n_fermions=n
  for n=5,6,7,8,9,10
  
Some parameter values that will not work for different reasons:
$ scamact query Hubbard,n_sites=0
$ scamact query Hubbard,n_sites=10,n_fermions=20
$ scamact query Hubbard,n_sites=30,n_fermions=15
  This should work.
$ scamact query Hubbard,n_sites=32,n_fermions=16
  This should not work.

None of the above calls has yet generated any matrix rows.
The following calls generate the entire matrix, row by row.

$ scamact stat Hubbard
  Prints pattern and value statistics.

$ scamact plot Hubbard
  Plots the sparsity pattern.

None of the above calls has stored the generated matrix.
The following calls have to keep the entire (sparse) matrix in memory.

$ scamact lanczos Hubbard
  Estimate the minimal and maximal eigenvalue [symmetric matrices only]

$ scamact spectrum Hubbard
  Compute the *entire* spectrum.
  Requires us to store the entire matrix as a dense (!) matrix in memory. 
  Only recommended for truly small matrices.

$ scamact output mm Hubbard
  Output matrix in MatrixMarket [mm] format.
