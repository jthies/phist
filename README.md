PHIST - a Pipelined Hybrid Parallel Iterative Solver Toolkit
============================================================

PHIST provides implementations of and interfaces to block iterative solvers for sparse linear and eigenvalue problems.
In contrast to other libraries we support multiple backends (e.g. Trilinos, PETSc and our own optimized kernels),
and interfaces in multiple languages such as C, C++, Fortran 2003 and Python. PHIST has a clear focus on 
portability and hardware performance: in particular we support row-major storage of block vectors and using GPUs (via 
GHOST or Trilinos/Tpetra).

Disclaimer
==========

This project is still under development, so **please do not expect anything to work!** Just 
kidding... But some things may, and probably will, be broken.
Nevertheless, please report any bugs, issues and feature requests to the [issue 
tracker](https://bitbucket.org/essex/phist/issues).


What is PHIST?
==============

The 'Pipelined Hybrid-parallel Iterative Solver Toolkit' was developed as a
sparse iterative solver framework within the [ESSEX project](http://blogs.fau.de/essex/) (Equipping Sparse Solvers for 
Exascale). The aim of PHIST is to provide an environment for developing iterative solvers for sparse linear 
systems and eigenvalue problems that can tackle the challenges of today's increasingly complex CPUs and
arithmetic coprocessors.

**Pipelined** indicates in the broadest sense that algorithms are optimized for processors with wide 'vector 
units' (SIMD/SIMT on GPUs). Standard schemes may not expose sufficient parallelism to allow performing the same 
operation on 4, 8 or 32 elements independently. Hence PHIST may for instance solve multiple systems with 
different shifts and right-hand sides (but the same matrix) simultaneously,
replacing those that converge. Likewise, pipelining of operations during block orthogonalization allows using faster 
fused kernels. Another 'pipelining' idea that will be exploited in future versions of PHIST is the reformulation of 
schemes like CG and GMRES for allowing overlapping communication and computation (see papers on pipelined GMRES by 
Ghysels et al).

**Hybrid parallel** means that we assume an 'MPI+X' programming model, where only MPI communication between processes is 
assumed and 'X' may be any additional accelerator, CPU or core level programming scheme. The 'X' depends on the 
underlying kernel library used, for instance, [ ghost ](https://bitbucket.org/essex/ghost) uses OpenMP, SIMD intrinsics 
and CUDA, whereas Epetra implements no additional parallelism beyond MPI 
(see [ this wiki page ](https://bitbucket.org/essex/phist/wiki/User_Guide/030 Kernel Libraries.wiki#!supported-kernel-libraries) 
for details on the supported kernel libraries).

More high-level information can be found in the [ wiki ](https://bitbucket.org/essex/phist/wiki/Home.wiki). The main source of information on implementation details should be the headers and source code themselves,
along with the examples found in phist/drivers/ and phist/examples/. It is possible to generate HTML documentation using Doxygen (just type 'make doc' in your build directory, and the output is written to
phist/doc/html). However, this documentation is not always complete and well-formatted or structured.


Trying it out
=============

The git repository can be checked out using the command

  git clone git@bitbucket.org:essex/phist

PHIST uses Cmake, and for trying it out you can just stick with the default settings.
Detailed instructions for compiling PHIST can be found [ here ](https://bitbucket.org/essex/phist/wiki/User_Guide/000 Installation.wiki)

There are two categories of programs being built: drivers are examples of high-level algorithms that can be used to e.g. 
compute some eigenvalues of a matrix. They are installed along with the headers and libraries of PHIST. The other 
category (examples) includes benchmarks and examples of kernels and core functionality which are more interesting to 
advanced users and developers than to end users.

A number of scripts that show how to use these drivers can be found in the exampleRuns/ directory. They may have to be 
adjusted to your system.

Finally there are small examples of how to use PHIST in an external application in the exampleProjects/ subdirectory.

Citing PHIST
============

In order to refer to this software, please cite at least this paper and specify the version of the software used.

- Thies, J. et al.: PHIST: a Pipelined, Hybrid-parallel Iterative Solver 
Toolkit. ACM Trans. Math. Software 46 (4), 2020 [author's version](https://elib.dlr.de/123323/1/phist2018.pdf)

A more extensive list of publications can be found on
[ this wiki page](https://bitbucket.org/essex/phist/wiki/WikiHome.wiki#!Related Publications).


Contact
=======

For questions or comments, please contact j.thies@tudelft.nl
