PHIST - a Pipelined Hybrid Parallel Iterative Solver Toolkit
============================================================

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
fused kernels. Another 'pipelining' idea that will be exploited in future versions of PHIST is the reformuation of 
schemes like CG and GMRES for allowing overlapping communication and computation (see papers on pipelined GMRES by 
Ghysels et al).

**Hybrid parallel** means that we assume an 'MPI+X' programming model, where only MPI communication between processes is 
assumed and 'X' may be any additional accelerator, CPU or core level programming scheme. The 'X' depends on the 
underlying kernel library used, for instance, [ ghost ](https://bitbucket.org/essex/ghost) uses OpenMP, SIMD intrinsics 
and CUDA, whereas Epetra implements no additional parallelism beyond MPI 
(see [ this wiki page ](https://bitbucket.org/essex/phist/wiki/User_Guide/030 Kernel Libraries.wiki#!supported-kernel-libraries) 
for details on the supported kernel libraries).

More information can be found in the [ wiki ](https://bitbucket.org/essex/phist/wiki/Home.wiki)


Trying it out
=============

The git repository can be checked out using the command

  git clone git@bitbucket.org:essex/phist

Instructions for compiling PHIST can be found [ here ](https://bitbucket.org/essex/phist/wiki/User_Guide/000 Installation.wiki)

Contact
=======

For questions or comments, please contact Jonas.Thies@DLR.de
