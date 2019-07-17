phist-1.7.5
    - patches for compiling on non-x86 architectures

phist.1.7.4
    - some new kernels and fixes to the Python interface so it can function as a backend PyNCT
    - fixes for xSDK 0.4.0 (dec'18) release
    - new cmake option PHIST_HOST_OPTIMIZE

phist-1.7.3
    - xSDK 0.4.0 release candidate (fall 2018)

phist-1.7.2
    - between-releases version for the xSDK 0.4.0 (fall 2018) release
    - some improved xSDK compatibility like allowing to turn off host-specific optimizations via a CMake flag

phist-1.6.1
    - performance fix for ghost
    - bug fixes to BiCGStab as a correction solver
    - fixes to high precision stuff
    - avoid reallocating temporary vectors in refactored jadaOp

phist-1.6.0
    - imporoved support for preconditioning,
    - a refactored jadaOp,
    - GPU support via ghost and tpetra
    - our own ortho manager for Anasazi and Belos

phist-1.5.2
    - development version (before refactoring of jadaOp)

phist-1.4.3
    - fix compilation with petsc+complex

phist-1.4.2     
    - various bug fixes concerning finding blas/lapack libs (in particular MKL and in combination with the spack package manager)

phist-1.4.1
    - minor bug fixes reg. cmake and hanging ghost tests

phist-1.4.0
    - updated Tpetra interface
    - improved support for generalized EVP
    - more linear solvers
    - a preconditioning interface

phist-1.2.0     
    - bug fixes
    
phist-1.0.0    
    - stable kernel and core interfaces
    - subspacejada as the primary eignesolver for exterior eigenpairs
    - iterative linear solvers GMRES and MINRES for use in Jacobi-Davidson
    - a simple preconditioning interface
    - CARP-CG as an alternative linear solver
    - support for GHOST 1.0.x (including hybrid CPU/GPU)
