##############################################################################
# Copyright (c) 2013-2017, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# This file is part of Spack.
# Created by Todd Gamblin, tgamblin@llnl.gov, All rights reserved.
# LLNL-CODE-647188
#
# For details, see https://github.com/spack/spack
# Please also see the NOTICE and LICENSE files for our notice and the LGPL.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License (as
# published by the Free Software Foundation) version 2.1, February 1999.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
# conditions of the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
##############################################################################

from spack import *


class Phist(CMakePackage):
    """ PHIST - a Pipelined, Hybrid-parallel Iterative Solver Toolkit.
        
PHIST provides implementations of and interfaces to block iterative solvers for sparse linear and eigenvalue problems.
In contrast to other libraries we support multiple backends (e.g. Trilinos, PETSc and our own optimized kernels),
and interfaces in multiple languages such as C, C++, Fortran 2003 and Python. PHIST has a clear focus on
portability and hardware performance: in particular we support row-major storage of block vectors and using GPUs (via
GHOST or Trilinos/Tpetra).
"""

    homepage = "https://bitbucket.org/essex/phist/"

    version('develop',
            git='https://bitbucket.org/essex/phist/phist.git', branch='devel')
    version('master',
            git='https://bitbucket.org/essex/phist/phist.git', branch='master')
    version('1.4.1', '53ca2c2c000a36790e1626ed1eb642e3',
            url='https://bitbucket.org/essex/phist/get/phist-1.4.1.tar.gz')

    def cmake_args(self):
        spec=self.spec
        outlev=2
        if '+quiet' in spec:
                outlev=1
        elif '+verbose' in spec:
                outlev=3
        else:
                outlev=2
        lapack_libs = spec['lapack'].libs + spec['blas'].libs
        kernel_lib=spec.variants['kernel_lib'].value
        if kernel_lib=='petsc+complex':
                kernel_lib='petsc'
                    
        args = ['-DPHIST_KERNEL_LIB=%s' % kernel_lib,
                '-DPHIST_OUTLEV=%d' % outlev,
                '-DTPL_LAPACK_LIBRARIES=%s' % lapack_libs,
                '-DPHIST_ENABLE_MPI:BOOL=%s' % ('ON' if '+mpi' in spec else 'OFF'),
                '-DBUILD_SHARED_LIBS:BOOL=%s' % ('ON' if '+shared' in spec else 'OFF'),
                ];
                        
        return args

    variant('shared',  default=True, description='Enables the build of shared libraries')
    variant('mpi', default=True, description='enable/disable MPI (note that the kernel library must also support this)')
    variant('verbose',default=False,description='be more verbose')
    variant('quiet',default=False,description='be less verbose')
    variant(name='kernel_lib',   default='builtin',
            description='select the kernel library (backend) for phist',
            values=['builtin','epetra','tpetra','petsc','petsc+complex','eigen','ghost'])

    # ###################### Dependencies ##########################

    # Everything should be compiled position independent (-fpic)
    depends_on('cmake@3.5:')
    depends_on('blas')
    depends_on('lapack')
    depends_on('mpi', when='+mpi')
    depends_on('trilinos@12:', when='kernel_lib=tpetra')
    # Epetra backend also works with older Trilinos versions
    depends_on('trilinos@11:12', when='kernel_lib=epetra')
    depends_on('petsc~complex', when='kernel_lib=petsc')
    depends_on('petsc+complex', when='kernel_lib=petsc+complex')
    depends_on('eigen', when='kernel_lib=eigen')


