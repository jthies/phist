'''Python wrapper for phist_kernels library, unit tests
'''

#--------------------------------------------------------------------------------
# imports
from __future__ import print_function
import ctypes as _ct

from phist_tools import *
from phist_kernels import *

try:
    import unittest2 as unittest
except ImportError:
    import unittest


#--------------------------------------------------------------------------------
# tests

class PhistKernelsTest(unittest.TestCase):

    def test_some_kernels(self):

        # initialize
        PYST_CHK_IERR(phist_kernels_init)

        # create comm
        myComm = comm_ptr()
        PYST_CHK_IERR(phist_comm_create, myComm)

        # check comm size and rank
        nprocs = c_int()
        me = c_int()
        PYST_CHK_IERR(phist_comm_get_size, myComm, nprocs)
        PYST_CHK_IERR(phist_comm_get_rank, myComm, me)
        print('nprocs', nprocs, 'me', me)

        # create map
        nglobal = gidx(10)
        myMap = map_ptr()
        PYST_CHK_IERR(phist_map_create, myMap, myComm, nglobal)

        # get number of local elements
        nlocal = lidx()
        PYST_CHK_IERR(phist_map_get_local_length, myMap, nlocal)
        print('global elements', nglobal, 'local elements', nlocal)

        # delete map
        PYST_CHK_IERR(phist_map_delete, myMap)

        # create a matrix
        A = DsparseMat_ptr()
        PYST_CHK_IERR(phist_Dcreate_matrix, A, myComm, b'spinSZ8')
        # get map
        PYST_CHK_IERR(phist_DsparseMat_get_domain_map, A, myMap)

        # create two mvecs
        vec1 = Dmvec_ptr()
        vec2 = Dmvec_ptr()
        vec3 = Dmvec_ptr()
        PYST_CHK_IERR(phist_Dmvec_create, vec1, myMap, 4)
        PYST_CHK_IERR(phist_Dmvec_create, vec2, myMap, 2)
        PYST_CHK_IERR(phist_Dmvec_view_block, vec1, vec3, 0,1 )

        # create some sdMats
        mat1 = DsdMat_ptr()
        mat2 = DsdMat_ptr()
        PYST_CHK_IERR(phist_DsdMat_create, mat1, 4, 2, myComm)
        PYST_CHK_IERR(phist_DsdMat_view_block, mat1, mat2, 0,1, 0,1)

        # put random data in all mvecs
        PYST_CHK_IERR(phist_Dmvec_random, vec1)
        PYST_CHK_IERR(phist_Dmvec_random, vec2)

        # QR-orthogonalize vec1
        try:
            PYST_CHK_IERR(phist_Dmvec_QR, vec2, mat2)
        except PYST_Exception as e:
            if e.value == -99:
                print('mvec_QR not available')
            else:
                raise

        # calculate vec1^T vec2
        PYST_CHK_IERR(phist_DmvecT_times_mvec, 1, vec1, vec2, 0, mat1)

        # calculate A*vec2
        PYST_CHK_IERR(phist_DsparseMat_times_mvec, 1, A, vec2, 0, vec3)

        # delete sdMats
        PYST_CHK_IERR(phist_DsdMat_delete, mat2)
        PYST_CHK_IERR(phist_DsdMat_delete, mat1)

        # delete mvecs
        PYST_CHK_IERR(phist_Dmvec_delete, vec3)
        PYST_CHK_IERR(phist_Dmvec_delete, vec2)
        PYST_CHK_IERR(phist_Dmvec_delete, vec1)

        # delete matrix
        PYST_CHK_IERR(phist_DsparseMat_delete, A)

        # delete comm
        PYST_CHK_IERR(phist_comm_delete, myComm)

        # finalize
        PYST_CHK_IERR(phist_kernels_finalize)

