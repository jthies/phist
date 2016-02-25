'''Python wrapper for phist_core library, unit tests
'''

#--------------------------------------------------------------------------------
# imports
from __future__ import print_function
import ctypes as _ct

from phist_tools import *
from phist_kernels import *
from phist_core import *

try:
    import unittest2 as unittest
except ImportError:
    import unittest


#--------------------------------------------------------------------------------
# tests
    
class PhistCoreTest(unittest.TestCase):
    
    def test_orthog(self):

        # initialize
        PYST_CHK_IERR(phist_kernels_init)

        # create comm
        myComm = comm_ptr()
        PYST_CHK_IERR(phist_comm_create, myComm)

        # create map
        nglobal = gidx(100)
        myMap = map_ptr()
        PYST_CHK_IERR(phist_map_create, myMap, myComm, nglobal)

        # create two mvecs
        vec1 = Dmvec_ptr()
        vec2 = Dmvec_ptr()
        PYST_CHK_IERR(phist_Dmvec_create, vec1, myMap, 10)
        PYST_CHK_IERR(phist_Dmvec_create, vec2, myMap, 2)

        # create some sdMats
        mat1 = DsdMat_ptr()
        mat2 = DsdMat_ptr()
        mat3 = DsdMat_ptr()
        PYST_CHK_IERR(phist_DsdMat_create, mat1, 10, 10, myComm)
        PYST_CHK_IERR(phist_DsdMat_create, mat2, 2, 2, myComm)
        PYST_CHK_IERR(phist_DsdMat_view_block, mat1, mat3, 0,9, 0,1)

        # put random data in all mvecs
        PYST_CHK_IERR(phist_Dmvec_random, vec1)
        PYST_CHK_IERR(phist_Dmvec_random, vec2)

        # QR-orthogonalize vec1
        rnkV = _ct.c_int()
        PYST_CHK_NEG_IERR(phist_Dorthog, None, vec1, None, mat1, None, 1, rnkV)

        # orthogonalize vec2 wrt. vec1
        rnkVW = _ct.c_int()
        PYST_CHK_NEG_IERR(phist_Dorthog, vec1, vec2, None, mat2, mat3, 3, rnkVW)

        # delete sdMats
        PYST_CHK_IERR(phist_DsdMat_delete, mat3)
        PYST_CHK_IERR(phist_DsdMat_delete, mat2)
        PYST_CHK_IERR(phist_DsdMat_delete, mat1)

        # delete mvecs
        PYST_CHK_IERR(phist_Dmvec_delete, vec2)
        PYST_CHK_IERR(phist_Dmvec_delete, vec1)

        # delete map
        PYST_CHK_IERR(phist_map_delete, myMap)

        # delete comm
        PYST_CHK_IERR(phist_comm_delete, myComm)

        # finalize
        PYST_CHK_IERR(phist_kernels_finalize)

