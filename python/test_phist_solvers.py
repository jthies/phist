'''Python wrapper for phist_solvers library, unit tests
'''

#--------------------------------------------------------------------------------
# imports
from __future__ import print_function
import ctypes as _ct

from phist_tools import *
from phist_kernels import *
from phist_core import *
from phist_solvers import *

try:
    import unittest2 as unittest
except ImportError:
    import unittest


#--------------------------------------------------------------------------------
# tests

class PhistSolversTest(unittest.TestCase):
    
#    @unittest.skip("this test is not running because of missing symbols like DlinearOp")
    def test_subspacejada_on_spinSZ(self):

        # initialize
        PYST_CHK_IERR(phist_kernels_init)

        # create comm
        myComm = comm_ptr()
        PYST_CHK_IERR(phist_comm_create, myComm)

        # create a matrix
        A = DsparseMat_ptr()
        PYST_CHK_IERR(phist_Dcreate_matrix, A, myComm, b'spinSZ8')

        # get map
        myMap = map_ptr()
        PYST_CHK_IERR(phist_DsparseMat_get_domain_map, A, myMap)

        # wrap A in operator
        opA = DlinearOp()
        opB = DlinearOp_ptr()
        PYST_CHK_IERR(phist_DlinearOp_wrap_sparseMat, opA, A)

        # create two mvecs
        v0 = Dmvec_ptr()
        vecQ = Dmvec_ptr()
        PYST_CHK_IERR(phist_Dmvec_create, v0, myMap, 1)
        PYST_CHK_IERR(phist_Dmvec_random, v0)
        PYST_CHK_IERR(phist_Dmvec_create, vecQ, myMap, 11)

        # create some sdMats
        matR = DsdMat_ptr()
        PYST_CHK_IERR(phist_DsdMat_create, matR, 11, 11, myComm)


        # call subspacejada
        jadaOpts = phist_jadaOpts_t()
        jadaOpts.numEigs = 10
        jadaOpts.v0 = v0
        jadaOpts.which = eigSort_SR
        jadaOpts.how = eigExtr_STANDARD
        jadaOpts.convTol = 1.e-8
        jadaOpts.maxIters = 50
        nIter = c_int()
        nConv = c_int()
        resNorm = (D*11)()
        ev = (Z*11)()
        res = PYST_CHK_NEG_IERR(phist_Dsubspacejada, opA, opB, jadaOpts, vecQ, matR, ev, resNorm, nConv, nIter)
        print('eigenvalues:  ', [ev[i][0] for i in range(10)])
        print('ritz residual:', [resNorm[i] for i in range(10)])


        # delete sdMats
        PYST_CHK_IERR(phist_DsdMat_delete, matR)

        # delete mvecs
        PYST_CHK_IERR(phist_Dmvec_delete, vecQ)
        PYST_CHK_IERR(phist_Dmvec_delete, v0)

        # delete matrix
        PYST_CHK_IERR(phist_DsparseMat_delete, A)

        # delete comm
        PYST_CHK_IERR(phist_comm_delete, myComm)

        # finalize
        PYST_CHK_IERR(phist_kernels_finalize)


