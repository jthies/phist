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
    
    def test_subspacejada_on_spinSZ(self):

        # initialize
        PYST_CHK_IERR(phist_kernels_init)

        # create comm
        myComm = comm_ptr()
        PYST_CHK_IERR(phist_comm_create, myComm)
        
        numEigs=c_int(8)
        nv1    =      8
        nv     =      12
        
        minBas=c_int(20)
        maxBas=c_int(30)
        
        blockSize1=c_int(1)
        blockSize2=c_int(4)

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
        #ask subspacejada to return a larger basis Q so we can restart with it
        PYST_CHK_IERR(phist_Dmvec_create, vecQ, myMap, minBas)

        # create some sdMats
        matR = DsdMat_ptr()
        PYST_CHK_IERR(phist_DsdMat_create, matR, minBas, minBas, myComm)


        # call subspacejada
        jadaOpts = phist_jadaOpts_t()
        jadaOpts.symmetry=matSym_HERMITIAN
        jadaOpts.numEigs = numEigs
        jadaOpts.blockSize=blockSize1
        jadaOpts.minBas=minBas
        jadaOpts.maxBas=maxBas
        jadaOpts.v0 = v0
        jadaOpts.which = eigSort_SR
        jadaOpts.how = eigExtr_STANDARD
        jadaOpts.convTol = 1.e-8
        jadaOpts.maxIters = 20
        jadaOpts.innerSolvType=linSolv_MINRES
        nIter = c_int()
        nConv = c_int()
        resNorm = (D*nv)()
        ev = (Z*nv)()
        res = PYST_CHK_IERR(phist_Dsubspacejada, opA, opB, jadaOpts, vecQ, matR, ev, resNorm, nConv, nIter)
        print('eigenvalues after first run:  ', [ev[i][0] for i in range(nv)])
        print('ritz residual:', [resNorm[i] for i in range(nv)])
        
        # Test that restarting with the computed Q gives all eigenvalues in a few iterations.
        # We need to do a few iteration because of the limited 'look ahead' strategy in the
        # solver, it only checks the residuals of the next2*blockSize eigenpairs or so. So to
        # speed things up a little we increase the block size to 4 for this experiment. 
        jadaOpts.maxIters=4
        jadaOpts.blockSize=blockSize2
        jadaOpts.innerSolvBlockSize=blockSize2
        jadaOpts.v0=vecQ
        ev2 = (Z*nv)()
        res = PYST_CHK_IERR(phist_Dsubspacejada, opA, opB, jadaOpts, vecQ, matR, ev2, resNorm, nConv, nIter)

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


