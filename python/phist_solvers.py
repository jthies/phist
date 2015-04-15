'''Python wrapper for phist_solvers library
'''

# imports
from __future__ import print_function
import ctypes as _ct
# we need to load kernels, tools and core
import phist_tools as _phist_tools
import phist_kernels as _phist_kernels
import phist_core as _phist_core


#--------------------------------------------------------------------------------
# load library
_phist_solvers = _ct.CDLL(name='libphist_solvers.so', mode=_ct.RTLD_GLOBAL)


#--------------------------------------------------------------------------------
# helper functions
from phist_kernels import _DeclareHelper
_declare = _DeclareHelper(lib=_phist_solvers)

from phist_kernels import _set


#--------------------------------------------------------------------------------
# some helper data types
c_int = _phist_kernels.c_int
c_int_p = _phist_kernels.c_int_p

#--------------------------------------------------------------------------------
# data type dependent routines
for _varT in ('S', 'D', 'C', 'Z'):
    _prefix = 'phist_'+_varT

    # scalar data types
    _ST_ = getattr(_phist_kernels, _varT)
    _MT_ = {'S': _phist_kernels.S, 'D': _phist_kernels.D, 'C': _phist_kernels.S, 'Z': _phist_kernels.D}[_varT]
    _CT_ = {'S': _phist_kernels.C, 'D': _phist_kernels.Z, 'C': _phist_kernels.C, 'Z': _phist_kernels.Z}[_varT]
    _ST_p = _ct.POINTER(_ST_)
    _CT_p = _ct.POINTER(_CT_)
    _MT_p = _ct.POINTER(_MT_)
    _ST_pp = _ct.POINTER(_ST_p)
    _MT_pp = _ct.POINTER(_MT_p)

    # use types from kernels
    _mvec_ptr = getattr(_phist_kernels, _varT+'mvec_ptr')
    _sdMat_ptr = getattr(_phist_kernels, _varT+'sdMat_ptr')

    # use types from core
    _op_ptr = getattr(_phist_kernels, _varT+'op_ptr')

    # from phist_subspacejada_decl.h
    #void SUBR(subspacejada)(
    #        TYPE(const_op_ptr) A_op,  TYPE(const_op_ptr) B_op,
    #        TYPE(const_mvec_ptr) v0,  eigSort_t which,
    #        _MT_ tol,                 int nEig,
    #        int* nIter,               int blockDim,
    #        int minBase,              int maxBase,
    #        int innerBlockDim,        int innerMaxBase,
    #        int initialShiftIter,     _ST_ initialShift,
    #        bool innerIMGS,           bool innerGMRESabortAfterFirstConverged,
    #        bool symmetric,
    #        TYPE(mvec_ptr) Q,         TYPE(sdMat_ptr) R,
    #        _CT_* ev,                 _MT_* resNorm,
    #        int* ierr);
    _declare(None, _prefix+'subspacejada', (_op_ptr, _op_ptr,
                                            _mvec_ptr, _phist_tools.eigSort_t,
                                            _MT_, c_int,
                                            c_int_p, c_int,
                                            c_int, c_int,
                                            c_int, c_int,
                                            c_int, _ST_,
                                            _ct.c_bool, _ct.c_bool,
                                            _ct.c_bool,
                                            _mvec_ptr, _sdMat_ptr,
                                            _CT_p, _MT_p, c_int_p), skip_if_missing=True)

    # from phist_jdqr_decl.h
    #void SUBR(jdqr)(TYPE(const_op_ptr) A_op, TYPE(const_op_ptr) B_op,
    #                TYPE(mvec_ptr) X, _ST_* evals, _MT_* resid, int* is_cmplx,
    #                phist_jadaOpts_t options, int* num_eigs, int* num_iters,
    #                int* ierr);
    #_declare(None, _prefix+'jdqr', (_op_ptr, _op_ptr,
    #                                _mvec_ptr, _ST_p, _MT_p, c_int_p,
    #                                phist_jadaOpts_t, c_int_p, c_int_p,
    #                                c_int_p), skip_if_missing=True)

#--------------------------------------------------------------------------------
# tests / sample code
if __name__ == '__main__':

    # also import phist_tools
    from phist_tools import *
    from phist_kernels import *
    from phist_core import *

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
    opA = Dop()
    opB = Dop_ptr()
    PYST_CHK_IERR(phist_Dop_wrap_sparseMat, opA, A)

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
    nIter = c_int(1000)
    resNorm = (D*11)()
    ev = (Z*11)()
    res = PYST_CHK_NEG_IERR(phist_Dsubspacejada, opA, opB, v0, eigSort_SR, 1.e-8, 10, nIter, 2, 20, 40, 2, 8, 0, 0., False, True, True, vecQ, matR, ev, resNorm)
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

