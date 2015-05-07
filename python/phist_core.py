'''Python wrapper for phist_core library
'''

# imports
from __future__ import print_function
import ctypes as _ct
# we need to load phist_kernels and tools before
import phist_tools as _phist_tools
import phist_kernels as _phist_kernels


#--------------------------------------------------------------------------------
# load library
_phist_core = _ct.CDLL(name='libphist_core.so', mode=_ct.RTLD_GLOBAL)


#--------------------------------------------------------------------------------
# helper functions
from phist_kernels import _DeclareHelper
_declare = _DeclareHelper(lib=_phist_core)

from phist_kernels import _set


#--------------------------------------------------------------------------------
# some helper data types
c_int = _phist_kernels.c_int
c_int_p = _phist_kernels.c_int_p
_map_ptr = _phist_kernels.map_ptr

#--------------------------------------------------------------------------------
# data type dependent routines
for _varT in ('S', 'D', 'C', 'Z'):
    _prefix = 'phist_'+_varT

    # scalar data types
    _ST_ = getattr(_phist_kernels, _varT)
    _MT_ = {'S': _phist_kernels.S, 'D': _phist_kernels.D, 'C': _phist_kernels.S, 'Z': _phist_kernels.D}[_varT]
    _ST_p = _ct.POINTER(_ST_)
    _MT_p = _ct.POINTER(_MT_)
    _ST_pp = _ct.POINTER(_ST_p)
    _MT_pp = _ct.POINTER(_MT_p)

    # use types from kernels
    _sparseMat_ptr = getattr(_phist_kernels, _varT+'sparseMat_ptr')
    _sparseMat_ptr_p = getattr(_phist_kernels, _varT+'sparseMat_ptr_p')
    _mvec_ptr = getattr(_phist_kernels, _varT+'mvec_ptr')
    _mvec_ptr_p = getattr(_phist_kernels, _varT+'mvec_ptr_p')
    _sdMat_ptr = getattr(_phist_kernels, _varT+'sdMat_ptr')
    _sdMat_ptr_p = getattr(_phist_kernels, _varT+'sdMat_ptr_p')

    # op
    class _op(_ct.Structure):
        pass
    _op._fields_ = [('A', _ct.c_void_p),
                    ('range_map', _map_ptr),
                    ('domain_map', _map_ptr),
                    ('apply', _ct.c_void_p),
                    ('applyT', _ct.c_void_p),
                    ('apply_shifted', _ct.c_void_p)]

    _set(_varT+'op', _op)
    _op_ptr = _ct.POINTER(_op)
    _set(_varT+'op_ptr', _op_ptr)
    _op_ptr_p = _ct.POINTER(_op_ptr)
    _set(_varT+'op_ptr_p', _op_ptr_p)

    # from phist_driver_operator_decl.h
    #void SUBR(op_wrap_sparseMat)(TYPE(op_ptr) op, TYPE(const_sparseMat_ptr) A, int* ierr);
    #void SUBR(op_identity)(TYPE(op_ptr) op, int* ierr);
    _declare(None, _prefix+'op_wrap_sparseMat', (_op_ptr, _sparseMat_ptr, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'op_identity', (_op_ptr, c_int_p), skip_if_missing=True)

    # from phist_orthog_decl.h
    #void SUBR(orthog)(TYPE(const_mvec_ptr) V, TYPE(mvec_ptr) W, TYPE(const_op_ptr) B, TYPE(sdMat_ptr) R1, TYPE(sdMat_ptr) R2, int numSweeps, int* rankVW, int* iflag);
    _declare(None, _prefix+'orthog', (_mvec_ptr, _mvec_ptr, _op_ptr, _sdMat_ptr, _sdMat_ptr, c_int, c_int_p, c_int_p), skip_if_missing=True)


#--------------------------------------------------------------------------------
# tests / sample code
if __name__ == '__main__':

    # also import phist_tools
    from phist_tools import *
    from phist_kernels import *
    import ctypes

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
    rnkV = ctypes.c_int()
    PYST_CHK_NEG_IERR(phist_Dorthog, None, vec1, None, mat1, None, 1, rnkV)

    # orthogonalize vec2 wrt. vec1
    rnkVW = ctypes.c_int()
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

