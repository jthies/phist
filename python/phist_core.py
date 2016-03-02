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
class _DeclareHelper:
    '''a helper class to set appropriate restype and argtypes and make the function available in the current scope'''

    def __init__(self, lib):
        self.lib = lib

    def __call__(self, restype, fcn_name, argtypes, skip_if_missing=False):
        '''actually sets restype and argtypes...'''
        if skip_if_missing and not hasattr(self.lib, fcn_name):
            return
        getattr(self.lib, fcn_name).restype = restype
        getattr(self.lib, fcn_name).argtypes = argtypes
        globals()[fcn_name] = getattr(self.lib, fcn_name)

# declare a function for all four data types, analogous to SUBR(func)
_declare = _DeclareHelper(lib=_phist_core)

# declare a global name
def _set(varName, value):
    globals()[varName] = value

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

    # linearOp
    class _linearOp(_ct.Structure):
        pass
    _linearOp._fields_ = [('A', _ct.c_void_p),
                    ('range_map', _map_ptr),
                    ('domain_map', _map_ptr),
                    ('aux',   _ct.c_void_p),
                    ('apply', _ct.c_void_p),
                    ('applyT', _ct.c_void_p),
                    ('apply_shifted', _ct.c_void_p)]

    _set(_varT+'linearOp', _linearOp)
    _linearOp_ptr = _ct.POINTER(_linearOp)
    _set(_varT+'linearOp_ptr', _linearOp_ptr)
    _linearOp_ptr_p = _ct.POINTER(_linearOp_ptr)
    _set(_varT+'linearOp_ptr_p', _linearOp_ptr_p)

    # from phist_driver_operator_decl.h
    #void SUBR(linearOp_wrap_sparseMat)(TYPE(linearOp_ptr) op, TYPE(const_sparseMat_ptr) A, int* ierr);
    #void SUBR(linearOp_identity)(TYPE(linearOp_ptr) op, int* ierr);
    _declare(None, _prefix+'linearOp_wrap_sparseMat', (_linearOp_ptr, _sparseMat_ptr, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'linearOp_identity', (_linearOp_ptr, c_int_p), skip_if_missing=True)

    # from phist_orthog_decl.h
    #void SUBR(orthog)(TYPE(const_mvec_ptr) V, TYPE(mvec_ptr) W, TYPE(const_linearOp_ptr) B, TYPE(sdMat_ptr) R1, TYPE(sdMat_ptr) R2, int numSweeps, int* rankVW, int* iflag);
    _declare(None, _prefix+'orthog', (_mvec_ptr, _mvec_ptr, _linearOp_ptr, _sdMat_ptr, _sdMat_ptr, c_int, c_int_p, c_int_p), skip_if_missing=True)


