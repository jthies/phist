'''Python wrapper for phist_tools library
'''

#--------------------------------------------------------------------------------
# imports
from __future__ import print_function
import ctypes as _ct

#--------------------------------------------------------------------------------
# load the library
# try to preload the MPI library, WARNING: this may fail if there are multiple
# MPI libraries in the PATH!
try:
    _mpi = _ct.CDLL(name='libmpi_mt.so', mode=_ct.RTLD_GLOBAL)
except:
    try:
        _mpi = _ct.CDLL(name='libmpi.so', mode=_ct.RTLD_GLOBAL)
    except:
        pass
_phist_tools = _ct.CDLL(name='libphist_tools.so', mode=_ct.RTLD_GLOBAL)



#--------------------------------------------------------------------------------
# from phist_tools.h
#retcode2str
_phist_tools.phist_retcode2str.argtypes = (_ct.c_int,)
_phist_tools.phist_retcode2str.restype = _ct.c_char_p
_phist_tools.phist_retcode2str.__doc__ = 'convert phist return code to string (e.g. "ERROR")'
phist_retcode2str = _phist_tools.phist_retcode2str
#phist_kernel_lib
_phist_tools.phist_kernel_lib.argtypes = list()
_phist_tools.phist_kernel_lib.restype = _ct.c_char_p
_phist_tools.phist_kernel_lib.__doc__ = 'return configured kernel library'
phist_kernel_lib = _phist_tools.phist_kernel_lib



#--------------------------------------------------------------------------------
# from phist_enums.h
# eigSort
eigSort_t = _ct.c_uint
eigSort_NO_EIGSORT    = _ct.c_uint(0)
eigSort_LM      = _ct.c_uint(1)
eigSort_SM      = _ct.c_uint(2)
eigSort_LR      = _ct.c_uint(3)
eigSort_SR      = _ct.c_uint(4)
eigSort_TARGET  = _ct.c_uint(5)
_phist_tools.eigSort2str.argtypes = (eigSort_t,)
_phist_tools.eigSort2str.resType = _ct.c_char_p
eigSort2str = _phist_tools.eigSort2str
_phist_tools.str2eigSort.argtypes = (_ct.c_char_p,)
_phist_tools.str2eigSort.resType = eigSort_t
str2eigSort = _phist_tools.str2eigSort

# eigExtr
eigExtr_t = _ct.c_uint
eigExtr_STANDARD = _ct.c_uint(0)
eigSort_HARMONIC = _ct.c_uint(1)
eigSort_INVALID_EIGEXTR_T = _ct.c_uint(99)
_phist_tools.eigExtr2str.argtypes = (eigExtr_t,)
_phist_tools.eigExtr2str.resType = _ct.c_char_p
eigExtr2str = _phist_tools.eigExtr2str
_phist_tools.str2eigExtr.argtypes = (_ct.c_char_p,)
_phist_tools.str2eigExtr.resType = eigExtr_t
str2eigExtr = _phist_tools.str2eigExtr

# linSolv
linSolv_t = _ct.c_uint
linSolv_DO_NOTHING   = _ct.c_uint(0)
linSolv_GMRES        = _ct.c_uint(1)
linSolv_MINRES       = _ct.c_uint(2)
linSolv_CARP_CG      = _ct.c_uint(3)
linSolv_USER_DEFINED = _ct.c_uint(98)
_phist_tools.linSolv2str.argtypes = (linSolv_t,)
_phist_tools.linSolv2str.resType = _ct.c_char_p
linSolv2str = _phist_tools.linSolv2str
_phist_tools.str2linSolv.argtypes = (_ct.c_char_p,)
_phist_tools.str2linSolv.resType = linSolv_t
str2linSolv = _phist_tools.str2linSolv
# matSym
matSym_t = _ct.c_uint
matSym_GENERAL            = _ct.c_uint(0)
matSym_HERMITIAN          = _ct.c_uint(1)
matSym_COMPLEX_SYMMETRIC  = _ct.c_uint(2)
matSym_PATTERN_SYMMETRIC  = _ct.c_uint(3)


#--------------------------------------------------------------------------------
# PYST_CHK_IERR function wrapper
# Exception raised by PYST_CHK_IERR 
class PYST_Exception(Exception):
    '''Exception raised by PYST_CHK_IERR'''
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return '%s (%s)' % (self.value, phist_retcode2str(self.value))

# useful wrapper for calling phist functions
def PYST_CHK_IERR(phist_fcn, *args):
    '''Run the phist function phist_fcn with given arguments and append and check the ierr flag automatically'''
    err = _ct.c_int()
    args_with_err = [v for v in args]
    args_with_err.append(_ct.pointer(err))
    phist_fcn(*args_with_err)
    if err.value != 0:
        raise PYST_Exception(err.value)

# for functions that may return positive values
def PYST_CHK_NEG_IERR(phist_fcn, *args):
    '''Run the phist function phist_fcn with given arguments and append, check and return the ierr flag'''
    err = _ct.c_int()
    args_with_err = [v for v in args]
    args_with_err.append(_ct.pointer(err))
    phist_fcn(*args_with_err)
    if err.value < 0:
        raise PYST_Exception(err.value)
    return err.value

#--------------------------------------------------------------------------------
# tests
if __name__ == '__main__':

    print('Just print some return codes')
    print(0, phist_retcode2str(0))
    print(PYST_Exception(-1))

    # test PYST_CHK_IERR
    def test_fcn(bla, ierr):
        ierr[0] = bla

    print('Run a function with PYST_CHK_IERR')
    PYST_CHK_IERR(test_fcn, 0)

    print('Run a function that fails with PYST_CHK_IERR and print the error')
    try:
        PYST_CHK_IERR(test_fcn, -1)
        raise RuntimeError('PYST_CHK_IERR should have raised a PYST_Exception!')
    except PYST_Exception as e:
        print(e)

    print('Get the configured phist kernel lib')
    print(phist_kernel_lib())
