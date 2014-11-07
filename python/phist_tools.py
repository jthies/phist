'''Python wrapper for phist_tools library
'''

#--------------------------------------------------------------------------------
# imports
import ctypes as _ct

#--------------------------------------------------------------------------------
# load the library
_phist_tools = _ct.CDLL(name='libphist_tools.so', mode=_ct.RTLD_GLOBAL)

#--------------------------------------------------------------------------------
# phist_retcode2str
_phist_tools.phist_retcode2str.argtypes = (_ct.c_int,)
_phist_tools.phist_retcode2str.restype = _ct.c_char_p
_phist_tools.phist_retcode2str.__doc__ = 'convert phist return code to string (e.g. "ERROR")'

phist_retcode2str = _phist_tools.phist_retcode2str



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


#--------------------------------------------------------------------------------
# tests
if __name__ == '__main__':

    print 'Just print some return codes'
    print 0, phist_retcode2str(0)
    print PYST_Exception(-1)

    # test PYST_CHK_IERR
    def test_fcn(bla, ierr):
        ierr[0] = bla

    print 'Run a function with PYST_CHK_IERR'
    PYST_CHK_IERR(test_fcn, 0)

    print 'Run a function that fails with PYST_CHK_IERR and print the error'
    try:
        PYST_CHK_IERR(test_fcn, -1)
        raise RuntimeError('PYST_CHK_IERR should have raised a PYST_Exception!')
    except PYST_Exception as e:
        print e
