'''Python wrapper for phist_tools library, unit tests
'''

#--------------------------------------------------------------------------------
# imports
from __future__ import print_function
import ctypes as _ct

from phist_tools import *

try:
    import unittest2 as unittest
except ImportError:
    import unittest


#--------------------------------------------------------------------------------
# tests

# for testing PYST_CHK_IERR
def some_fcn(bla, ierr):
    ierr[0] = bla
    
class PhistToolsTest(unittest.TestCase):
    
    @unittest.skip("this test is skipped (to show how that can be done)")
    def test_skipped(self):
        self.fail("shouldn't happen")

    def test_print_return_code(self):
        print('Just print some return codes')
        print(0, phist_retcode2str(0))
        print(PYST_Exception(-1))

    def test_check_return_code(self):
        print('Run a function with PYST_CHK_IERR')
        PYST_CHK_IERR(some_fcn, 0)

    def test_check_failing_function(self):
        print('Run a function that fails with PYST_CHK_IERR and print the error')
        try:
            PYST_CHK_IERR(some_fcn, -1)
            self.fail('PYST_CHK_IERR should have raised a PYST_Exception!')
        except PYST_Exception as e:
            print(e)

    def test_get_kernel_lib(self):
        print('Get the configured phist kernel lib')
        print(phist_kernel_lib())
