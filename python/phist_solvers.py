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
c_double = _ct.c_double
c_void_p = _ct.c_void_p

#--------------------------------------------------------------------------------
# data type independent structs/routines
# from phist_jadaOpts.h
#typedef struct phist_jadaOpts_t { ... } phist_jadaOpts_t;
class phist_jadaOpts_t(_ct.Structure):
    '''This struct can be used to consistently pass
       parameters to our various Jacobi-Davidson methods.


       // what do you want to compute?
       int numEigs; //! howmany eigenpairs are sought?
       eigSort_t which; //! LM, SM, LR, SR, or TARGET
       double convTol; //! convergence tolerance for eigenvalues
       matSym_t symmetry; //! Symmetry properties of the matrix
       // JaDa configuration
       int maxIters; //! maximum iterations allowed
       int blockSize; //! only for block methods (subspacejada)
       int minBas; //! number of vectors retained upon restart
       int maxBas; //! maximum number of vectors allowed in the basis
       // how should JaDa start up?
       void* v0; //! can be used to pass in a start-up vector(-space) (can have any number of
                 //! columns). v0 is assumed to be orthonormal.
       int arno; //! 0: no Arnoldi steps. 1: minBas Arnoldi steps to start up.
       double initialShift_r; //! can be used to start with an initial shift
                              //! (ignored if arno!=0)
       double initialShift_i; //! imaginary part of initial shift
       
       int initialShiftIters; // perform given number of iterations with a fixed shift
       // inner solver configuration
       linSolv_t innerSolvType; /*! GMRES, MINRES, CARP_CG, USER_DEFINED currently supported.
                                 * If set to USER_DEFINED, you have to provide the customSolver*
                                 * interface below.
                                 */
       int innerSolvBlockSize;
       int innerSolvMaxBas;
       int innerSolvMaxIters;
       int innerSolvRobust; /*! extra effort to get good jada updates
                             * (in practice this may mean a more accurate orthogonalization etc.)
                             */
       int innerSolvStopAfterFirstConverged;
       
       
         //! pointer to solver object if innerSolvType==USER_DEFINED
         void* customSolver;
       
         //! this function is used instead of phist_jadaCorrectionSolver_run if innerSolvType is USER_DEFINED.
         //! For jdqr or subspacejada with block size 1 it is enough to implement the simpler interface below,
         //! we first check in those cases if that interface is set before checking for this one.
         //! note that the scalar arguments are always passed in as doubles so that a single function pointer can
         //! be used in this untyped struct.
         void (*customSolver_run)(         void*       customSolverData,
                                           void const* A_op,      void const*    B_op,
                                           void const* Qtil,      void const*    BQtil,
                                           const double sigma_r[], const double sigma_i[],
                                           void const* res,        const int resIndex[],
                                           const double tol[],       int maxIter,
                                           void* t,                int robust,
                                           int abortAfterFirstConvergedInBlock,
                                           int * iflag);
       
         //! simplified interface if only single-vector jdqr or subspacejada is used.
         void (*customSolver_run1)(        void*  customSolverData,
                                           void const* A_op,     void const*    B_op,
                                           void const* Qtil,     void const*    BQtil,
                                           const double sigma,   double sigma_i,
                                           void const* res,      double tol,
                                           int maxIter,          void* t,
                                           int robust,
                                           int * iflag);
    '''

    _fields_ = [("numEigs",             c_int),
                ("which",               _phist_tools.eigSort_t),
                ("convTol",             c_double),
                ("symmetry",            _phist_tools.matSym_t),
                ("maxIters",            c_int),
                ("blockSize",           c_int),
                ("minBas",              c_int),
                ("maxBas",              c_int),
                ("v0",                  c_void_p),
                ("arno",                c_int),
                ("initialShift_r",      c_double),
                ("initialShift_i",      c_double),
                ("initialShiftIters",   c_int),
                ("innerSolvType",       _phist_tools.linSolv_t),
                ("innerSolvBlockSize",  c_int),
                ("innerSolvMaxBas",     c_int),
                ("innerSolvMaxIters",   c_int),
                ("innerSolvRobust",     c_int),
                ("innerSolvStopAfterFirstConverged", c_int),
                ("customSolver",        c_void_p),
                ("customSolver_run",    c_void_p),
                ("customSolver_run1",  c_void_p),]

    def __init__(self):
        super(_ct.Structure,self).__init__()
        _phist_solvers.phist_jadaOpts_setDefaults(_ct.byref(self))


# from phist_jadaOpts.h
phist_jadaOpts_t_p = _ct.POINTER(phist_jadaOpts_t)
#void phist_jadaOpts_setDefaults(phist_jadaOpts_t *opts);
_phist_solvers.phist_jadaOpts_setDefaults.argtypes = (phist_jadaOpts_t_p,)
_phist_solvers.phist_jadaOpts_setDefaults.restype = None
phist_jadaOpts_setDefaults = _phist_solvers.phist_jadaOpts_setDefaults


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
    #void SUBR(subspacejada)( TYPE(const_op_ptr) A_op,  TYPE(const_op_ptr) B_op,
    #                         phist_jadaOpts_t opts,
    #                         TYPE(mvec_ptr) Q,         TYPE(sdMat_ptr) R,
    #                         _CT_* ev,                 _MT_* resNorm,
    #                         int *nConv,                int *nIter,
    #                         int* iflag);
    _declare(None, _prefix+'subspacejada', (_op_ptr, _op_ptr,
                                            phist_jadaOpts_t,
                                            _mvec_ptr, _sdMat_ptr,
                                            _CT_p, _MT_p, 
                                            c_int_p, c_int_p,
                                            c_int_p), skip_if_missing=True)

    # from phist_jdqr_decl.h
    #void SUBR(jdqr)(TYPE(const_op_ptr) A_op, TYPE(const_op_ptr) B_op,
    #                TYPE(mvec_ptr) X, TYPE(mvec_ptr) Q, TYPE(sdMat_ptr) R,
    #                _ST_* evals, _MT_* resid, int* is_cmplx,
    #                phist_jadaOpts_t options, int* num_eigs, int* num_iters,
    #                int* iflag);
    _declare(None, _prefix+'jdqr', (_op_ptr, _op_ptr,
                                    _mvec_ptr, _mvec_ptr, _sdMat_ptr,
                                    _ST_p, _MT_p, c_int_p,
                                    phist_jadaOpts_t, c_int_p, c_int_p,
                                    c_int_p), skip_if_missing=True)

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
    jadaOpts = phist_jadaOpts_t()
    jadaOpts.numEigs = 10
    jadaOpts.v0 = v0
    jadaOpts.which = eigSort_SR
    jadaOpts.convTol = 1.e-8
    jadaOpts.maxIters = 1000
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


