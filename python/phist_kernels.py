'''Python wrapper for phist_kernels library
'''

# imports
import ctypes as _ct
# we need to load phist_tools before
import phist_tools as _phist_tools

#--------------------------------------------------------------------------------
# load library
print 'Using ', _phist_tools.phist_kernel_lib(), ' kernels'
_phist_kernels = _ct.CDLL(name='libphist_kernels_'+_phist_tools.phist_kernel_lib()+'.so', mode=_ct.RTLD_GLOBAL)

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

_declare = _DeclareHelper(lib=_phist_kernels)

def _set(varName, value):
    globals()[varName] = value


#--------------------------------------------------------------------------------
# data types (mostly void pointers)
c_int = _ct.c_int
c_int_p = _ct.POINTER(_ct.c_int)
c_char_p = _ct.c_char_p
class comm_ptr(_ct.c_void_p):
    pass
comm_ptr_p = _ct.POINTER(comm_ptr)
class map_ptr(_ct.c_void_p):
    pass
map_ptr_p = _ct.POINTER(map_ptr)

# from phist_driver_utils.h
# int phist_sizeof_lidx()
_declare(c_int, 'phist_sizeof_lidx', list())
# int phist_sizeof_gidx()
_declare(c_int, 'phist_sizeof_gidx', list())

# detect gidx_t and lidx_t
def _c_integer_type_from_size(size):
    if size == 4:
        return _ct.c_int32
    elif size == 8:
        return _ct.c_int64
    raise NotImplementedError('Unknown integer size', size)

_lidx_size = phist_sizeof_lidx()
_gidx_size = phist_sizeof_gidx()
print 'lidx size', _lidx_size, 'gidx size', _gidx_size

class lidx(_c_integer_type_from_size(_lidx_size)):
    def __repr__(self):
        return 'lidx(%s)' % (self.value)

class gidx(_c_integer_type_from_size(_gidx_size)):
    def __repr__(self):
        return 'gidx(%s)' % (self.value)

lidx_p = _ct.POINTER(lidx)
gidx_p = _ct.POINTER(gidx)


#--------------------------------------------------------------------------------
# data type independent routines
# from phist_kernels.h

#void phist_kernels_init(int *argc, char*** argv, int* ierr);
_phist_kernels.phist_kernels_init.restype = None
_phist_kernels.phist_kernels_init.argtypes = (c_int_p, _ct.POINTER(_ct.POINTER(_ct.c_char_p)), c_int_p)
def phist_kernels_init(ierr):
    # construct dummy args, possibly dangerours!!
    argc = c_int(1)
    argv0 = _ct.c_char_p('dummyarg')
    argv = _ct.pointer(argv0)
    _phist_kernels.phist_kernels_init(_ct.byref(argc), _ct.byref(argv), ierr)

#void phist_kernels_finalize(int* ierr);
_declare(None, 'phist_kernels_finalize', (c_int_p,))

#void phist_comm_create(comm_ptr_t* comm, int* ierr);
_declare(None, 'phist_comm_create', (comm_ptr_p, c_int_p))

#void phist_comm_delete(comm_ptr_t comm, int* ierr);
_declare(None, 'phist_comm_delete', (comm_ptr, c_int_p))

#void phist_comm_get_rank(const_comm_ptr_t comm, int* rank, int* ierr);
_declare(None, 'phist_comm_get_rank', (comm_ptr, c_int_p, c_int_p))

#void phist_comm_get_size(const_comm_ptr_t comm, int* size, int* ierr);
_declare(None, 'phist_comm_get_size', (comm_ptr, c_int_p, c_int_p))

#void phist_map_create(map_ptr_t* map, const_comm_ptr_t comm, gidx_t nglob, int *ierr);
_declare(None, 'phist_map_create', (map_ptr_p, comm_ptr, gidx, c_int_p))

#void phist_map_delete(map_ptr_t map, int *ierr);
_declare(None, 'phist_map_delete', (map_ptr, c_int_p))

#void phist_map_get_comm(const_map_ptr_t map, const_comm_ptr_t* comm, int* ierr);
_declare(None, 'phist_map_get_comm', (map_ptr, comm_ptr_p, c_int_p))

#void phist_map_get_local_length(const_map_ptr_t map, lidx_t* nloc, int* ierr);
_declare(None, 'phist_map_get_local_length', (map_ptr, lidx_p, c_int_p))

#void phist_map_get_ilower(const_map_ptr_t map, gidx_t* ilower, int* ierr);
_declare(None, 'phist_map_get_ilower', (map_ptr, lidx_p, c_int_p))

#void phist_map_get_iupper(const_map_ptr_t map, gidx_t* iupper, int* ierr);
_declare(None, 'phist_map_get_iupper', (map_ptr, lidx_p, c_int_p))


#--------------------------------------------------------------------------------
# data type dependent routines
S = _ct.c_float
D = _ct.c_double
C = (S*2)
Z = (D*2)
S_p = _ct.POINTER(S)
D_p = _ct.POINTER(D)
C_p = _ct.POINTER(C)
Z_p = _ct.POINTER(Z)
S_pp = _ct.POINTER(S_p)
D_pp = _ct.POINTER(D_p)
C_pp = _ct.POINTER(C_p)
Z_pp = _ct.POINTER(Z_p)
for _varT in ('S', 'D', 'C', 'Z'):
    _prefix = 'phist_'+_varT

    # scalar data types
    _ST_ = globals()[_varT]
    _MT_ = {'S': S, 'D': D, 'C': S, 'Z': D}[_varT]
    _ST_p = _ct.POINTER(_ST_)
    _MT_p = _ct.POINTER(_MT_)
    _ST_pp = _ct.POINTER(_ST_p)
    _MT_pp = _ct.POINTER(_MT_p)

    # sparseMat
    class _sparseMat_ptr(_ct.c_void_p):
        pass
    _set(_varT+'sparseMat_ptr', _sparseMat_ptr)
    _sparseMat_ptr_p = _ct.POINTER(_sparseMat_ptr)
    _set(_varT+'sparseMat_ptr_p', _sparseMat_ptr_p)

    # mvec
    class _mvec_ptr(_ct.c_void_p):
        pass
    _set(_varT+'mvec_ptr', _mvec_ptr)
    _mvec_ptr_p = _ct.POINTER(_mvec_ptr)
    _set(_varT+'mvec_ptr_p', _mvec_ptr_p)

    # sdMat
    class _sdMat_ptr(_ct.c_void_p):
        pass
    _set(_varT+'sdMat_ptr', _sdMat_ptr)
    _sdMat_ptr_p = _ct.POINTER(_sdMat_ptr)
    _set(_varT+'sdMat_ptr_p', _sdMat_ptr_p)

    # from phist_driver_utils.h
    #void  SUBR(sparseMat_read)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t comm, char* filename, int* ierr);
    _declare(None, _prefix+'sparseMat_read', (_sparseMat_ptr_p, comm_ptr, c_char_p, c_int_p), skip_if_missing=True)

    #void  SUBR(create_matrix)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t comm, char* filename, int* ierr);
    _declare(None, _prefix+'create_matrix', (_sparseMat_ptr_p, comm_ptr, c_char_p, c_int_p), skip_if_missing=True)


    # from phist_kernels_decl.h
    #void SUBR(type_avail)(int* ierr);
    _declare(None, _prefix+'type_avail', (c_int_p,), skip_if_missing=True)

    #void  SUBR(sparseMat_read_mm)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t comm, const char* filename, int* ierr);
    #void  SUBR(sparseMat_read_hb)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t comm, const char* filename, int* ierr);
    #void  SUBR(sparseMat_read_bin)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t comm, const char* filename, int* ierr);
    _declare(None, _prefix+'sparseMat_read_mm', (_sparseMat_ptr_p, comm_ptr, c_char_p, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sparseMat_read_hb', (_sparseMat_ptr_p, comm_ptr, c_char_p, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sparseMat_read_bin', (_sparseMat_ptr_p, comm_ptr, c_char_p, c_int_p), skip_if_missing=True)

    #void SUBR(sparseMat_get_row_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* ierr);
    #void SUBR(sparseMat_get_col_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* ierr);
    #void SUBR(sparseMat_get_range_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* ierr);
    #void SUBR(sparseMat_get_domain_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* ierr);
    _declare(None, _prefix+'sparseMat_get_row_map', (_sparseMat_ptr, map_ptr_p, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sparseMat_get_col_map', (_sparseMat_ptr, map_ptr_p, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sparseMat_get_range_map', (_sparseMat_ptr, map_ptr_p, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sparseMat_get_domain_map', (_sparseMat_ptr, map_ptr_p, c_int_p), skip_if_missing=True)

    #void SUBR(mvec_create)(TYPE(mvec_ptr)* V, const_map_ptr_t map, int nvec, int* ierr);
    #void SUBR(mvec_create_view)(TYPE(mvec_ptr)* V, const_map_ptr_t map, _ST_* values, lidx_t lda, int nvec, int* ierr);
    #void SUBR(sdMat_create)(TYPE(sdMat_ptr)* M, int nrows, int ncols, const_comm_ptr_t comm, int* ierr);
    _declare(None, _prefix+'mvec_create', (_mvec_ptr_p, map_ptr, c_int, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'mvec_create_view', (_mvec_ptr_p, map_ptr, _ST_p, lidx, c_int, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sdMat_create', (_sdMat_ptr_p, c_int, c_int, comm_ptr, c_int_p), skip_if_missing=True)

    #void SUBR(sparseMat_delete)(TYPE(sparseMat_ptr) A, int* ierr);
    #void SUBR(mvec_delete)(TYPE(mvec_ptr) V, int* ierr);
    #void SUBR(sdMat_delete)(TYPE(sdMat_ptr) M, int* ierr);
    _declare(None, _prefix+'sparseMat_delete', (_sparseMat_ptr, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'mvec_delete', (_mvec_ptr, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sdMat_delete', (_sdMat_ptr, c_int_p), skip_if_missing=True)

    #void SUBR(mvec_my_length)(TYPE(const_mvec_ptr) V, lidx_t* len, int* ierr);
    #void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) V, const_map_ptr_t* map, int* ierr);
    #void SUBR(mvec_get_comm)(TYPE(const_mvec_ptr) V, const_comm_ptr_t* comm, int* ierr);
    #void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) V, int* nvec, int* ierr);
    _declare(None, _prefix+'mvec_my_length', (_mvec_ptr, lidx_p, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'mvec_get_map', (_mvec_ptr, map_ptr_p, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'mvec_get_comm', (_mvec_ptr, comm_ptr_p, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'mvec_num_vectors', (_mvec_ptr, c_int_p, c_int_p), skip_if_missing=True)


    #void SUBR(mvec_extract_view)(TYPE(mvec_ptr) V, _ST_** val, lidx_t* lda, int* ierr);
    #void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) M, int* nrows, int* ierr);
    #void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) M, int* ncols, int* ierr);
    #void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) M, _ST_** val, lidx_t* lda, int* ierr);
    _declare(None, _prefix+'mvec_extract_view', (_mvec_ptr, _ST_pp, lidx_p, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sdMat_get_nrows', (_sdMat_ptr, c_int_p, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sdMat_get_ncols', (_sdMat_ptr, c_int_p, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sdMat_extract_view', (_sdMat_ptr, _ST_pp, lidx_p, c_int_p), skip_if_missing=True)

    #void SUBR(mvec_to_device)(TYPE(mvec_ptr) V, int* ierr);
    #void SUBR(mvec_from_device)(TYPE(mvec_ptr) V, int* ierr);
    #void SUBR(sdMat_to_device)(TYPE(sdMat_ptr) M, int* ierr);
    #void SUBR(sdMat_from_device)(TYPE(sdMat_ptr) M, int* ierr);
    _declare(None, _prefix+'mvec_to_device', (_mvec_ptr, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'mvec_from_device', (_mvec_ptr, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sdMat_to_device', (_sdMat_ptr, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sdMat_from_device', (_sdMat_ptr, c_int_p), skip_if_missing=True)

    #void SUBR(mvec_to_mvec)(TYPE(const_mvec_ptr) v_in, TYPE(mvec_ptr) v_out, int* ierr);
    #void SUBR(mvec_view_block)(TYPE(mvec_ptr) V,TYPE(mvec_ptr)* Vblock,int jmin, int jmax, int* ierr);
    #void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) V, TYPE(mvec_ptr) Vblock, int jmin, int jmax, int* ierr);
    #void SUBR(mvec_set_block)(TYPE(mvec_ptr) V, TYPE(const_mvec_ptr) Vblock, int jmin, int jmax, int* ierr);
    _declare(None, _prefix+'mvec_to_mvec', (_mvec_ptr, _mvec_ptr, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'mvec_view_block', (_mvec_ptr, _mvec_ptr_p, c_int, c_int, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'mvec_get_block', (_mvec_ptr, _mvec_ptr, c_int, c_int, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'mvec_set_block', (_mvec_ptr, _mvec_ptr, c_int, c_int, c_int_p), skip_if_missing=True)

    #void SUBR(sdMat_view_block)(TYPE(sdMat_ptr) M, TYPE(sdMat_ptr)* Mblock, int imin, int imax, int jmin, int jmax, int* ierr);
    #void SUBR(sdMat_get_block)(TYPE(const_sdMat_ptr) M, TYPE(sdMat_ptr) Mblock, int imin, int imax, int jmin, int jmax, int* ierr);
    #void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) M, TYPE(const_sdMat_ptr) Mblock, int imin, int imax, int jmin, int jmax, int* ierr);
    _declare(None, _prefix+'sdMat_view_block', (_sdMat_ptr, _sdMat_ptr_p, c_int, c_int, c_int, c_int, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sdMat_get_block', (_sdMat_ptr, _sdMat_ptr, c_int, c_int, c_int, c_int, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sdMat_set_block', (_sdMat_ptr, _sdMat_ptr, c_int, c_int, c_int, c_int, c_int_p), skip_if_missing=True)

    #void SUBR(mvec_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* ierr);
    #void SUBR(sdMat_put_value)(TYPE(sdMat_ptr) V, _ST_ value, int* ierr);
    #void SUBR(mvec_random)(TYPE(mvec_ptr) V, int* ierr);
    #void SUBR(sdMat_random)(TYPE(sdMat_ptr) V, int* ierr);
    #void SUBR(mvec_print)(TYPE(const_mvec_ptr) V, int* ierr);
    #void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) M, int* ierr);
    _declare(None, _prefix+'mvec_put_value', (_mvec_ptr, _ST_, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sdMat_put_value', (_sdMat_ptr, _ST_, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'mvec_random', (_mvec_ptr, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sdMat_random', (_sdMat_ptr, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'mvec_print', (_mvec_ptr, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sdMat_print', (_sdMat_ptr, c_int_p), skip_if_missing=True)

    #void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) V, _MT_* vnrm, int *ierr);
    #void SUBR(mvec_normalize)(TYPE(mvec_ptr) V, _MT_* vnrm, int* ierr);
    #void SUBR(mvec_scale)(TYPE(mvec_ptr) V, _ST_ scalar, int* ierr);
    #void SUBR(mvec_vscale)(TYPE(mvec_ptr) V, const _ST_* scalar, int* ierr);
    #void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) X, _ST_ beta,  TYPE(mvec_ptr) Y, int* ierr);
    #void SUBR(mvec_vadd_mvec)(const _ST_ alpha[], TYPE(const_mvec_ptr) X, const _ST_ beta,  TYPE(mvec_ptr) Y, int* ierr);
    #void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) V, TYPE(const_mvec_ptr) W, _ST_* vw, int* ierr);
    #void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) V, TYPE(const_mvec_ptr) W, _ST_ beta, TYPE(sdMat_ptr) C, int* ierr);
    _declare(None, _prefix+'mvec_norm2', (_mvec_ptr, _MT_p, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'mvec_normalize', (_mvec_ptr, _MT_p, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'mvec_scale', (_mvec_ptr, _ST_, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'mvec_vscale', (_mvec_ptr, _ST_p, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'mvec_add_mvec', (_ST_, _mvec_ptr, _ST_, _mvec_ptr, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'mvec_vadd_mvec', (_ST_p, _mvec_ptr, _ST_, _mvec_ptr, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'mvec_dot_mvec', (_mvec_ptr, _mvec_ptr, _ST_p, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'mvecT_times_mvec', (_ST_, _mvec_ptr, _mvec_ptr, _ST_, _sdMat_ptr, c_int_p), skip_if_missing=True)

    #void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) V, TYPE(const_sdMat_ptr) C, _ST_ beta,  TYPE(mvec_ptr) W, int* ierr);
    #void SUBR(mvec_times_sdMat_inplace)(TYPE(mvec_ptr) V, TYPE(const_sdMat_ptr) M, int *ierr);
    #void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) A, _ST_ beta,  TYPE(sdMat_ptr) B, int* ierr);
    #void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, TYPE(const_sdMat_ptr) W, _ST_ beta, TYPE(sdMat_ptr) C, int* ierr);
    #void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, TYPE(const_sdMat_ptr) W, _ST_ beta, TYPE(sdMat_ptr) C, int* ierr);
    _declare(None, _prefix+'mvec_times_sdMat', (_ST_, _mvec_ptr, _sdMat_ptr, _ST_, _mvec_ptr, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'mvec_times_sdMat_inplce', (_mvec_ptr, _sdMat_ptr, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sdMat_add_sdMat', (_ST_, _sdMat_ptr, _ST_, _sdMat_ptr, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sdMat_times_sdMat', (_ST_, _sdMat_ptr, _ST_, _sdMat_ptr, _sdMat_ptr, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sdMatT_times_sdMat', (_ST_, _sdMat_ptr, _ST_, _sdMat_ptr, _sdMat_ptr, c_int_p), skip_if_missing=True)

    #void SUBR(sparseMat_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* ierr);
    #void SUBR(sparseMatT_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* ierr);
    #void SUBR(sparseMat_times_mvec_vadd_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, const _ST_ shifts[], TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* ierr);
    _declare(None, _prefix+'sparseMat_times_mvec', (_ST_, _sparseMat_ptr, _mvec_ptr, _ST_, _mvec_ptr, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sparseMatT_times_mvec', (_ST_, _sparseMat_ptr, _mvec_ptr, _ST_, _mvec_ptr, c_int_p), skip_if_missing=True)
    _declare(None, _prefix+'sparseMat_times_mvec_vadd_mvec', (_ST_, _sparseMat_ptr, _ST_p, _mvec_ptr, _ST_, _mvec_ptr, c_int_p), skip_if_missing=True)

    #void SUBR(mvec_QR)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, int* ierr);
    _declare(None, _prefix+'mvec_QR', (_mvec_ptr, _sdMat_ptr, c_int_p), skip_if_missing=True)

    # NOT IMPLEMENTED
    #void SUBR(sparseMat_create_fromRowFunc)(TYPE(sparseMat_ptr) *A, const_comm_ptr_t comm, gidx_t nrows, gidx_t ncols, lidx_t maxnne, int (*rowFunPtr)(ghost_gidx_t,ghost_lidx_t*,ghost_gidx_t*,void*), int *ierr);

#void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, Dmvec_t* reV, Dmvec_t* imV, int *ierr);
#void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, Smvec_t* reV, Smvec_t* imV, int *ierr);
_declare(None, 'phist_Cmvec_split', (Cmvec_ptr, Smvec_ptr, Smvec_ptr, c_int_p), skip_if_missing=True)
_declare(None, 'phist_Zmvec_split', (Zmvec_ptr, Dmvec_ptr, Dmvec_ptr, c_int_p), skip_if_missing=True)





#--------------------------------------------------------------------------------
# tests / sample code
if __name__ == '__main__':

    # also import phist_tools
    from phist_tools import *

    # initialize
    PYST_CHK_IERR(phist_kernels_init)

    # create comm
    myComm = comm_ptr()
    PYST_CHK_IERR(phist_comm_create, myComm)

    # check comm size and rank
    nprocs = c_int()
    me = c_int()
    PYST_CHK_IERR(phist_comm_get_size, myComm, nprocs)
    PYST_CHK_IERR(phist_comm_get_rank, myComm, me)
    print 'nprocs', nprocs, 'me', me

    # create map
    nglobal = gidx(10)
    myMap = map_ptr()
    PYST_CHK_IERR(phist_map_create, myMap, myComm, nglobal)

    # get number of local elements
    nlocal = lidx()
    PYST_CHK_IERR(phist_map_get_local_length, myMap, nlocal)
    print 'global elements', nglobal, 'local elements', nlocal

    # delete map
    PYST_CHK_IERR(phist_map_delete, myMap)

    # create a matrix
    A = DsparseMat_ptr()
    PYST_CHK_IERR(phist_Dcreate_matrix, A, myComm, 'spinSZ8')
    # get map
    PYST_CHK_IERR(phist_DsparseMat_get_domain_map, A, myMap)

    # create two mvecs
    vec1 = Dmvec_ptr()
    vec2 = Dmvec_ptr()
    vec3 = Dmvec_ptr()
    PYST_CHK_IERR(phist_Dmvec_create, vec1, myMap, 4)
    PYST_CHK_IERR(phist_Dmvec_create, vec2, myMap, 2)
    PYST_CHK_IERR(phist_Dmvec_view_block, vec1, vec3, 0,1 )

    # create some sdMats
    mat1 = DsdMat_ptr()
    mat2 = DsdMat_ptr()
    PYST_CHK_IERR(phist_DsdMat_create, mat1, 4, 2, myComm)
    PYST_CHK_IERR(phist_DsdMat_view_block, mat1, mat2, 0,1, 0,1)

    # put random data in all mvecs
    PYST_CHK_IERR(phist_Dmvec_random, vec1)
    PYST_CHK_IERR(phist_Dmvec_random, vec2)

    # QR-orthogonalize vec1
    PYST_CHK_IERR(phist_Dmvec_QR, vec2, mat2)

    # calculate vec1^T vec2
    PYST_CHK_IERR(phist_DmvecT_times_mvec, 1, vec1, vec2, 0, mat1)

    # calculate A*vec2
    PYST_CHK_IERR(phist_DsparseMat_times_mvec, 1, A, vec2, 0, vec3)

    # delete sdMats
    PYST_CHK_IERR(phist_DsdMat_delete, mat2)
    PYST_CHK_IERR(phist_DsdMat_delete, mat1)

    # delete mvecs
    PYST_CHK_IERR(phist_Dmvec_delete, vec3)
    PYST_CHK_IERR(phist_Dmvec_delete, vec2)
    PYST_CHK_IERR(phist_Dmvec_delete, vec1)

    # delete matrix
    PYST_CHK_IERR(phist_DsparseMat_delete, A)

    # delete comm
    PYST_CHK_IERR(phist_comm_delete, myComm)

    # finalize
    PYST_CHK_IERR(phist_kernels_finalize)

