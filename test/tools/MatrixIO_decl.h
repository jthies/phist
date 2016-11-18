
// given a file base like xyz, tries to find a matrix file [D,S,Z,C]xyz.[mm,bin,...]
// that the kernel library can read and reads the sparse matrix.
void SUBR(read_mat)(const char* filebase,phist_const_comm_ptr comm,int nglob, int mglob, 
        TYPE(sparseMat_ptr) *ptr, int* iflag);

//! \name some functions for initializing sparseMats and mvecs
//!@{

  //! create identity matrix. If called with row=-1, cols[0] and cols[1] are interpreted as global number of
  //! rows and columns, respectively (gnrows and gncols).
  //! If gnrows>gncols, creates [speye(gncols);zeros(gnrows-gncols,gncols)].
  //! If gncols>gnrows, creates [speye(gnrows),zeros(gnrows,gncols-gnrows)].
  int PHIST_TG_PREFIX(idfunc)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  //! create some sparse matrix without specific properties
  int PHIST_TG_PREFIX(some_rowFunc)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);
  
  //! creates a simple tridiagonal Hermitian and positive definite matrix. Before passing this function to
  //! sparseMat_create_fromRowFunc, it must be initialized by calling it with row=-1 and cols[0] containing
  //! the global number of rows (gnrows). The maximum row length is 3.
  int PHIST_TG_PREFIX(hpd_tridiag)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  //! fill an mvec with characteristic values [1 2 3 ...] so we can more easily check permutation functionality
  int PHIST_TG_PREFIX(mvec123func)(ghost_gidx i, ghost_lidx j, void* val, void* last_arg);

  //! reverse vector of mvec123func
  int PHIST_TG_PREFIX(mvec321func)(ghost_gidx i, ghost_lidx j, void* val, void* last_arg);

