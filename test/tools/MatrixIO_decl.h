
// given a file base like xyz, tries to find a matrix file [D,S,Z,C]xyz.[mm,bin,...]
// that the kernel library can read and reads the sparse matrix.
void SUBR(read_mat)(const char* filebase,phist_const_comm_ptr comm,int nglob,TYPE(sparseMat_ptr) *ptr, int* iflag);

//! \name some functions for initializing sparseMats and mvecs
//!@{

  //!
  int PHIST_TG_PREFIX(idfunc)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  //!
  int PHIST_TG_PREFIX(some_rowFunc)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  //!
  int PHIST_TG_PREFIX(mvec123func)(ghost_gidx i, ghost_lidx j, void* val, void* last_arg);

  //!
  int PHIST_TG_PREFIX(mvec321func)(ghost_gidx i, ghost_lidx j, void* val, void* last_arg);

