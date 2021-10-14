/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/

// given a file base like xyz, tries to find a matrix file [D,S,Z,C]xyz.[mm,bin,...]
// that the kernel library can read and reads the sparse matrix.
void SUBR(read_mat)(const char* filebase,phist_const_comm_ptr comm,int nglob, int mglob, 
        TYPE(sparseMat_ptr) *ptr, int* iflag);

//! \name some functions for initializing sparseMats and mvecs
//!@{

  //! create identity matrix. 
  
  //! If called with row=-1, cols[0] and cols[1] are interpreted as global number of
  //! rows and columns, respectively (gnrows and gncols).
  //! If called with row=-2, the function is "un-initialized" (the global sizes are reset).
  //! If gnrows>gncols, creates [speye(gncols);zeros(gnrows-gncols,gncols)].
  //! If gncols>gnrows, creates [speye(gnrows),zeros(gnrows,gncols-gnrows)].
  //! If arg!=NULL, the entries are scaled by *((_ST_*)arg)
  int PHIST_TG_PREFIX(idfunc)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  //! some row function that uses a workspace. Will only create the identity matrix if
  //! initialized using idfunc_init_workspace is called. Requires arg to be a struct with
  //! a pointer to idfunc_with_workspace_arg.
  int PHIST_TG_PREFIX(idfunc_with_workspace)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  //! 'constructor'
  int PHIST_TG_PREFIX(idfunc_init_workspace)(void *arg, void** work);
  
  //! argument for idfunc_with_workspace
  typedef struct
  {
    phist_gidx gnrows;
    phist_gidx gncols;
    _ST_ scale;
  } PHIST_TG_PREFIX(idfunc_with_workspace_arg);

  typedef struct
  {
    _ST_* data;
    PHIST_TG_PREFIX(idfunc_with_workspace_arg)* arg;
  } PHIST_TG_PREFIX(idfunc_workspace);

  //! creates a simple tridiagonal Hermitian and positive definite matrix. Before passing this function to
  //! sparseMat_create_fromRowFunc, it must be initialized by calling it with row=-1 and cols[0] containing
  //! the global number of rows (gnrows). The maximum row length is 3.
  //! If called with row=-2, the function is "un-initialized" (the global sizes are reset).
  int PHIST_TG_PREFIX(hpd_tridiag)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  // tridiagonal matrix with entries [-1 2 -1] or [-i 2 -i] in the complex case (hpd)
  int PHIST_TG_PREFIX(lapl_tridiag)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  //! creates a simple tridiagonal non-Hermitian but positive definite matrix. For usage info, see hpd_tridiag.
  int PHIST_TG_PREFIX(nhpd_tridiag)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  //! creates a simple tridiagonal Hermitian and indefinite matrix. For usage info, see hpd_tridiag.
  int PHIST_TG_PREFIX(hid_tridiag)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  //! creates a simple tridiagonal non-Hermitian and indefinite matrix. For usage info, see hpd_tridiag.
  int PHIST_TG_PREFIX(nhid_tridiag)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  //! creates an approximate inverse of hpd_tridiag (the inverse of the 2x2 block diagonal approximation of A)
  int PHIST_TG_PREFIX(hpd_tridiag_ainv)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  //! creates an approximate inverse of lapl_tridiag (the inverse of the 2x2 block diagonal approximation of A)
  int PHIST_TG_PREFIX(lapl_tridiag_ainv)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  //! creates an approximate inverse of nhpd_tridiag (the inverse of the 2x2 block diagonal approximation of A)
  int PHIST_TG_PREFIX(nhpd_tridiag_ainv)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  //! creates an approximate inverse of hid_tridiag (the inverse of the 2x2 block diagonal approximation of A)
  int PHIST_TG_PREFIX(hid_tridiag_ainv)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  //! creates an approximate inverse of nhid_tridiag (the inverse of the 2x2 block diagonal approximation of A)
  int PHIST_TG_PREFIX(nhid_tridiag_ainv)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  //! defines a matrix with only a constant subdiagonal and an entry (1,N)
  //! that defines a periodic "right shift" of vector elements.
  //! To initialize and finalize, call with row=-1 and -2, resp. (see hpd_tridiag above)
  int PHIST_TG_PREFIX(right_shift_perio)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  //! defines a matrix with only a constant subdiagonal that defines a non-periodic "right shift" of vector elements.
  //! To initialize and finalize, call with row=-1 and -2, resp. (see hpd_tridiag above)
  int PHIST_TG_PREFIX(right_shift)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  //! defines a matrix with only a constant superdiagonal and an entry (N,1)
  //! that defines a periodic "left shift" of vector elements.
  //! To initialize and finalize, call with row=-1 and -2, resp. (see hpd_tridiag above)
  int PHIST_TG_PREFIX(left_shift_perio)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  //! defines a matrix with only a constant superdiagonal that defines a non-periodic "left shift" of vector elements.
  //! To initialize and finalize, call with row=-1 and -2, resp. (see hpd_tridiag above)
  int PHIST_TG_PREFIX(left_shift)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);

  //! fill an mvec with characteristic values [1 2 3 ...] so we can more easily check permutation functionality.
  //! last_arg must point to an array of two ints containing the number of rows and columns in the target mvec.
  int PHIST_TG_PREFIX(mvec123func)(ghost_gidx i, ghost_lidx j, void* val, void* last_arg);

  //! reverse vector of mvec123func
  //! last_arg must point to an array of two ints containing the number of rows and columns in the target mvec.
  int PHIST_TG_PREFIX(mvec321func)(ghost_gidx i, ghost_lidx j, void* val, void* last_arg);

