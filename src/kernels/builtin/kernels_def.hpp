/*! \file builtin/kernels_def.hpp
 * included by builtin/kernels.cpp
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 * \author "Jonas Thies <Jonas.Thies@DLR.de>
 *
*/

#ifndef PHIST_HAVE_MPI  
#error "builtin kernels don't work without MPI"
#endif

// Declaration of Fortran implemented functions
extern "C" {
  void SUBR(crsMat_create_fromRowFunc_f)(TYPE(sparseMat_ptr)*,const_comm_ptr_t comm,gidx_t,gidx_t, 
      lidx_t, void (*)(ghost_gidx_t,ghost_lidx_t*,ghost_gidx_t*,void*), int*);
  void SUBR(crsMat_delete_f)(TYPE(sparseMat_ptr) A, int* iflag);
  void SUBR(crsMat_get_map_f)(TYPE(const_sparseMat_ptr),const_map_ptr_t*,int*);
  void SUBR(crsMat_read_mm_f)(void*A,const_comm_ptr_t comm, int fname_len, const char* fname, int* iflag);
  void SUBR(mvecT_times_mvec_f)(_ST_,TYPE(const_mvec_ptr),TYPE(const_mvec_ptr),_ST_,TYPE(sdMat_ptr),int*);
  void SUBR(mvec_QR_f)(TYPE(mvec_ptr),TYPE(sdMat_ptr),int*);
  void SUBR(mvec_add_mvec_f)(_ST_,TYPE(const_mvec_ptr),_ST_,TYPE(mvec_ptr),int*);
  void SUBR(mvec_create_f)(TYPE(mvec_ptr)*,const_map_ptr_t,lidx_t,int*);
  void SUBR(mvec_delete_f)(TYPE(mvec_ptr),int*);
  void SUBR(mvec_dot_mvec_f)(TYPE(const_mvec_ptr),TYPE(const_mvec_ptr),_ST_*,int*);
  void SUBR(mvec_extract_view_f)(TYPE(mvec_ptr),_ST_**,lidx_t*,int*);
  void SUBR(mvec_gather_mvecs_f)(TYPE(mvec_ptr),TYPE(const_mvec_ptr) W[], int, int*);
  void SUBR(mvec_get_block_f)(TYPE(const_mvec_ptr),TYPE(mvec_ptr),int,int,int*);
  void SUBR(mvec_get_map_f)(TYPE(const_mvec_ptr),const_map_ptr_t*,int*);
  void SUBR(mvec_my_length_f)(TYPE(const_mvec_ptr),lidx_t*,int*);
  void SUBR(mvec_norm2_f)(TYPE(const_mvec_ptr),_MT_*,int*);
  void SUBR(mvec_num_vectors_f)(TYPE(const_mvec_ptr),int*,int*);
  void SUBR(mvec_print_f)(TYPE(const_mvec_ptr),int*);
  void SUBR(mvec_put_value_f)(TYPE(mvec_ptr),_ST_,int*);
  void SUBR(mvec_put_func_f)(TYPE(mvec_ptr),void(*)(ghost_gidx_t,ghost_lidx_t,void*),int*);
  void SUBR(mvec_random_f)(TYPE(mvec_ptr),int*);
  void SUBR(mvec_scale_f)(TYPE(mvec_ptr),_ST_,int*);
  void SUBR(mvec_scatter_mvecs_f)(TYPE(const_mvec_ptr),TYPE(mvec_ptr) W[], int, int*);
  void SUBR(mvec_set_block_f)(TYPE(mvec_ptr),TYPE(const_mvec_ptr),int,int,int*);
  void SUBR(mvec_times_sdMat_f)(_ST_,TYPE(const_mvec_ptr),TYPE(const_sdMat_ptr),_ST_,TYPE(mvec_ptr),int*);
  void SUBR(mvec_times_sdMat_inplace_f)(TYPE(mvec_ptr) V, TYPE(const_sdMat_ptr) M, int*);
  void SUBR(mvec_to_mvec_f)(TYPE(const_mvec_ptr),TYPE(mvec_ptr),int*);
  void SUBR(mvec_vadd_mvec_f)(const _ST_*,TYPE(const_mvec_ptr),_ST_,TYPE(mvec_ptr),int*);
  void SUBR(mvec_view_block_f)(TYPE(mvec_ptr),TYPE(mvec_ptr)*,int,int,int*);
  void SUBR(mvec_vscale_f)(TYPE(mvec_ptr),const _ST_*,int*);
  void SUBR(sdMatT_times_sdMat_f)(_ST_,TYPE(const_sdMat_ptr),TYPE(const_sdMat_ptr),_ST_,TYPE(sdMat_ptr),int*);
  void SUBR(sdMat_times_sdMatT_f)(_ST_,TYPE(const_sdMat_ptr),TYPE(const_sdMat_ptr),_ST_,TYPE(sdMat_ptr),int*);
#ifdef PHIST_HIGH_PRECISION_KERNELS
  void SUBR(sdMat_cholesky_f)(TYPE(sdMat_ptr),int*,int*,int*);
  void SUBR(sdMat_backwardSubst_sdMat_f)(const TYPE(sdMat_ptr), const int*, int, TYPE(sdMat_ptr), int*);
  void SUBR(sdMat_forwardSubst_sdMat_f)(const TYPE(sdMat_ptr), const int*, int, TYPE(sdMat_ptr), int*);
#endif
  void SUBR(sdMat_add_sdMat_f)(_ST_,TYPE(const_sdMat_ptr),_ST_,TYPE(sdMat_ptr),int*);
  void SUBR(sdMat_create_f)(TYPE(sdMat_ptr)*,int,int,const_comm_ptr_t,int*);
  void SUBR(sdMat_create_view_f)(TYPE(sdMat_ptr)*,const_comm_ptr_t,void*, lidx_t, int,int,int*);
  void SUBR(sdMat_delete_f)(TYPE(sdMat_ptr),int*);
  void SUBR(sdMat_extract_view_f)(TYPE(sdMat_ptr),_ST_**,lidx_t*,int*);
  void SUBR(sdMat_get_block_f)(TYPE(const_mvec_ptr),TYPE(mvec_ptr),int,int,int,int,int*);
  void SUBR(sdMat_get_ncols_f)(TYPE(const_sdMat_ptr),int*,int*);
  void SUBR(sdMat_get_nrows_f)(TYPE(const_sdMat_ptr),int*,int*);
  void SUBR(sdMat_print_f)(TYPE(const_sdMat_ptr),int*);
  void SUBR(sdMat_put_value_f)(TYPE(sdMat_ptr),_ST_,int*);
  void SUBR(sdMat_random_f)(TYPE(sdMat_ptr),int*);
  void SUBR(sdMat_identity_f)(TYPE(sdMat_ptr),int*);
  void SUBR(sdMat_set_block_f)(TYPE(mvec_ptr),TYPE(const_mvec_ptr),int,int,int,int,int*);
  void SUBR(sdMat_times_sdMat_f)(_ST_,TYPE(const_sdMat_ptr),TYPE(const_sdMat_ptr),_ST_,TYPE(sdMat_ptr),int*);
  void SUBR(sdMat_view_block_f)(TYPE(mvec_ptr),TYPE(mvec_ptr)*,int,int,int,int,int*);
}

extern "C" void SUBR(type_avail)(int *iflag)
{
  *iflag=0;
}

extern "C" void SUBR(sparseMat_read_mm)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t vcomm,
        const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  if (filename==NULL)
  {
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
  PHIST_CHK_IERR(SUBR(crsMat_read_mm_f)(A,vcomm,strlen(filename),filename,iflag),*iflag);
}

extern "C" void SUBR(sparseMat_read_bin)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t vcomm,
const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_hb)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t vcomm,
const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_get_row_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR( SUBR(crsMat_get_map_f) (A,map,iflag), *iflag);
}

extern "C" void SUBR(sparseMat_get_col_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR( SUBR(sparseMat_get_row_map) (A,map,iflag), *iflag);
}

extern "C" void SUBR(sparseMat_get_domain_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR( SUBR(sparseMat_get_row_map) (A,map,iflag), *iflag);
}

extern "C" void SUBR(sparseMat_get_range_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR( SUBR(sparseMat_get_row_map) (A,map,iflag), *iflag);
}

extern "C" void SUBR(mvec_create)(TYPE(mvec_ptr)* V, 
    const_map_ptr_t map, lidx_t nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR( SUBR(mvec_create_f) (V,map,nvec,iflag), *iflag);
}

extern "C" void SUBR(mvec_create_view)(TYPE(mvec_ptr)* V, const_map_ptr_t map, 
    _ST_* values, lidx_t lda, int nvec,
    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_create)(TYPE(sdMat_ptr)* M, 
    int nrows, int ncols, const_comm_ptr_t comm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_create_f)(M,nrows,ncols,comm,iflag),*iflag);
}

extern "C" void SUBR(sdMat_create_view)(TYPE(sdMat_ptr)* M, const_comm_ptr_t comm,
        _ST_* values, lidx_t lda, int nrows, int ncols,
        int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_create_view_f)(M,comm,(void*)values,lda,nrows,ncols,iflag),*iflag);
}
                  

extern "C" void SUBR(mvec_my_length)(TYPE(const_mvec_ptr) V, lidx_t* len, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_my_length_f)(V,len,iflag),*iflag);
}

extern "C" void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) V, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_get_map_f)(V,map,iflag),*iflag);
}

extern "C" void SUBR(mvec_get_comm)(TYPE(const_mvec_ptr) V, const_comm_ptr_t* comm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  const_map_ptr_t map;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(V,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_comm(map,comm,iflag),*iflag);
}

extern "C" void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) V, int* nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors_f)(V,nvec,iflag),*iflag);
}

extern "C" void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) M, int* nrows, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows_f)(M,nrows,iflag),*iflag);
}

extern "C" void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) M, int* ncols, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols_f)(M,ncols,iflag),*iflag);
}

extern "C" void SUBR(mvec_extract_view)(TYPE(mvec_ptr) V, _ST_** val, lidx_t* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_extract_view_f)(V,val,lda,iflag),*iflag);
}

extern "C" void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) V, _ST_** val, lidx_t* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view_f)(V,val,lda,iflag),*iflag);
}

extern "C" void SUBR(mvec_to_mvec)(TYPE(const_mvec_ptr) v_in, TYPE(mvec_ptr) v_out, int* iflag)
{
#ifndef PHIST_HAVE_PARMETIS
  // we don't need this when no reordering occurs!
  *iflag=PHIST_NOT_IMPLEMENTED;
#else
  PHIST_CHK_IERR(SUBR(mvec_to_mvec_f)(v_in, v_out, iflag),*iflag);
#endif
}

extern "C" void SUBR(mvec_view_block)(TYPE(mvec_ptr) V,
    TYPE(mvec_ptr)* Vblock,
    int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_view_block_f)(V,Vblock,jmin,jmax,iflag),*iflag);
}

extern "C" void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) V,
    TYPE(mvec_ptr) Vblock,
    int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_get_block_f)(V,Vblock,jmin,jmax,iflag),*iflag);
}

extern "C" void SUBR(mvec_set_block)(TYPE(mvec_ptr) V,
    TYPE(const_mvec_ptr) Vblock,
    int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_set_block_f)(V,Vblock,jmin,jmax,iflag),*iflag);
}

extern "C" void SUBR(sdMat_view_block)(TYPE(mvec_ptr) M, 
    TYPE(mvec_ptr)* Mblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_view_block_f)(M,Mblock,imin,imax,jmin,jmax,iflag),*iflag);
}

extern "C" void SUBR(sdMat_get_block)(TYPE(const_mvec_ptr) M, 
    TYPE(mvec_ptr) Mblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_get_block_f)(M,Mblock,imin,imax,jmin,jmax,iflag),*iflag);
}

extern "C" void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) M, 
    TYPE(const_sdMat_ptr) Mblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_set_block_f)(M,Mblock,imin,imax,jmin,jmax,iflag),*iflag);
}

extern "C" void SUBR(sparseMat_delete)(TYPE(sparseMat_ptr) A, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(crsMat_delete_f)(A,iflag), *iflag);
}

extern "C" void SUBR(mvec_delete)(TYPE(mvec_ptr) V, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR( SUBR(mvec_delete_f) (V,iflag),*iflag);
}

extern "C" void SUBR(sdMat_delete)(TYPE(sdMat_ptr) M, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_delete_f)(M,iflag),*iflag);
}

extern "C" void SUBR(mvec_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#ifdef PHIST_PERFCHECK
  lidx_t nlocal;
  int nV;
  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&nlocal,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag);
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,0,nlocal, STREAM_STORE(nV*nlocal*sizeof(_ST_)));
#endif
  PHIST_CHK_IERR(SUBR(mvec_put_value_f)(V,value,iflag),*iflag);
}

extern "C" void SUBR(mvec_put_func)(TYPE(mvec_ptr) V,
        int (*funPtr)(ghost_gidx_t,ghost_lidx_t,void*), int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_put_func_f)(V,(void(*)(ghost_gidx_t,ghost_lidx_t,void*))funPtr,iflag),*iflag);
}

extern "C" void SUBR(sdMat_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_put_value_f)(V,value,iflag),*iflag);
}

extern "C" void SUBR(mvec_random)(TYPE(mvec_ptr) V, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_random_f)(V,iflag),*iflag);
}

extern "C" void SUBR(mvec_print)(TYPE(const_mvec_ptr) V, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_print_f)(V,iflag),*iflag);
}

extern "C" void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) M, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_print_f)(M,iflag),*iflag);
}

extern "C" void SUBR(sdMat_random)(TYPE(sdMat_ptr) M, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_random_f)(M,iflag),*iflag);
}

extern "C" void SUBR(sdMat_identity)(TYPE(sdMat_ptr) M, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_identity_f)(M,iflag),*iflag);
}

extern "C" void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) V,
    _MT_* vnrm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_norm2_f)(V,vnrm,iflag),*iflag);
}

extern "C" void SUBR(mvec_normalize)(TYPE(mvec_ptr) V,
    _MT_* vnrm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"  
  PHIST_CHK_IERR(SUBR(mvec_norm2_f)(V,vnrm,iflag),*iflag);
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nvec,iflag),*iflag);
  _ST_*scale = new _ST_[nvec];
  for(int i = 0; i < nvec; i++)
    scale[i] = st::one()/vnrm[i];
  PHIST_CHK_IERR(SUBR(mvec_vscale_f)(V,scale,iflag),*iflag);
  delete[] scale;
}

extern "C" void SUBR(mvec_scale)(TYPE(mvec_ptr) V, 
    _ST_ scalar, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_scale_f)(V,scalar,iflag),*iflag);
}

extern "C" void SUBR(mvec_vscale)(TYPE(mvec_ptr) V, 
    const _ST_* scalar, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_vscale_f)(V,scalar,iflag),*iflag);
}

extern "C" void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) X,
    _ST_ beta,  TYPE(mvec_ptr)       Y, 
    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_add_mvec_f)(alpha,X,beta,Y,iflag),*iflag);
}

extern "C" void SUBR(mvec_vadd_mvec)(const _ST_ alpha[], TYPE(const_mvec_ptr) X,
    _ST_ beta,  TYPE(mvec_ptr)       Y, 
    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_vadd_mvec_f)(alpha,X,beta,Y,iflag),*iflag);
}

extern "C" void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) A,
    _ST_ beta,  TYPE(sdMat_ptr)       B, 
    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_add_sdMat_f)(alpha,A,beta,B,iflag),*iflag);
}

extern "C" void SUBR(sparseMat_times_mvec_communicate)(TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag = 0;
}

extern "C" void SUBR(sparseMat_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
    TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

#ifdef PHIST_TIMEMONITOR
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(x, &nvec, iflag), *iflag);
  for(int i = 0; i < nvec; i++)
    phist_totalMatVecCount();
#endif

  void SUBR(crsMat_times_mvec_f)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
      TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag);
  PHIST_CHK_IERR(SUBR(crsMat_times_mvec_f)(alpha,A,x,beta,y,iflag),*iflag);
}

extern "C" void SUBR(sparseMatT_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
    TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

//! y[i]=alpha*(A*x[i]+shifts[i]*x[i]) + beta*y[i]
extern "C" void SUBR(sparseMat_times_mvec_vadd_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A,
        const _ST_ shifts[], TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

#ifdef PHIST_TIMEMONITOR
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(x, &nvec, iflag), *iflag);
  for(int i = 0; i < nvec; i++)
    phist_totalMatVecCount();
#endif

  void SUBR(crsMat_times_mvec_vadd_mvec_f)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
      const _ST_ shifts[], TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag);
  PHIST_CHK_IERR(SUBR(crsMat_times_mvec_vadd_mvec_f)(alpha,A,shifts,x,beta,y,iflag),*iflag);
}

extern "C" void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) v, 
    TYPE(const_mvec_ptr) w, 
    _ST_* s, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#ifdef PHIST_PERFCHECK
  lidx_t nlocal;
  int nV;
  PHIST_CHK_IERR(SUBR(mvec_my_length)(v,&nlocal,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(v,&nV,iflag),*iflag);
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,0,nlocal, STREAM_LOAD(2*nV*nlocal*sizeof(_ST_)));
#endif
  PHIST_CHK_IERR(SUBR(mvec_dot_mvec_f)(v,w,s,iflag),*iflag);
}

extern "C" void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
    TYPE(const_sdMat_ptr) C, 
    _ST_ beta, TYPE(mvec_ptr) W, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#ifdef PHIST_PERFCHECK
  lidx_t nlocal;
  int nV, nW;
  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&nlocal,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&nW,iflag),*iflag);
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,nW,nlocal, STREAM_TRIAD((nV+2*nW)*nlocal*sizeof(_ST_)));
#endif
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat_f)(alpha,V,C,beta,W,iflag),*iflag);
}

extern "C" void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
    TYPE(const_sdMat_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_times_sdMat_f)(alpha,V,W,beta,C,iflag),*iflag);
}

extern "C" void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
    TYPE(const_sdMat_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMatT_times_sdMat_f)(alpha,V,W,beta,C,iflag),*iflag);
}

extern "C" void SUBR(sdMat_times_sdMatT)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
    TYPE(const_sdMat_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_times_sdMatT_f)(alpha,V,W,beta,C,iflag),*iflag);
}

extern "C" void SUBR(sdMat_cholesky)(TYPE(sdMat_ptr) C, int* perm, int* rank, int* iflag)
{
#ifdef PHIST_HIGH_PRECISION_KERNELS
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_cholesky_f)(C,perm,rank,iflag),*iflag);
#else
  *iflag=PHIST_NOT_IMPLEMENTED;
#endif
}

extern "C" void SUBR(sdMat_backwardSubst_sdMat)(const TYPE(sdMat_ptr) R, int* perm, int rank, TYPE(sdMat_ptr) X, int* iflag)
{
#ifdef PHIST_HIGH_PRECISION_KERNELS
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_backwardSubst_sdMat_f)(R,perm,rank,X,iflag),*iflag);
#else
  *iflag=PHIST_NOT_IMPLEMENTED;
#endif
}

extern "C" void SUBR(sdMat_forwardSubst_sdMat)(const TYPE(sdMat_ptr) R, int* perm, int rank, TYPE(sdMat_ptr) X, int* iflag)
{
#ifdef PHIST_HIGH_PRECISION_KERNELS
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_forwardSubst_sdMat_f)(R,perm,rank,X,iflag),*iflag);
#else
  *iflag=PHIST_NOT_IMPLEMENTED;
#endif
}


extern "C" void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
    TYPE(const_mvec_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#ifdef PHIST_PERFCHECK
  lidx_t nlocal;
  int nV, nW;
  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&nlocal,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nV,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&nW,iflag),*iflag);
  PHIST_PERFCHECK_VERIFY(__FUNCTION__,nV,nW,nlocal, STREAM_LOAD((nV+nW)*nlocal*sizeof(_ST_)));
#endif
  PHIST_CHK_IERR(SUBR(mvecT_times_mvec_f)(alpha,V,W,beta,C,iflag),*iflag);
}

extern "C" void SUBR(mvec_QR)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  // trilinos tsqr would be better but seems to expect column-major mvecs
  PHIST_CHK_NEG_IERR(SUBR(mvec_QR_f)(V,R,iflag),*iflag);
}

extern "C" void SUBR(mvec_gather_mvecs)(TYPE(mvec_ptr) V, TYPE(const_mvec_ptr) W[], int nblocks, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_MARK_AS_EXPERIMENTAL(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_gather_mvecs_f)(V,W,nblocks,iflag),*iflag);
}

extern "C" void SUBR(mvec_scatter_mvecs)(TYPE(const_mvec_ptr) V, TYPE(mvec_ptr) W[], int nblocks, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_MARK_AS_EXPERIMENTAL(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_scatter_mvecs_f)(V,W,nblocks,iflag),*iflag);
}

extern "C" void SUBR(mvec_times_sdMat_inplace)(TYPE(mvec_ptr) V, TYPE(const_sdMat_ptr) M, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace_f)(V, M, iflag), *iflag);
}

extern "C" void SUBR(sparseMat_create_fromRowFunc)(TYPE(sparseMat_ptr) *A, const_comm_ptr_t vcomm,
        gidx_t nrows, gidx_t ncols, lidx_t maxnne, 
                int (*rowFunPtr)(ghost_gidx_t,ghost_lidx_t*,ghost_gidx_t*,void*), 
                int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(crsMat_create_fromRowFunc_f)(A, vcomm, nrows, ncols, maxnne, 
        (void(*)(ghost_gidx_t,ghost_lidx_t*,ghost_gidx_t*,void*))rowFunPtr, iflag), *iflag);
}

#include "../kernels_nogpu.c"


