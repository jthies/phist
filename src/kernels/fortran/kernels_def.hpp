extern "C" {

void SUBR(type_avail)(int *ierr)
{
  *ierr=0;
}

void SUBR(crsMat_read_mm)(TYPE(crsMat_ptr)* A, const char* filename,int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(crsMat_read_mm_f)(void*A,int fname_len, const char* fname, int* ierr);
  PHIST_CHK_IERR(SUBR(crsMat_read_mm_f)(A,strlen(filename),filename,ierr),*ierr);
}

void SUBR(crsMat_read_bin)(TYPE(crsMat_ptr)* A, const char* filename,int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=-99;
}

void SUBR(crsMat_read_hb)(TYPE(crsMat_ptr)* A, const char* filename,int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=-99;
}

void SUBR(crsMat_get_row_map)(TYPE(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(crsMat_get_map_f)(TYPE(const_crsMat_ptr),const_map_ptr_t*,int*);
  PHIST_CHK_IERR( SUBR(crsMat_get_map_f) (A,map,ierr), *ierr);
}

void SUBR(crsMat_get_col_map)(TYPE(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  PHIST_CHK_IERR( SUBR(crsMat_get_row_map) (A,map,ierr), *ierr);
}

void SUBR(crsMat_get_domain_map)(TYPE(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  PHIST_CHK_IERR( SUBR(crsMat_get_row_map) (A,map,ierr), *ierr);
}

void SUBR(crsMat_get_range_map)(TYPE(const_crsMat_ptr) A, const_map_ptr_t* map, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  PHIST_CHK_IERR( SUBR(crsMat_get_row_map) (A,map,ierr), *ierr);
}

void SUBR(mvec_create)(TYPE(mvec_ptr)* V, 
    const_map_ptr_t map, lidx_t nvec, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"  
  void SUBR(mvec_create_f)(TYPE(mvec_ptr)*,const_map_ptr_t,lidx_t,int*);
  PHIST_CHK_IERR( SUBR(mvec_create_f) (V,map,nvec,ierr), *ierr);
  // initialize with zero...
  PHIST_CHK_IERR( SUBR(mvec_put_value)(*V,st::zero(),ierr),*ierr);
}

void SUBR(mvec_create_view)(TYPE(mvec_ptr)* V, const_map_ptr_t map, 
    _ST_* values, lidx_t lda, int nvec,
    int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=-99;
}

void SUBR(sdMat_create)(TYPE(sdMat_ptr)* M, 
    int nrows, int ncols, const_comm_ptr_t comm, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"  
  void SUBR(sdMat_create_f)(TYPE(sdMat_ptr)*,int,int,const_comm_ptr_t,int*);
  PHIST_CHK_IERR(SUBR(sdMat_create_f)(M,nrows,ncols,comm,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(sdMat_put_value)(*M,st::zero(),ierr),*ierr);
}

void SUBR(mvec_my_length)(TYPE(const_mvec_ptr) V, lidx_t* len, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(mvec_my_length_f)(TYPE(const_mvec_ptr),lidx_t*,int*);
  PHIST_CHK_IERR(SUBR(mvec_my_length_f)(V,len,ierr),*ierr);
}

void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) V, const_map_ptr_t* map, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(mvec_get_map_f)(TYPE(const_mvec_ptr),const_map_ptr_t*,int*);
  PHIST_CHK_IERR(SUBR(mvec_get_map_f)(V,map,ierr),*ierr);
}

void SUBR(mvec_get_comm)(TYPE(const_mvec_ptr) V, const_comm_ptr_t* comm, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  const_map_ptr_t map;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(V,&map,ierr),*ierr);
  PHIST_CHK_IERR(phist_map_get_comm(map,comm,ierr),*ierr);
}

void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) V, int* nvec, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(mvec_num_vectors_f)(TYPE(const_mvec_ptr),int*,int*);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors_f)(V,nvec,ierr),*ierr);
}

void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) M, int* nrows, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(sdMat_get_nrows_f)(TYPE(const_sdMat_ptr),int*,int*);
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows_f)(M,nrows,ierr),*ierr);
}

void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) M, int* ncols, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(sdMat_get_ncols_f)(TYPE(const_sdMat_ptr),int*,int*);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols_f)(M,ncols,ierr),*ierr);
}

void SUBR(mvec_extract_view)(TYPE(mvec_ptr) V, _ST_** val, lidx_t* lda, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(mvec_extract_view_f)(TYPE(mvec_ptr),_ST_**,lidx_t*,int*);
  PHIST_CHK_IERR(SUBR(mvec_extract_view_f)(V,val,lda,ierr),*ierr);
}

void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) V, _ST_** val, lidx_t* lda, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(sdMat_extract_view_f)(TYPE(sdMat_ptr),_ST_**,lidx_t*,int*);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view_f)(V,val,lda,ierr),*ierr);
}

void SUBR(mvec_view_block)(TYPE(mvec_ptr) V,
    TYPE(mvec_ptr)* Vblock,
    int jmin, int jmax, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(mvec_view_block_f)(TYPE(mvec_ptr),TYPE(mvec_ptr)*,int,int,int*);
  PHIST_CHK_IERR(SUBR(mvec_view_block_f)(V,Vblock,jmin,jmax,ierr),*ierr);
}

void SUBR(mvec_view_scattered)(TYPE(mvec_ptr) V, TYPE(mvec_ptr)* Vscat,
    int* cols, int ncols, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=-99;
}

void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) V,
    TYPE(mvec_ptr) Vblock,
    int jmin, int jmax, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(mvec_get_block_f)(TYPE(const_mvec_ptr),TYPE(mvec_ptr),int,int,int*);
  PHIST_CHK_IERR(SUBR(mvec_get_block_f)(V,Vblock,jmin,jmax,ierr),*ierr);
}

void SUBR(mvec_set_block)(TYPE(mvec_ptr) V,
    TYPE(const_mvec_ptr) Vblock,
    int jmin, int jmax, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(mvec_set_block_f)(TYPE(mvec_ptr),TYPE(const_mvec_ptr),int,int,int*);
  PHIST_CHK_IERR(SUBR(mvec_set_block_f)(V,Vblock,jmin,jmax,ierr),*ierr);
}

void SUBR(sdMat_view_block)(TYPE(mvec_ptr) M, 
    TYPE(mvec_ptr)* Mblock,
    int imin, int imax, int jmin, int jmax, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(sdMat_view_block_f)(TYPE(mvec_ptr),TYPE(mvec_ptr)*,int,int,int,int,int*);
  PHIST_CHK_IERR(SUBR(sdMat_view_block_f)(M,Mblock,imin,imax,jmin,jmax,ierr),*ierr);
}

void SUBR(sdMat_get_block)(TYPE(const_mvec_ptr) M, 
    TYPE(mvec_ptr) Mblock,
    int imin, int imax, int jmin, int jmax, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(sdMat_get_block_f)(TYPE(const_mvec_ptr),TYPE(mvec_ptr),int,int,int,int,int*);
  PHIST_CHK_IERR(SUBR(sdMat_get_block_f)(M,Mblock,imin,imax,jmin,jmax,ierr),*ierr);
}

void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) M, 
    TYPE(const_sdMat_ptr) Mblock,
    int imin, int imax, int jmin, int jmax, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(sdMat_set_block_f)(TYPE(mvec_ptr),TYPE(const_mvec_ptr),int,int,int,int,int*);
  PHIST_CHK_IERR(SUBR(sdMat_set_block_f)(M,Mblock,imin,imax,jmin,jmax,ierr),*ierr);
}

void SUBR(crsMat_delete)(TYPE(crsMat_ptr) A, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(crsMat_delete_f)(TYPE(crsMat_ptr) A, int* ierr);
  PHIST_CHK_IERR(SUBR(crsMat_delete_f)(A,ierr), *ierr);
}

void SUBR(mvec_delete)(TYPE(mvec_ptr) V, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(mvec_delete_f)(TYPE(mvec_ptr),int*);
  PHIST_CHK_IERR( SUBR(mvec_delete_f) (V,ierr),*ierr);
}

void SUBR(sdMat_delete)(TYPE(sdMat_ptr) M, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(sdMat_delete_f)(TYPE(sdMat_ptr),int*);
  PHIST_CHK_IERR(SUBR(sdMat_delete_f)(M,ierr),*ierr);
}

void SUBR(mvec_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(mvec_put_value_f)(TYPE(mvec_ptr),_ST_,int*);
  PHIST_CHK_IERR(SUBR(mvec_put_value_f)(V,value,ierr),*ierr);
}

void SUBR(sdMat_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(sdMat_put_value_f)(TYPE(sdMat_ptr),_ST_,int*);
  PHIST_CHK_IERR(SUBR(sdMat_put_value_f)(V,value,ierr),*ierr);
}

void SUBR(mvec_random)(TYPE(mvec_ptr) V, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(mvec_random_f)(TYPE(mvec_ptr),int*);
  PHIST_CHK_IERR(SUBR(mvec_random_f)(V,ierr),*ierr);
}

void SUBR(mvec_print)(TYPE(const_mvec_ptr) V, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(mvec_print_f)(TYPE(const_mvec_ptr),int*);
  PHIST_CHK_IERR(SUBR(mvec_print_f)(V,ierr),*ierr);
}

void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) M, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(sdMat_print_f)(TYPE(const_sdMat_ptr),int*);
  PHIST_CHK_IERR(SUBR(sdMat_print_f)(M,ierr),*ierr);
}

void SUBR(sdMat_random)(TYPE(sdMat_ptr) M, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(sdMat_random_f)(TYPE(sdMat_ptr),int*);
  PHIST_CHK_IERR(SUBR(sdMat_random_f)(M,ierr),*ierr);
}

void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) V,
    _MT_* vnrm, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(mvec_norm2_f)(TYPE(const_mvec_ptr),_MT_*,int*);
  PHIST_CHK_IERR(SUBR(mvec_norm2_f)(V,vnrm,ierr),*ierr);
}

void SUBR(mvec_normalize)(TYPE(mvec_ptr) V,
    _MT_* vnrm, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"  
  PHIST_CHK_IERR(SUBR(mvec_norm2)(V,vnrm,ierr),*ierr);
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nvec,ierr),*ierr);
  _ST_*scale = new _ST_[nvec];
  for(int i = 0; i < nvec; i++)
    scale[i] = st::one()/vnrm[i];
  PHIST_CHK_IERR(SUBR(mvec_vscale)(V,scale,ierr),*ierr);
  delete[] scale;
}

void SUBR(mvec_scale)(TYPE(mvec_ptr) V, 
    _ST_ scalar, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(mvec_scale_f)(TYPE(mvec_ptr),_ST_,int*);
  PHIST_CHK_IERR(SUBR(mvec_scale_f)(V,scalar,ierr),*ierr);
}

void SUBR(mvec_vscale)(TYPE(mvec_ptr) V, 
    const _ST_* scalar, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(mvec_vscale_f)(TYPE(mvec_ptr),const _ST_*,int*);
  PHIST_CHK_IERR(SUBR(mvec_vscale_f)(V,scalar,ierr),*ierr);
}

void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) X,
    _ST_ beta,  TYPE(mvec_ptr)       Y, 
    int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(mvec_add_mvec_f)(_ST_,TYPE(const_mvec_ptr),_ST_,TYPE(mvec_ptr),int*);
  PHIST_CHK_IERR(SUBR(mvec_add_mvec_f)(alpha,X,beta,Y,ierr),*ierr);
}

void SUBR(mvec_vadd_mvec)(const _ST_ alpha[], TYPE(const_mvec_ptr) X,
    _ST_ beta,  TYPE(mvec_ptr)       Y, 
    int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(mvec_vadd_mvec_f)(const _ST_*,TYPE(const_mvec_ptr),_ST_,TYPE(mvec_ptr),int*);
  PHIST_CHK_IERR(SUBR(mvec_vadd_mvec_f)(alpha,X,beta,Y,ierr),*ierr);
}

void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) A,
    _ST_ beta,  TYPE(sdMat_ptr)       B, 
    int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(sdMat_add_sdMat_f)(_ST_,TYPE(const_sdMat_ptr),_ST_,TYPE(sdMat_ptr),int*);
  PHIST_CHK_IERR(SUBR(sdMat_add_sdMat_f)(alpha,A,beta,B,ierr),*ierr);
}

void SUBR(crsMat_times_mvec)(_ST_ alpha, TYPE(const_crsMat_ptr) A, 
    TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(crsMat_times_mvec_f)(_ST_ alpha, TYPE(const_crsMat_ptr) A, 
      TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* ierr);
  PHIST_CHK_IERR(SUBR(crsMat_times_mvec_f)(alpha,A,x,beta,y,ierr),*ierr);
}

void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) v, 
    TYPE(const_mvec_ptr) w, 
    _ST_* s, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(mvec_dot_mvec_f)(TYPE(const_mvec_ptr),TYPE(const_mvec_ptr),_ST_*,int*);
  PHIST_CHK_IERR(SUBR(mvec_dot_mvec_f)(v,w,s,ierr),*ierr);
}

void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
    TYPE(const_sdMat_ptr) C, 
    _ST_ beta, TYPE(mvec_ptr) W, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(mvec_times_sdMat_f)(_ST_,TYPE(const_mvec_ptr),TYPE(const_sdMat_ptr),_ST_,TYPE(mvec_ptr),int*);
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat_f)(alpha,V,C,beta,W,ierr),*ierr);
}

void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
    TYPE(const_sdMat_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(sdMat_times_sdMat_f)(_ST_,TYPE(const_sdMat_ptr),TYPE(const_sdMat_ptr),_ST_,TYPE(sdMat_ptr),int*);
  PHIST_CHK_IERR(SUBR(sdMat_times_sdMat_f)(alpha,V,W,beta,C,ierr),*ierr);
}

void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
    TYPE(const_sdMat_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(sdMatT_times_sdMat_f)(_ST_,TYPE(const_sdMat_ptr),TYPE(const_sdMat_ptr),_ST_,TYPE(sdMat_ptr),int*);
  PHIST_CHK_IERR(SUBR(sdMatT_times_sdMat_f)(alpha,V,W,beta,C,ierr),*ierr);
}


void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
    TYPE(const_mvec_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  void SUBR(mvecT_times_mvec_f)(_ST_,TYPE(const_mvec_ptr),TYPE(const_mvec_ptr),_ST_,TYPE(sdMat_ptr),int*);
  PHIST_CHK_IERR(SUBR(mvecT_times_mvec_f)(alpha,V,W,beta,C,ierr),*ierr);
}

void SUBR(mvec_QR)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  // trilinos tsqr would be better but seems to expect column-major mvecs
  void SUBR(mvec_QR_f)(TYPE(mvec_ptr),TYPE(sdMat_ptr),int*);
  PHIST_CHK_NEG_IERR(SUBR(mvec_QR_f)(V,R,ierr),*ierr);
}

}