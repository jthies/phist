/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
// some kernel functions are straightforwardly implemented based on others,
// so we provide a central implementation for these non-performance critical
// things.

// get global sparseMat size (number of rows) \ingroup crsmat
extern "C" void SUBR(sparseMat_global_nrows)(TYPE(sparseMat_ptr) A, phist_gidx* s, int* iflag)
{
  phist_const_map_ptr map=NULL;
  PHIST_CHK_IERR(SUBR(sparseMat_get_range_map)(A,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(map,s,iflag),*iflag);
}

// get global sparseMat size (number of columns) \ingroup crsmat
extern "C" void SUBR(sparseMat_global_ncols)(TYPE(sparseMat_ptr) A, phist_gidx* s, int* iflag)
{
  phist_const_map_ptr map=NULL;
  PHIST_CHK_IERR(SUBR(sparseMat_get_domain_map)(A,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(map,s,iflag),*iflag);
}

// retrieve local length of the vectors in V \ingroup mvec
extern "C" void SUBR(mvec_my_length)(TYPE(const_mvec_ptr) V, phist_lidx* len, int* iflag)
{
  phist_const_map_ptr map=NULL;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(V,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_local_length(map,len,iflag),*iflag);
}

// retrieve global length of the vectors in V \ingroup mvec
extern "C" void SUBR(mvec_global_length)(TYPE(const_mvec_ptr) V, phist_gidx* len, int* iflag)
{
  phist_const_map_ptr map=NULL;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(V,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(map,len,iflag),*iflag);
}

// retrieve the comm used for MPI communication in V \ingroup mvec
extern "C" void SUBR(mvec_get_comm)(TYPE(const_mvec_ptr) V, phist_const_comm_ptr* comm, int* iflag)
{
  phist_const_map_ptr map=NULL;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(V,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_comm(map,comm,iflag),*iflag);
}

// y[i]=alpha*(A*x+shift*x) + beta*y
extern "C" void SUBR(sparseMat_times_mvec_add_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A,
        _ST_ shift, TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  int nv;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(x,&nv,iflag),*iflag);
  _ST_* shifts=(_ST_*)malloc(nv*sizeof(_ST_));
  for (int i=0;i<nv;i++) shifts[i]=shift;
  SUBR(sparseMat_times_mvec_vadd_mvec)(alpha,A,shifts,x,beta,y,iflag);
  free(shifts);
}

// create a new mvec with the same dimensions (number of rows and columns) and
// distribution (map)  as an existing one. The values of the new object are not
// initialized explicitly, so if you want to clone the vector contents as well,
// you will have to call mvec_set_block afterwards (or similar). *V must be NULL
// on input.
extern "C" void SUBR(mvec_clone_shape)(TYPE(mvec_ptr)* V, TYPE(const_mvec_ptr) V_in, int* iflag)
{
  int iflag_in=*iflag;
  *iflag=0;
  phist_const_map_ptr map=NULL;
PHIST_CHK_IERR(SUBR(mvec_get_map)(V_in,&map,iflag),*iflag);
int nvecs;
PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V_in,&nvecs,iflag),*iflag);
*iflag=iflag_in;
PHIST_CHK_IERR(SUBR(mvec_create)(V,map,nvecs,iflag),*iflag);
}

#ifdef PHIST_BUILTIN_RNG

int PHIST_TG_PREFIX(copyDataFunc)(ghost_gidx i, ghost_lidx j, void* vval,void* vdata)
{
  dwrap* wrap=(dwrap*)vdata;
  int lda = wrap->lda;
  int ii = i - wrap->ilower;
  _MT_* val = (_MT_*)vval;
  
  if (ii>=wrap->lnrows || j>=wrap->lncols)
  {
    return -1; // index out of bounds;
  }
  
#ifdef IS_COMPLEX
  val[0]=(_MT_)wrap->data[ii*lda+2*j];
  val[1]=(_MT_)wrap->data[ii*lda+2*j+1];
#else
  //PHIST_SOUT(PHIST_INFO,"copyDataFunc %d %d",(int)i,(int)j);
  val[0]=(_MT_)wrap->data[ii*lda+j];
  //PHIST_SOUT(PHIST_INFO," %8.4e\n", val[0]);
#endif
return 0;
}

extern "C" void SUBR(mvec_random)(TYPE(mvec_ptr) V, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  phist_const_comm_ptr comm;
  phist_gidx gnrows,ilower,iupper,pre_skip,post_skip;
  phist_const_map_ptr map=NULL;
  phist_lidx lnrows,nvec;
  
  bool is_linear_map;
  
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nvec,iflag),*iflag);

  PHIST_CHK_IERR(SUBR(mvec_get_map)(V,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_local_length(map,&lnrows,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(map,&gnrows,iflag),*iflag);
  // these may return iflag=1 if the map is not a linear map. In that
  // case we can still create a reproducible sequence by skipping a little 
  // further in the random number stream (see next if statement)
  PHIST_CHK_NEG_IERR(phist_map_get_ilower(map,&ilower,iflag),*iflag);
  is_linear_map=(*iflag==0);
  PHIST_CHK_NEG_IERR(phist_map_get_iupper(map,&iupper,iflag),*iflag);
  is_linear_map&=(*iflag==0);
  
#ifdef IS_COMPLEX
  const int nelem=2;
#else
  const int nelem=1;
#endif  

  // deal with non-standard (linear) maps
  if (is_linear_map==false)
  {
    int rank, size;
    PHIST_CHK_IERR(phist_map_get_comm(map,&comm,iflag),*iflag);
    PHIST_CHK_IERR(phist_comm_get_rank(comm,&rank,iflag),*iflag);
    PHIST_CHK_IERR(phist_comm_get_size(comm,&size,iflag),*iflag);
    pre_skip  = rank            * gnrows           *nvec*nelem;
    post_skip = ((size-rank)    * gnrows - lnrows) *nvec*nelem;
  }
  else
  {
    pre_skip = ilower*nvec*nelem;
    post_skip= (gnrows-iupper)*nvec*nelem;
  }
  
    
  // we use the most robust way of implementing this, which should work for
  // any situation (row/col major, GPU/CPU etc.): generate row-major clone data
  // and set the vector elements using mvec_put_func.
  double *randbuf;
  size_t sz=lnrows*nvec*nelem;
  *iflag = posix_memalign((void**)&randbuf, 64, sz*sizeof(double));

  if (*iflag!=0)
  {
    *iflag=PHIST_MEM_ALLOC_FAILED;
    return;
  }
                      
  drandom_1(sz, randbuf,(int64_t)pre_skip, (int64_t)post_skip);
 
  /*
 for (int i=0; i<lnrows; i++)
 {
   PHIST_SOUT(PHIST_INFO,"%d",i);
   for (int j=0; j<lda; j++)
   {
     PHIST_SOUT(PHIST_INFO,"  %8.4e",randbuf[i*lda+j]);
   }
 PHIST_SOUT(PHIST_INFO,"\n");
 }
 */
  dwrap wrap;
  wrap.lda=nvec*nelem;
  wrap.lnrows=lnrows;
  wrap.lncols=nvec;
  wrap.ilower=ilower;
  wrap.data=randbuf;
  if (is_linear_map)
  {
    PHIST_CHK_IERR(SUBR(mvec_put_func)(V,&PHIST_TG_PREFIX(copyDataFunc),&wrap,iflag),*iflag);
  }
  else
  {
    // since put_func with copyDataFunc won't work because we can't reconstruct the local index of the given
    // global one by subtracting ilower. So instead we copy the data manually. This is not a safe way of doing
    // it because this may be a "device rank", but so far we only support GPUs with GHOST, and GHOST only has
    // linearly distributed indices.
    _ST_* V_raw=NULL;
    phist_lidx ldV;
    PHIST_CHK_IERR(SUBR(mvec_extract_view)(V,&V_raw,&ldV,iflag),*iflag);
    ghost_lidx pos_V, pos_buf;
    for (ghost_lidx i=0; i<lnrows; i++)
    {
      for (int j=0; j<nvec; j++)
      {
        pos_buf=i*wrap.lda+j;
#ifdef PHIST_MVECS_ROW_MAJOR
        pos_V=i*ldV+j;
#else
        pos_V=j*ldV+i;
#endif
        V_raw[pos_V] = (_ST_)((_MT_)wrap.data[pos_buf]);
#ifdef IS_COMPLEX
        V_raw[pos_V] += (_ST_)((_MT_)wrap.data[pos_buf+1]*st::cmplx_I());
#endif
      }
    }
  }
  free(randbuf);
}

extern "C" void SUBR(sdMat_random)(TYPE(sdMat_ptr) M, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  phist_lidx nrows,ncols,lda;
  _ST_* M_raw = NULL;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(M,&nrows,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(M,&ncols,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(M,&M_raw,&lda,iflag),*iflag);
      
#ifdef IS_COMPLEX
  const int nelem=2;
#else
  const int nelem=1;
#endif

  // generate random doubles and copy them by hand, then upload to GPU if applicable.
  double *randbuf = NULL;
  *iflag = posix_memalign((void**)&randbuf, 64, nrows*ncols*nelem*sizeof(double));
  if (*iflag!=0)
  {
    *iflag=PHIST_MEM_ALLOC_FAILED;
    return;
  }
  phist_Drandom_number(nelem*nrows*ncols, randbuf);

  for (int j=0; j<ncols; j++)
  {
    for (int i=0; i<nrows; i++)
    {
#ifdef PHIST_SDMATS_ROW_MAJOR
# ifdef IS_COMPLEX
      M_raw[lda*i+j] = (_MT_)randbuf[nelem*(j*nrows+i)] + st::cmplx_I()*(_MT_)randbuf[nelem*(j*nrows+i)+1];
# else
      M_raw[lda*i+j] = (_MT_)randbuf[nelem*(j*nrows+i)];
# endif
#else
# ifdef IS_COMPLEX
      M_raw[lda*j+i] = (_MT_)randbuf[nelem*(j*nrows+i)] + st::cmplx_I()*(_MT_)randbuf[nelem*(j*nrows+i)+1];
# else
      M_raw[lda*j+i] = (_MT_)randbuf[nelem*(j*nrows+i)];
# endif
#endif
    }
  }
 
  PHIST_CHK_IERR(SUBR(sdMat_to_device)(M,iflag),*iflag);
  free(randbuf);
}
#endif
