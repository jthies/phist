// some kernel functions are straightforwardly implemented based on others,
// so we provide a central implementation for these non-performance critical
// things.

//! get global sparseMat size (number of rows) \ingroup crsmat
extern "C" void SUBR(sparseMat_global_nrows)(TYPE(sparseMat_ptr) A, gidx_t* s, int* iflag)
{
  const_map_ptr_t map=NULL;
  PHIST_CHK_IERR(SUBR(sparseMat_get_range_map)(A,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(map,s,iflag),*iflag);
}

//! get global sparseMat size (number of columns) \ingroup crsmat
extern "C" void SUBR(sparseMat_global_ncols)(TYPE(sparseMat_ptr) A, gidx_t* s, int* iflag)
{
  const_map_ptr_t map=NULL;
  PHIST_CHK_IERR(SUBR(sparseMat_get_domain_map)(A,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(map,s,iflag),*iflag);
}

//! retrieve local length of the vectors in V \ingroup mvec
extern "C" void SUBR(mvec_my_length)(TYPE(const_mvec_ptr) V, lidx_t* len, int* iflag)
{
  const_map_ptr_t map=NULL;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(V,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_local_length(map,len,iflag),*iflag);
}

//! retrieve global length of the vectors in V \ingroup mvec
extern "C" void SUBR(mvec_global_length)(TYPE(const_mvec_ptr) V, gidx_t* len, int* iflag)
{
  const_map_ptr_t map=NULL;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(V,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(map,len,iflag),*iflag);
}

//! retrieve the comm used for MPI communication in V \ingroup mvec
extern "C" void SUBR(mvec_get_comm)(TYPE(const_mvec_ptr) V, const_comm_ptr_t* comm, int* iflag)
{
  const_map_ptr_t map=NULL;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(V,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_comm(map,comm,iflag),*iflag);
}

//! y[i]=alpha*(A*x+shift*x) + beta*y
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
#ifdef PHIST_BUILTIN_RNG

int PREFIX(copyDataFunc)(ghost_gidx_t i, ghost_lidx_t j, void* vval,void* vdata)
{
  dwrap* wrap=(dwrap*)vdata;
  int lda = wrap->lda;
  int ii = i - wrap->ilower;
  double* data = (double*)vdata;
  _MT_* val = (_MT_*)vval;
  
  if (ii>=wrap->lnrows || j>=wrap->lncols)
  {
    return -1; // index out of bounds;
  }
  
#ifdef IS_COMPLEX
  val[0]=data[i*lda+2*j];
  val[1]=(_MT_)data[ii*lda+2*j+1];
#else
  val[0]=(_MT_)data[ii*lda+2*j];
#endif
return 0;
}

extern "C" void SUBR(mvec_random)(TYPE(mvec_ptr) V, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  gidx_t gnrows,ilower,iupper,pre_skip,post_skip;
  const_map_ptr_t map=NULL;
  lidx_t lnrows,nvec;
  
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nvec,iflag),*iflag);

  PHIST_CHK_IERR(SUBR(mvec_get_map)(V,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_local_length(map,&lnrows,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(map,&gnrows,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_ilower(map,&ilower,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_iupper(map,&iupper,iflag),*iflag);
  
  pre_skip = ilower*nvec;
  post_skip= (gnrows-iupper)*nvec;  
  
#ifdef IS_COMPLEX
  const int nelem=2;
#else
  const int nelem=1;
#endif  
  
  // we use the most robust way of implementing this, which should work for
  // any situation (row/col major, GPU/CPU etc.): generate row-major clone data
  // and set the vector elements using mvec_put_func.
  lidx_t lda=lnrows*nelem;
  double *randbuf;
  *iflag = posix_memalign(&randbuf, 64, nvec*lda*sizeof(double));
  if (*iflag!=0)
  {
    *iflag=PHIST_MEM_ALLOC_FAILED;
    return;
  }
                      
  drandom_general(nelem*nvec,(int)lnrows, randbuf,(int)lda,(int64_t)pre_skip, (int64_t)post_skip);
  dwrap wrap;
  wrap.lda=lda;
  wrap.lnrows=lnrows;
  wrap.lncols=nvec;
  wrap.ilower=ilower;
  wrap.data=randbuf;
  PHIST_CHK_IERR(SUBR(mvec_put_func)(V,&PREFIX(copyDataFunc),&wrap,iflag),*iflag);
  free(randbuf);
}
#endif
