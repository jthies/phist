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

// augmented spMVM with single shift
extern "C" void SUBR(sparseMat_times_mvec_aug)(_ST_ alpha, TYPE(const_sparseMat_ptr) A,
        _ST_ shift, TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, 
        _ST_ a, _ST_ b, TYPE(mvec_ptr) z,
        _ST_* dot_xx, _ST_* dot_xy, _ST_* dot_yy, int* iflag)
{
  int nv;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(x,&nv,iflag),*iflag);
  _ST_* shifts=(_ST_*)malloc(nv*sizeof(_ST_));
  for (int i=0;i<nv;i++) shifts[i]=shift;
  SUBR(sparseMat_times_mvec_vaug)(alpha,A,shifts,x,beta,y,
        a,b,z,dots_xx,dots_xy,dots_yy,iflag);
  free(shifts);
}
#ifdef PHIST_BUILTIN_RNG

int PREFIX(copyDataFunc)(ghost_gidx_t i, ghost_lidx_t j, void* vval,void* vdata)
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
  
#ifdef IS_COMPLEX
  const int nelem=2;
#else
  const int nelem=1;
#endif  

  pre_skip = ilower*nvec*nelem;
  post_skip= (gnrows-iupper)*nvec*nelem;
    
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
  PHIST_CHK_IERR(SUBR(mvec_put_func)(V,&PREFIX(copyDataFunc),&wrap,iflag),*iflag);
  free(randbuf);
}

extern "C" void SUBR(sdMat_random)(TYPE(sdMat_ptr) M, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  lidx_t nrows,ncols,lda;
  _ST_* M_raw;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(M,&nrows,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(M,&ncols,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(M,&M_raw,&lda,iflag),*iflag);
      
#ifdef IS_COMPLEX
  const int nelem=2;
#else
  const int nelem=1;
#endif

  // generate random doubles and copy them by hand, then upload to GPU if applicable.
  double *randbuf;
  *iflag = posix_memalign((void**)&randbuf, 64, nrows*ncols*nelem*sizeof(double));
  if (*iflag!=0)
  {
    *iflag=PHIST_MEM_ALLOC_FAILED;
    return;
  }
  phist_Drandom_number(nelem*nrows*ncols, randbuf);
  _MT_ *mM_raw=(_MT_*)M_raw;
  for (int j=0; j<ncols; j++)
  {
    for (int i=0; i<nrows; i++)
    {
#ifdef PHIST_SDMATS_ROW_MAJOR
      mM_raw[lda*nelem*i+j]=(_MT_)randbuf[j*nelem*nrows+i];
# ifdef IS_COMPLEX
      mM_raw[lda*nelem*i+j+1]=(_MT_)randbuf[j*nelem*nrows+i+1];
# endif
#else
      mM_raw[lda*nelem*j+i]=(_MT_)randbuf[j*nelem*nrows+i];
# ifdef IS_COMPLEX
      mM_raw[lda*nelem*j+i+1]=(_MT_)randbuf[j*nelem*nrows+i+1];
# endif
#endif
    }
  }
 
  PHIST_CHK_IERR(SUBR(sdMat_to_device)(M,iflag),*iflag);
  free(randbuf);
}
#endif
