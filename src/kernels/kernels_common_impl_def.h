// some kernel functions are straightforwardly implemented based on others,
// so we provide a central implementation for these non-performance critical
// things.

//! get global sparseMat size (number of rows) \ingroup crsmat
void SUBR(sparseMat_global_nrows)(TYPE(sparseMat_ptr) A, gidx_t* s, int* iflag)
{
  const_map_ptr_t map=NULL;
  PHIST_CHK_IERR(SUBR(sparseMat_get_range_map)(A,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(map,s,iflag),*iflag);
}

//! get global sparseMat size (number of columns) \ingroup crsmat
void SUBR(sparseMat_global_ncols)(TYPE(sparseMat_ptr) A, gidx_t* s, int* iflag)
{
  const_map_ptr_t map=NULL;
  PHIST_CHK_IERR(SUBR(sparseMat_get_domain_map)(A,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(map,s,iflag),*iflag);
}

//! retrieve local length of the vectors in V \ingroup mvec
void SUBR(mvec_my_length)(TYPE(const_mvec_ptr) V, lidx_t* len, int* iflag)
{
  const_map_ptr_t map=NULL;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(V,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_local_length(map,len,iflag),*iflag);
}

//! retrieve global length of the vectors in V \ingroup mvec
void SUBR(mvec_global_length)(TYPE(const_mvec_ptr) V, gidx_t* len, int* iflag)
{
  const_map_ptr_t map=NULL;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(V,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(map,len,iflag),*iflag);
}

//! retrieve the comm used for MPI communication in V \ingroup mvec
void SUBR(mvec_get_comm)(TYPE(const_mvec_ptr) V, const_comm_ptr_t* comm, int* iflag)
{
  const_map_ptr_t map=NULL;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(V,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_comm(map,comm,iflag),*iflag);
}

//! y[i]=alpha*(A*x+shift*x) + beta*y
void SUBR(sparseMat_times_mvec_add_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A,
        _ST_ shift, TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  int nv;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(x,&nv,iflag),*iflag);
  _ST_* shifts=(_ST_*)malloc(nv*sizeof(_ST_));;
  for (int i=0;i<nv;i++) shifts[i]=shift;
  SUBR(sparseMat_times_mvec_vadd_mvec)(alpha,A,shifts,x,beta,y,iflag);
  free(shifts);
}
