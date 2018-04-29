/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
/*! \file Eigen/kernels_def.hpp
 * included by Eigen/kernels.cpp
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 * \author "Jonas Thies <Jonas.Thies@DLR.de>
 *
*/

extern "C" void SUBR(type_avail)(int *iflag)
{
  *iflag=0;
}

extern "C" void SUBR(sparseMat_read_mm)(TYPE(sparseMat_ptr)* A, phist_const_comm_ptr vcomm,
        const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  if (filename==NULL)
  {
    *iflag=PHIST_INVALID_INPUT;
    return;
  }

  std::string line;
  std::ifstream infile(filename);
  getline(infile, line);
#ifdef IS_COMPLEX
  if( line != "%%MatrixMarket matrix coordinate real general" &&  line != "%%MatrixMarket matrix coordinate real symmetric" &&
      line != "%%MatrixMarket matrix coordinate complex general" &&  line != "%%MatrixMarket matrix coordinate complex symmetric" )
#else
  if( line != "%%MatrixMarket matrix coordinate real general" &&  line != "%%MatrixMarket matrix coordinate real symmetric" )
#endif
  {
    PHIST_SOUT(PHIST_ERROR, "unexpected first line in .mm file: %s\n", line.c_str());
    *iflag = PHIST_NOT_IMPLEMENTED;
    return;
  }

  while(infile.peek() == '%')
    getline(infile, line);

  phist_gidx globalRows, globalCols, globalLines;
  infile >> globalRows >> globalCols >> globalLines;
  phist_map_ptr map;
  phist_context_ptr ctx;
  PHIST_CHK_IERR(phist_map_create(&map,vcomm,globalRows,iflag),*iflag);
  PHIST_CHK_IERR(phist_context_create(&ctx,map,NULL,map,map,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sparseMat_read_mm_with_context)(A, ctx, filename, iflag), *iflag);

  *iflag=PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_read_bin)(TYPE(sparseMat_ptr)* A, phist_const_comm_ptr vcomm,
const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_hb)(TYPE(sparseMat_ptr)* A, phist_const_comm_ptr vcomm,
const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_mm_with_context)(TYPE(sparseMat_ptr)* vA, phist_const_context_ptr vctx,
        const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  if (filename==NULL)
  {
    *iflag=PHIST_INVALID_INPUT;
    return;
  }

  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sparseMat_t*,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(phist::internal::default_context,ctx,vctx,*iflag);

  std::string line;
  std::ifstream infile(filename);
  getline(infile, line);
  bool symm;
  bool cmplx = false;
  if( line == "%%MatrixMarket matrix coordinate real general" )
    symm = false;
  else if( line == "%%MatrixMarket matrix coordinate real symmetric" )
    symm = true;
#ifdef IS_COMPLEX
  else if( line == "%%MatrixMarket matrix coordinate complex general" )
  {
    symm = false;
    cmplx = true;
  }
  else if( line == "%%MatrixMarket matrix coordinate complex symmetric" )
  {
    symm = true;
    cmplx = true;
  }
#endif
  else
  {
    PHIST_SOUT(PHIST_ERROR, "unexpected first line in .mm file: %s\n", line.c_str());
    *iflag = PHIST_NOT_IMPLEMENTED;
    return;
  }

  while(infile.peek() == '%')
    getline(infile, line);

  phist_gidx globalRows, globalCols, globalLines;
  infile >> globalRows >> globalCols >> globalLines;
  phist_gidx globalEntries = globalLines;
  if( symm )
    globalEntries += globalLines - globalRows;
  phist_lidx avg_nne = globalEntries/globalRows;

  phist_lidx nlocal;
  phist_gidx ilower, iupper, nglob;
  phist_const_comm_ptr comm;
  phist_const_map_ptr map=ctx->row_map;
  PHIST_CHK_IERR(phist_map_get_local_length(map, &nlocal, iflag), *iflag);
  PHIST_CHK_IERR(phist_map_get_ilower(map, &ilower, iflag), *iflag);
  PHIST_CHK_IERR(phist_map_get_iupper(map, &iupper, iflag), *iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(map, &nglob, iflag), *iflag);
  PHIST_CHK_IERR(phist_map_get_comm(map, &comm, iflag), *iflag);

  *A = new Traits<_ST_>::sparseMat_t;
  (*A)->m.resize(nlocal, nglob);
  (*A)->m.reserve(avg_nne*nlocal);
  (*A)->map = map;
  PHIST_SOUT(PHIST_INFO, "reading mm mat: %d x %d\n", nglob, nglob);

  for (phist_gidx i = 0; i < globalLines; i++)
  {
    _ST_ val;
    phist_gidx row, col;
#ifdef IS_COMPLEX
    _MT_ val_r, val_i = 0;
    if( cmplx )
      infile >> row >> col >> val_r >> val_i;
    else
      infile >> row >> col >> val_r;
    val = std::complex<_MT_>(val_r,val_i);
#else
    infile >> row >> col >> val;
#endif
    // count from zero...
    row--, col--;
    if( ilower <= row && row <= iupper )
    {
      (*A)->m.insert(row,col)=val;
    }
    if( symm && col != row && ilower <= col && col <= iupper )
    {
      (*A)->m.insert(col,row)=val;
    }
  }

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_read_bin_with_context)(TYPE(sparseMat_ptr)* A, phist_const_context_ptr ctx,
        const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_hb_with_context)(TYPE(sparseMat_ptr)* A, phist_const_context_ptr ctx,
        const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_local_nnz)(TYPE(const_sparseMat_ptr) vA, int64_t* local_nnz, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t,A,vA,*iflag);
  *local_nnz=A->m.nonZeros();
  *iflag=0;
}

extern "C" void SUBR(sparseMat_global_nnz)(TYPE(const_sparseMat_ptr) vA, int64_t* global_nnz, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t,A,vA,*iflag);
  *global_nnz=A->m.nonZeros();
  *iflag=0;
}

extern "C" void SUBR(sparseMat_get_row_map)(TYPE(const_sparseMat_ptr) vA, phist_const_map_ptr* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t,A,vA,*iflag);

  *map = A->map;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_get_col_map)(TYPE(const_sparseMat_ptr) vA, phist_const_map_ptr* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t,A,vA,*iflag);

  *map = A->map;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_get_domain_map)(TYPE(const_sparseMat_ptr) vA, phist_const_map_ptr* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t,A,vA,*iflag);

  *map = A->map;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_get_range_map)(TYPE(const_sparseMat_ptr) vA, phist_const_map_ptr* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t,A,vA,*iflag);

  *map = A->map;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_create)(TYPE(mvec_ptr)* vV, 
    phist_const_map_ptr map, int nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t*,V,vV,*iflag);

  phist_lidx nlocal;
  phist_gidx nglob;
  phist_const_comm_ptr comm;
  PHIST_CHK_IERR(phist_map_get_local_length(map,&nlocal,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(map,&nglob,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_comm(map,&comm,iflag),*iflag);

  *V = new Traits<_ST_>::mvec_t;
  (*V)->v_storage.resize(nlocal, nvec);
  (*V)->map = map;
  // placement new (construct new object!)
  new (&((*V)->v)) Traits<_ST_>::mvec_t::Eigen_Ref((*V)->v_storage);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_create)(TYPE(sdMat_ptr)* vM, 
    int nrows, int ncols, phist_const_comm_ptr comm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t*,M,vM,*iflag);

  *M = new Traits<_ST_>::sdMat_t;
  (*M)->m_storage.resize(nrows, ncols);
  (*M)->comm = comm;
  // placement new (construct new object!)
  new (&((*M)->m)) Traits<_ST_>::sdMat_t::Eigen_Ref((*M)->m_storage);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) vV, phist_const_map_ptr* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  *map = V->map;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) vV, int* nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);

  *nvec = V->v.cols();
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) vM, int* nrows, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,M,vM,*iflag);

  *nrows = M->m.rows();
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) vM, int* ncols, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,M,vM,*iflag);

  *ncols = M->m.cols();
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_extract_view)(TYPE(mvec_ptr) vV, _ST_** val, phist_lidx* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);

  *val = V->v.data();
  *lda = V->v.outerStride();
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) vM, _ST_** val, phist_lidx* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);

  *val = M->m.data();
  *lda = M->m.outerStride();
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_to_mvec)(TYPE(const_mvec_ptr) vV, TYPE(mvec_ptr) vW, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,W,vW,*iflag);

  W->v.noalias() = V->v;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_view_block)(TYPE(mvec_ptr) vV,
    TYPE(mvec_ptr)* vVblock,
    int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t*,Vblock,vVblock,*iflag);

  // placement new (construct new object!)
  phist_lidx nvec = jmax-jmin+1;
  if( *Vblock == NULL )
    *Vblock = new Traits<_ST_>::mvec_t;
  (*Vblock)->map = V->map;
  new (&((*Vblock)->v)) Traits<_ST_>::mvec_t::Eigen_Ref(V->v.block(0,jmin,V->v.rows(),nvec));
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) vV,
    TYPE(mvec_ptr) vVblock,
    int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,Vblock,vVblock,*iflag);

  int nvec = jmax-jmin+1;
  Vblock->v = V->v.block(0,jmin,V->v.rows(),nvec);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_set_block)(TYPE(mvec_ptr) vV,
    TYPE(const_mvec_ptr) vVblock,
    int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,Vblock,vVblock,*iflag);

  int nvec = jmax-jmin+1;
  V->v.block(0,jmin,Vblock->v.rows(),nvec) = Vblock->v;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_view_block)(TYPE(sdMat_ptr) vM, 
    TYPE(sdMat_ptr)* vMblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t*,Mblock,vMblock,*iflag);

  phist_lidx nrows = imax-imin+1;
  phist_lidx ncols = jmax-jmin+1;
  if( *Mblock == NULL )
    (*Mblock) = new Traits<_ST_>::sdMat_t;
  new (&((*Mblock)->m)) Traits<_ST_>::sdMat_t::Eigen_Ref(M->m.block(imin,jmin,nrows,ncols));
  (*Mblock)->comm = M->comm;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_get_block)(TYPE(const_sdMat_ptr) vM, 
    TYPE(sdMat_ptr) vMblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,M,vM,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,Mblock,vMblock,*iflag);

  int nrows = imax-imin+1;
  int ncols = jmax-jmin+1;
  Mblock->m = M->m.block(imin,jmin,nrows,ncols);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) vM, 
    TYPE(const_sdMat_ptr) vMblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,Mblock,vMblock,*iflag);

  int nrows = imax-imin+1;
  int ncols = jmax-jmin+1;
  M->m.block(imin,jmin,nrows,ncols) = Mblock->m;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_delete)(TYPE(sparseMat_ptr) vA, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  if( vA == NULL )
  {
    *iflag = PHIST_SUCCESS;
    return;
  }
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sparseMat_t,A,vA,*iflag);

  delete A;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_delete)(TYPE(mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  if( vV == NULL )
  {
    *iflag = PHIST_SUCCESS;
    return;
  }
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);

  delete V;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_delete)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  if( vM == NULL )
  {
    *iflag = PHIST_SUCCESS;
    return;
  }
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);

  delete M;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_put_value)(TYPE(mvec_ptr) vV, _ST_ value, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);
  
  V->v = Traits<_ST_>::mvec_t::Eigen_MVec::Constant(V->v.rows(),V->v.cols(),value);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_put_func)(TYPE(mvec_ptr) vV,
        phist_mvec_elemFunc funPtr,void* last_arg, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);

  // no real support for this in Eigen...
  int nvec;
  phist_lidx nlocal;
  phist_gidx ilower;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vV, &nvec, iflag), *iflag);
  PHIST_CHK_IERR(phist_map_get_local_length(V->map,&nlocal,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_ilower(V->map,&ilower,iflag),*iflag);
  for(phist_lidx j = 0; j < nvec; j++)
  {
    for(phist_lidx i = 0; i < nlocal; i++)
    {
      PHIST_CHK_IERR(*iflag=funPtr(ilower+i,j,(void*)&(V->v(i,j)),last_arg),*iflag);
    }
  }
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_put_value)(TYPE(sdMat_ptr) vM, _ST_ value, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);

  M->m = Traits<_ST_>::sdMat_t::Eigen_SdMat::Constant(M->m.rows(),M->m.cols(),value);
  *iflag = PHIST_SUCCESS;
}

#ifndef PHIST_BUILTIN_RNG
extern "C" void SUBR(mvec_random)(TYPE(mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);

  V->v = Traits<_ST_>::mvec_t::Eigen_MVec::Random(V->v.rows(),V->v.cols());
  *iflag = PHIST_SUCCESS;
}
#endif

extern "C" void SUBR(mvec_print)(TYPE(const_mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  
  std::cout << V->v << std::endl;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,M,vM,*iflag);

  std::cout << M->m << std::endl;
  *iflag = PHIST_SUCCESS;
}

#ifndef PHIST_BUILTIN_RNG
extern "C" void SUBR(sdMat_random)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);

  M->m = Traits<_ST_>::sdMat_t::Eigen_SdMat::Random(M->m.rows(),M->m.cols());
  *iflag = PHIST_SUCCESS;
}
#endif

extern "C" void SUBR(sdMat_identity)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);

  M->m = Traits<_ST_>::sdMat_t::Eigen_SdMat::Identity(M->m.rows(),M->m.cols());
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) vV,
    _MT_* vnrm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);

  Eigen::Map< Eigen::Array<_MT_,Eigen::Dynamic,1> > nrm(vnrm,V->v.cols());
  nrm = V->v.colwise().norm();
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_normalize)(TYPE(mvec_ptr) vV,
    _MT_* vnrm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);

  Eigen::Map< Eigen::Array<_MT_,Eigen::Dynamic,1> > nrm(vnrm,V->v.cols());
  nrm = V->v.colwise().norm();
  Eigen::Matrix<_MT_,Eigen::Dynamic,1> inv_nrm(V->v.cols());
  inv_nrm = 1./nrm;
  V->v *= inv_nrm.asDiagonal();
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_scale)(TYPE(mvec_ptr) vV, 
    _ST_ scalar, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);

  V->v *= scalar;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_vscale)(TYPE(mvec_ptr) vV, 
    const _ST_* scalar, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);

  Eigen::Map< const Eigen::Matrix<_ST_,Eigen::Dynamic,1> > scale(scalar,V->v.cols());
  V->v *= scale.asDiagonal();
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) vV,
                                    _ST_ beta,  TYPE(mvec_ptr)       vW,
                                    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,W,vW,*iflag);

  W->v = alpha*V->v + beta*W->v;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_vadd_mvec)(const _ST_ alpha[], TYPE(const_mvec_ptr) vV,
                                    _ST_ beta,  TYPE(mvec_ptr)       vW,
                                    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,W,vW,*iflag);

  Eigen::Map< const Eigen::Matrix<_ST_,Eigen::Dynamic,1> > alpha_(alpha,V->v.cols());
  W->v = V->v*alpha_.asDiagonal() + beta*W->v;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA,
                                      _ST_ beta,  TYPE(sdMat_ptr)       vB, 
                                      int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,B,vB,*iflag);

  B->m = alpha*A->m + beta*B->m;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMatT_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA,
    _ST_ beta,  TYPE(sdMat_ptr)       vB, 
    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,B,vB,*iflag);

  B->m = alpha*A->m.adjoint() + beta*B->m;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_times_mvec_communicate)(TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag = 0;
}

extern "C" void SUBR(sparseMat_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) vA, 
    TYPE(const_mvec_ptr) vV, _ST_ beta, TYPE(mvec_ptr) vW, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,W,vW,*iflag);

  W->v = alpha*A->m*V->v + beta*W->v;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMatT_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) vA, 
    TYPE(const_mvec_ptr) vV, _ST_ beta, TYPE(mvec_ptr) vW, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,W,vW,*iflag);

  W->v = alpha*A->m.adjoint()*V->v + beta*W->v;
  *iflag = PHIST_SUCCESS;
}

//! y[i]=alpha*(A*x[i]+shifts[i]*x[i]) + beta*y[i]
extern "C" void SUBR(sparseMat_times_mvec_vadd_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) vA,
        const _ST_ shifts[], TYPE(const_mvec_ptr) vV, _ST_ beta, TYPE(mvec_ptr) vW, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,W,vW,*iflag);

  const Eigen::Map< const Eigen::Matrix<_ST_,Eigen::Dynamic,1> > shifts_(shifts,V->v.cols());
  W->v = alpha*(A->m*V->v+V->v*shifts_.asDiagonal()) + beta*W->v;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) vV, 
                                    TYPE(const_mvec_ptr) vW, 
                                    _ST_* s, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,W,vW,*iflag);

  Eigen::Map< Eigen::Array<_ST_,Eigen::Dynamic,1> > dot(s,V->v.cols());
  dot = (V->v.adjoint() * W->v).diagonal();
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) vV, 
    TYPE(const_sdMat_ptr) vM, 
    _ST_ beta, TYPE(mvec_ptr) vW, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,M,vM,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,W,vW,*iflag);

  W->v = alpha*V->v*M->m + beta*W->v;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA, 
    TYPE(const_sdMat_ptr) vB, 
    _ST_ beta, TYPE(sdMat_ptr) vC, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,B,vB,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,C,vC,*iflag);

  C->m = alpha*A->m*B->m + beta*C->m;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA, 
    TYPE(const_sdMat_ptr) vB, 
    _ST_ beta, TYPE(sdMat_ptr) vC, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,B,vB,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,C,vC,*iflag);

  C->m = alpha*A->m.adjoint()*B->m + beta*C->m;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_times_sdMatT)(_ST_ alpha, TYPE(const_sdMat_ptr) vA,
    TYPE(const_sdMat_ptr) vB,
    _ST_ beta, TYPE(sdMat_ptr) vC, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,B,vB,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,C,vC,*iflag);

  C->m = alpha*A->m*B->m.adjoint() + beta*C->m;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) vV, 
                                       TYPE(const_mvec_ptr) vW, 
                                       _ST_ beta, TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,W,vW,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);

  M->m = alpha*V->v.adjoint()*W->v + beta*M->m;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_QR)(TYPE(mvec_ptr) vV, TYPE(sdMat_ptr) vR, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,R,vR,*iflag);

  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_gather_mvecs)(TYPE(mvec_ptr) V, TYPE(const_mvec_ptr) W[], int nblocks, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_MARK_AS_EXPERIMENTAL(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_scatter_mvecs)(TYPE(const_mvec_ptr) V, TYPE(mvec_ptr) W[], int nblocks, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_MARK_AS_EXPERIMENTAL(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_create_fromRowFuncAndContext)(TYPE(sparseMat_ptr) *vA,
        phist_const_context_ptr vctx,
        phist_lidx maxnne,phist_sparseMat_rowFunc rowFunPtr,void* last_arg,
        int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sparseMat_t*,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(phist::internal::default_context,ctx,vctx,*iflag);
  phist_const_map_ptr map=ctx->row_map;

  phist_lidx nlocal;
  phist_gidx ilower, nglob;
  phist_const_comm_ptr comm;
  PHIST_CHK_IERR(phist_map_get_local_length(map, &nlocal, iflag), *iflag);
  PHIST_CHK_IERR(phist_map_get_ilower(map, &ilower, iflag), *iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(map, &nglob, iflag), *iflag);
  PHIST_CHK_IERR(phist_map_get_comm(map, &comm, iflag), *iflag);

  (*A) = new Traits<_ST_>::sparseMat_t;
  (*A)->map = map;
  (*A)->m.resize(nlocal,nglob);
  (*A)->m.reserve(maxnne*nlocal);

  std::vector<_ST_> vals(maxnne);
  std::vector<phist_gidx> cols(maxnne);

  for (phist_lidx i = 0; i < nlocal; i++)
  {
    phist_gidx row = ilower+i;
    phist_lidx row_nnz = 0;
    PHIST_CHK_IERR( *iflag = rowFunPtr(row,&row_nnz,cols.data(),vals.data(),last_arg), *iflag);
    for(phist_lidx j = 0; j < row_nnz; j++)
      (*A)->m.insert(row,cols[j]) = vals[j];
  }

  *iflag = PHIST_SUCCESS;
}


// NOTE: see the description of sparseMat_read_mm on how we treat input flags for this function
extern "C" void SUBR(sparseMat_create_fromRowFunc)(TYPE(sparseMat_ptr) *A, phist_const_comm_ptr vcomm,
        phist_gidx nrows, phist_gidx ncols, phist_lidx maxnne,
                phist_sparseMat_rowFunc rowFunPtr, void* last_arg,
                int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR( *iflag = (nrows == ncols) ? PHIST_SUCCESS : PHIST_NOT_IMPLEMENTED, *iflag);

  phist_map_ptr map = NULL;
  PHIST_CHK_IERR(phist_map_create(&map, vcomm, nrows, iflag), *iflag);
  phist_context_ptr ctx = NULL;
  PHIST_CHK_IERR(phist_context_create(&ctx,map,NULL,map,map,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFuncAndContext)(A, ctx, maxnne, rowFunPtr, last_arg, iflag), *iflag);

  *iflag = PHIST_SUCCESS;
}

#ifdef IS_COMPLEX
# ifdef IS_DOUBLE
extern "C" void SUBR(mvec_split)(TYPE(const_mvec_ptr) v_V, phist_Dmvec* v_reV, phist_Dmvec* v_imV, int *iflag)
# else
extern "C" void SUBR(mvec_split)(TYPE(const_mvec_ptr) v_V, phist_Smvec* v_reV, phist_Smvec* v_imV, int *iflag)
# endif
{
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,v_V,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_MT_>::mvec_t,reV,v_reV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_MT_>::mvec_t,imV,v_imV,*iflag);

  reV->v = V->v.real();
  imV->v = V->v.imag();
  
  *iflag = PHIST_SUCCESS;
}

# ifdef IS_DOUBLE
extern "C" void SUBR(mvec_combine)(TYPE(mvec_ptr) v_V, phist_Dconst_mvec_ptr v_reV, phist_Dconst_mvec_ptr v_imV, int *iflag)
#else
extern "C" void SUBR(mvec_combine)(TYPE(mvec_ptr) v_V, phist_Sconst_mvec_ptr v_reV, phist_Sconst_mvec_ptr v_imV, int *iflag)
#endif
{
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,v_V,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_MT_>::mvec_t,reV,v_reV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_MT_>::mvec_t,imV,v_imV,*iflag);
  
  V->v = reV->v + imV->v*std::complex<_MT_>(0,1);
  *iflag = PHIST_SUCCESS;
}
#endif
