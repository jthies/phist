#include "phist_config.h"
#ifdef PHIST_USE_SELL
#include "ghost/sell.h"
#endif

// we implement only the double precision real type D
extern "C" void SUBR(type_avail)(int* ierr)
{
  *ierr=0;
}

// \name Matrix input from a file
//@{


//! read a matrix from a MatrixMarket (ASCII) file
extern "C" void SUBR(crsMat_read_mm)(TYPE(crsMat_ptr)* vA, const_comm_ptr_t vcomm,
const char* filename,int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  TOUCH(vA);
  TOUCH(filename);
  *ierr = -99; // not implemented in ghost, use converter script to bin crs
}

//! read a matrix from a Ghost CRS (binary) file.
extern "C" void SUBR(crsMat_read_bin)(TYPE(crsMat_ptr)* vA, const_comm_ptr_t vcomm,
const char* filename,int* ierr)
{
  *ierr=0;
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
PHIST_GHOST_TASK_BEGIN

  ghost_sparsemat_t* mat;
  ghost_context_t *ctx;

  ghost_sparsemat_traits_t *mtraits=new ghost_sparsemat_traits_t;
        *mtraits=(ghost_sparsemat_traits_t)GHOST_SPARSEMAT_TRAITS_INITIALIZER;
#ifdef PHIST_USE_SELL
        mtraits->format = GHOST_SPARSEMAT_SELL;
        ghost_sell_aux_t aux = GHOST_SELL_AUX_INITIALIZER;
        aux.C = PHIST_SELL_C;
        mtraits->sortScope = PHIST_SELL_SIGMA;
        mtraits->aux = &aux;
        mtraits->flags = (ghost_sparsemat_flags_t)(GHOST_SPARSEMAT_DEFAULT|GHOST_SPARSEMAT_PERMUTE);
#else
#warning "GHOST interface compiled to use the reference CRS format, will probably not yield optimal performance"
        mtraits->format = GHOST_SPARSEMAT_CRS;
        mtraits->flags = (ghost_sparsemat_flags_t)(GHOST_SPARSEMAT_DEFAULT);
#endif
        mtraits->datatype = st::ghost_dt;
        char* cfname=const_cast<char*>(filename);
// TODO - check ghost return codes everywhere like this
  PHIST_CHK_GERR(ghost_context_create(&ctx,0,0,
        GHOST_CONTEXT_DEFAULT,cfname,GHOST_SPARSEMAT_SRC_FILE,MPI_COMM_WORLD,1.0),*ierr);
  PHIST_CHK_GERR(ghost_sparsemat_create(&mat,ctx,mtraits,1),*ierr);                               
  PHIST_CHK_GERR(mat->fromFile(mat,cfname),*ierr);
#if PHIST_OUTLEV >= PHIST_VERBOSE
  char *str;
  ghost_context_string(&str,ctx);
  PHIST_SOUT(PHIST_VERBOSE,"%s\n",str);
  free(str); str = NULL;
  ghost_sparsemat_string(&str,mat);
  PHIST_SOUT(PHIST_VERBOSE,"%s\n",str);
  free(str); str = NULL;
#endif
  *vA = (TYPE(crsMat_ptr))mat;
PHIST_GHOST_TASK_END
}

//! read a matrix from a Harwell-Boeing (HB) file
extern "C" void SUBR(crsMat_read_hb)(TYPE(crsMat_ptr)* vA, const_comm_ptr_t vcomm,
const char* filename,int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  TOUCH(vA);
  TOUCH(filename);
  *ierr = -99; // not implemented in ghost, use converter script to bin crs
}

//!@}

//! \name get information about the data distribution in a matrix (maps)

//!@{
//! get the row distribution of the matrix
extern "C" void SUBR(crsMat_get_row_map)(TYPE(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  CAST_PTR_FROM_VOID(const ghost_sparsemat_t,A,vA,*ierr);
  ghost_map_t* map = new ghost_map_t;
  map->ctx = A->context;
  map->vtraits_template=phist_default_vtraits();
  *vmap = (const_map_ptr_t)map;
}

//! get column distribution of a matrix
//! we currently treat all maps as the same as we don't allow any fancy
//! operations using them anyway and ghost can handle both halo'd (colmap)
//! and standard (rowmap) vectors in the mvm.
extern "C" void SUBR(crsMat_get_col_map)(TYPE(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  CAST_PTR_FROM_VOID(const ghost_sparsemat_t,A,vA,*ierr);
  ghost_map_t* map = new ghost_map_t;
  map->ctx = A->context;
  map->vtraits_template=phist_default_vtraits();
  *vmap = (const_map_ptr_t)map;
}

//! get the map for vectors x in y=A*x
//! we currently treat all maps as the same as we don't allow any fancy
//! operations using them anyway and ghost can handle both halo'd (colmap)
//! and standard (rowmap) vectors in the mvm.
extern "C" void SUBR(crsMat_get_domain_map)(TYPE(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  SUBR(crsMat_get_col_map)(vA,vmap,ierr);
}

//! get the map for vectors y in y=A*x
//! we currently treat all maps as the same as we don't allow any fancy
//! operations using them anyway and ghost can handle both halo'd (colmap)
//! and standard (rowmap) vectors in the mvm.
extern "C" void SUBR(crsMat_get_range_map)(TYPE(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  SUBR(crsMat_get_row_map)(vA,vmap,ierr);
}
//@}

//! \name constructors

//@{
//! create a block-vector. The entries are stored contiguously
//! at val in column major ordering.
extern "C" void SUBR(mvec_create)(TYPE(mvec_ptr)* vV, 
        const_map_ptr_t vmap, int nvec, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
PHIST_GHOST_TASK_BEGIN
  *ierr=0;
  CAST_PTR_FROM_VOID(const ghost_map_t, map,vmap,*ierr);
  ghost_densemat_t* result;
  ghost_densemat_traits_t vtraits = map->vtraits_template;/*ghost_cloneVtraits(map->vtraits_template);*/
        vtraits.ncols=nvec;
        vtraits.datatype = st::ghost_dt;
        vtraits.flags = (ghost_densemat_flags_t)(vtraits.flags & ~GHOST_DENSEMAT_VIEW);
  PHIST_CHK_GERR(ghost_densemat_create(&result,map->ctx,vtraits),*ierr);
  ST zero = st::zero();
  // this allocates the vector and fills it with zeros
  PHIST_CHK_GERR(result->fromScalar(result,&zero),*ierr);
  PHIST_DEB("mvec nrows: %" PRlidx "\n",result->traits.nrows);
  *vV=(TYPE(mvec_ptr))(result);
PHIST_GHOST_TASK_END
}

//! create a block-vector as view of raw data. The map tells the object
//! how many rows it should 'see' in the data (at most lda, the leading
//! dimension of the 2D array values). CAVEAT: This function only works
//! if nrowshalo==nrows in the map, which is in general only the case for
//! if there is only 1 MPI process or the matrix is trivially parallel.
extern "C" void SUBR(mvec_create_view)(TYPE(mvec_ptr)* vV, const_map_ptr_t vmap, 
        _ST_* values, lidx_t lda, int nvec,
        int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);

  CAST_PTR_FROM_VOID(const ghost_map_t, map,vmap,*ierr);
  ghost_densemat_t* result;
  ghost_densemat_traits_t vtraits = map->vtraits_template;/*ghost_cloneVtraits(map->vtraits_template);*/
        vtraits.ncols=nvec;
        vtraits.datatype = st::ghost_dt;

  PHIST_CHK_GERR(ghost_densemat_create(&result,map->ctx,vtraits),*ierr);

  if ((result->traits.nrows!=result->traits.nrowshalo)||(result->traits.nrowshalo!=lda))
  {
    PHIST_OUT(PHIST_ERROR,"viewing plain data as ghost_vec only works "
                          "for node-local or trivially parallel matrices!\n");
    PHIST_OUT(PHIST_ERROR,"nrows=%" PRlidx ", nrowshalo=%" PRlidx ", lda=%" PRlidx "\n",
        vtraits.nrows,vtraits.nrowshalo,lda);
    *ierr=-1;
    return;
  }

  PHIST_CHK_GERR(result->viewPlain(result,(void*)values,0,0,lda),*ierr);
  *vV=(TYPE(mvec_ptr))(result);
  return;
}


//! create a serial dense n x m matrix on all procs, with column major
//! ordering.
extern "C" void SUBR(sdMat_create)(TYPE(sdMat_ptr)* vM, int nrows, int ncols, 
        const_comm_ptr_t vcomm, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  CAST_PTR_FROM_VOID(MPI_Comm,comm,vcomm,*ierr);
#include "phist_std_typedefs.hpp"
PHIST_GHOST_TASK_BEGIN
  *ierr=0;
  ghost_densemat_t* result;
  ghost_densemat_traits_t dmtraits=GHOST_DENSEMAT_TRAITS_INITIALIZER;
        dmtraits.nrows=nrows;
        dmtraits.nrowshalo=nrows;
        dmtraits.nrowspadded=nrows;
        dmtraits.ncols=ncols;
        dmtraits.ncolspadded=ncols;
        dmtraits.datatype=st::ghost_dt;
#ifdef PHIST_SDMATS_ROW_MAJOR
        dmtraits.storage=GHOST_DENSEMAT_ROWMAJOR;
#else
        dmtraits.storage=GHOST_DENSEMAT_COLMAJOR;
#endif
  // I think the sdMat should not have a context
  ghost_context_t* ctx=NULL;
  PHIST_CHK_GERR(ghost_context_create(&ctx,nrows, ncols, GHOST_CONTEXT_REDUNDANT, 
        NULL, GHOST_SPARSEMAT_SRC_NONE, *comm, 1.0),*ierr);
  ghost_densemat_create(&result,ctx,dmtraits);
  ST zero = st::zero();
  result->fromScalar(result,&zero);
  *vM=(TYPE(sdMat_ptr))result;
PHIST_GHOST_TASK_END
}

//@}

//! retrieve local length of the vectors in V
extern "C" void SUBR(mvec_my_length)(TYPE(const_mvec_ptr) vV, lidx_t* len, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;
  CAST_PTR_FROM_VOID(const ghost_densemat_t,V,vV,*ierr);
  PHIST_CHK_IERR(*ierr=check_local_size(V->traits.nrows),*ierr);
  *len = V->traits.nrows;
}

//! retrieve the map of the vectors in V
extern "C" void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) vV, const_map_ptr_t* vmap, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  CAST_PTR_FROM_VOID(const ghost_densemat_t,V,vV,*ierr);
  ghost_map_t* map = new ghost_map_t;
  map->ctx=V->context; 
  map->vtraits_template=V->traits;
  *vmap=(const_map_ptr_t)map;
}

//! retrieve the comm used for MPI communication in V
extern "C" void SUBR(mvec_get_comm)(TYPE(const_mvec_ptr) vV, const_comm_ptr_t* vcomm, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  CAST_PTR_FROM_VOID(const ghost_densemat_t,V,vV,*ierr);

  if (V->context!=NULL)
  {
    *vcomm=(const_comm_ptr_t)&(V->context->mpicomm);
  }
  else
  {
    PHIST_OUT(PHIST_WARNING,"in mvec_get_comm: ghost_densemat_t without context!\n");
    MPI_Comm* comm = new MPI_Comm;
    *comm=MPI_COMM_SELF;
    *vcomm=(const_comm_ptr_t)(comm);
  }
  return;
}

//! retrieve number of vectors/columns in V
extern "C" void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) vV, int* nvec, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;
  CAST_PTR_FROM_VOID(const ghost_densemat_t,V,vV,*ierr);
  PHIST_CHK_IERR(*ierr=check_local_size(V->traits.ncols),*ierr);
  *nvec = V->traits.ncols;
}

//! get number of rows in local dense matrix
extern "C" void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) vM, int* nrows, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  CAST_PTR_FROM_VOID(const ghost_densemat_t,M,vM,*ierr);
  PHIST_CHK_IERR(*ierr=check_local_size(M->traits.nrows),*ierr);
  *nrows = M->traits.nrows;
}
  
//! get number of cols in local dense matrix
extern "C" void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) vM, int* ncols, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  CAST_PTR_FROM_VOID(const ghost_densemat_t,M,vM,*ierr);
  *ncols = M->traits.ncols;
}


extern "C" void SUBR(mvec_extract_view)(TYPE(mvec_ptr) vV, _ST_** val, lidx_t* lda, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"

  CAST_PTR_FROM_VOID(ghost_densemat_t,V, vV, *ierr);
  if (V->traits.flags & GHOST_DENSEMAT_SCATTERED)
  {
    PHIST_OUT(PHIST_ERROR,"%s: cannot view data with non-constant stride using "
        "this function (file %s, line %d)\n", __FUNCTION__, __FILE__, __LINE__);
    *ierr=-1;
    return;
  }
  if (V->val==NULL)
  {
    PHIST_OUT(PHIST_ERROR,"%s, pointer is NULL\n",__FUNCTION__);
    *ierr=-2;
    return;
  }
  PHIST_CHK_GERR(ghost_densemat_valptr(V,(void**)val),*ierr);
  PHIST_CHK_IERR(*ierr=check_local_size(V->traits.nrowspadded),*ierr);

#ifdef PHIST_MVECS_ROW_MAJOR
  *lda = V->traits.ncolspadded;
#else
  *lda = V->traits.nrowspadded;
#endif
}

extern "C" void SUBR(mvec_to_device)(TYPE(mvec_ptr) vV, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  CAST_PTR_FROM_VOID(ghost_densemat_t,V, vV, *ierr);
  PHIST_SOUT(PHIST_DEBUG,"ghost densemat upload\n"
                         "nrows=%" PRlidx ", ncols=%" PRlidx "\n"
                         "nrowshalo=%" PRlidx "\n"
                         "nrowspadded=%" PRlidx ", ncolspadded=%" PRlidx "\n",
                         V->traits.nrows, V->traits.ncols, 
                         V->traits.nrowshalo,
                         V->traits.nrowspadded, V->traits.ncolspadded);
  PHIST_SOUT(PHIST_DEBUG,"V flags: %d\n",(int)V->traits.flags);
  PHIST_CHK_GERR(V->upload(V),*ierr);
}

extern "C" void SUBR(mvec_from_device)(TYPE(mvec_ptr) vV, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  CAST_PTR_FROM_VOID(ghost_densemat_t,V, vV, *ierr);
  PHIST_CHK_GERR(V->download(V),*ierr);
}

extern "C" void SUBR(sdMat_to_device)(TYPE(sdMat_ptr) vM, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  CAST_PTR_FROM_VOID(ghost_densemat_t,M, vM, *ierr);
  PHIST_CHK_GERR(M->upload(M),*ierr);
}

extern "C" void SUBR(sdMat_from_device)(TYPE(sdMat_ptr) vM, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  CAST_PTR_FROM_VOID(ghost_densemat_t,M, vM, *ierr);
  PHIST_CHK_GERR(M->download(M),*ierr);
}

extern "C" void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) vM, _ST_** val, lidx_t* lda, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  CAST_PTR_FROM_VOID(ghost_densemat_t,M, vM, *ierr);

  if (M->traits.flags & GHOST_DENSEMAT_SCATTERED)
  {
    PHIST_OUT(PHIST_ERROR,"%s: cannot view data with non-constant stride using "
        "this function (file %s, line %d)\n", __FUNCTION__, __FILE__, __LINE__);
    *ierr=-1; 
    return;
  }

  PHIST_CHK_GERR(ghost_densemat_valptr(M,(void**)val),*ierr);

  PHIST_CHK_IERR(*ierr=check_local_size(M->traits.nrowspadded),*ierr);

#ifdef PHIST_SDMATS_ROW_MAJOR
  *lda = M->traits.ncolspadded;
#else
  *lda = M->traits.nrowspadded;
#endif

/*
  if (M->traits.storage==GHOST_DENSEMAT_ROWMAJOR)
  {
    *lda = M->traits.nrowspadded;
  }
else
  {

    
    *lda = M->traits.ncolspadded;
  }*/
}

extern "C" void SUBR(mvec_to_mvec)(TYPE(const_mvec_ptr) v_in, TYPE(mvec_ptr) v_out, int* ierr)
{
  CAST_PTR_FROM_VOID(ghost_densemat_t,V_in,v_in,*ierr);
  CAST_PTR_FROM_VOID(ghost_densemat_t,V_out,v_out,*ierr);
  
  if (V_in->context!=V_out->context)
  {
    PHIST_SOUT(PHIST_WARNING,"function %s only implemented for simple permutation operations\n"
                             "where the result and input vectors have the same context and  \n"
                             "one of them may be permuted\n",__FUNCTION__);
    *ierr=PHIST_NOT_IMPLEMENTED;
    return;
  }
  
  bool resultPermuted=false;//V_out->flags && IS_PERMUTED;
  bool inputPermuted=false;//V_out->flags && IS_PERMUTED;
  
  // first copy the data
  PHIST_CHK_GERR(V_out->fromVec(V_out,V_in,0,0),*ierr);
  // check permutation state
  if (resultPermuted==inputPermuted) return;
  if (resultPermuted)
  {
    PHIST_CHK_GERR(V_out->permute(V_out,V_out->context->permutation,GHOST_PERMUTATION_ORIG2PERM),*ierr);
  }
  else
  {
    PHIST_CHK_GERR(V_out->permute(V_out,V_out->context->permutation,GHOST_PERMUTATION_PERM2ORIG),*ierr);
  }
  return;
}

//! get a new vector that is a view of some columns of the original one,
//! Vblock = V(:,jmin:jmax). The new object Vblock is created but does not
//! allocate memory for the vector entries, instead using the entries from V
//! directly. When mvec_delete(Vblock) is called, the library has to take care
//! that the value array is not deleted 
extern "C" void SUBR(mvec_view_block)(TYPE(mvec_ptr) vV,
                             TYPE(mvec_ptr)* vVblock,
                             int jmin, int jmax, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*ierr);
  ghost_densemat_t *Vblock;
  V->viewCols(V, &Vblock, jmax-jmin+1, jmin);

  if (*vVblock!=NULL)
  {
    CAST_PTR_FROM_VOID(ghost_densemat_t,tmp,*vVblock,*ierr);
    //PHIST_DEB("destroying previous vector (view)\n");
    tmp->destroy(tmp);
  }
  PHIST_CHK_IERR(*ierr=((Vblock->traits.flags&GHOST_DENSEMAT_VIEW)-GHOST_DENSEMAT_VIEW),*ierr);
  *vVblock = (TYPE(mvec_ptr))Vblock;
}

//! get a new vector that is a copy of some columns of the original one,  
//! Vblock = V(:,jmin:jmax). The object Vblock must be created beforehand 
//! and the corresponding columns of V are copied into the value array    
//! of Vblock. V is not modified.
extern "C" void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) vV,
                             TYPE(mvec_ptr) vVblock,
                             int jmin, int jmax, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
PHIST_GHOST_TASK_BEGIN
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*ierr);
  CAST_PTR_FROM_VOID(ghost_densemat_t,Vblock,vVblock,*ierr);
  *ierr=0;
#ifdef TESTING
// nonzero error code if #vectors in Vblock too small or large
  PHIST_CHK_IERR(*ierr=(jmax-jmin+1)-Vblock->traits.ncols,*ierr);
#else
  TOUCH(jmax);
#endif  
  //TODO check bounds of Vblock
  Vblock->fromVec(Vblock,V,0,jmin);
PHIST_GHOST_TASK_END
}

//! given a multi-vector Vblock, set V(:,jmin:jmax)=Vblock by copying the corresponding
//! vectors. Vblock is not modified.
extern "C" void SUBR(mvec_set_block)(TYPE(mvec_ptr) vV,
                             TYPE(const_mvec_ptr) vVblock,
                             int jmin, int jmax, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
PHIST_GHOST_TASK_BEGIN
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*ierr);
  CAST_PTR_FROM_VOID(ghost_densemat_t,Vblock,vVblock,*ierr);
  // TODO - bounds checking
  // create a view of the requested columns of V
  ghost_densemat_t *Vcols;
  V->viewCols(V,&Vcols,jmax-jmin+1,jmin);
  // copy the data
  Vcols->fromVec(Vcols,Vblock,0,0);
  // delete the view
  Vcols->destroy(Vcols);
PHIST_GHOST_TASK_END
}

//! get a new sdMat that is a view of some rows and columns of the original one,
//! Mblock = M(imin:imax,jmin:jmax). The new object Vblock is created but does not
//! allocate memory for the vector entries, instead using the entries from V
//! directly.
extern "C" void SUBR(sdMat_view_block)(TYPE(mvec_ptr) vM, TYPE(mvec_ptr)* vMblock,
                             int imin, int imax, int jmin, int jmax, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,M,vM,*ierr);

  //TODO: we only view the host side of the vector here, this function should
  //      eventually be moved into ghost and the accelerator stuff added.

  // first just create a view of the corresponding columns
  ghost_densemat_t *Mblock;
  M->viewVec(M, &Mblock, imax-imin+1,imin,jmax-jmin+1, jmin);

  if (*vMblock!=NULL)
  {
    //PHIST_DEB("deleting previous object in %s\n",__FUNCTION__);
    CAST_PTR_FROM_VOID(ghost_densemat_t,tmp,*vMblock,*ierr);
    tmp->destroy(tmp);
  }
  PHIST_CHK_IERR(*ierr=((Mblock->traits.flags&GHOST_DENSEMAT_VIEW)-GHOST_DENSEMAT_VIEW),*ierr);
  *vMblock = (TYPE(sdMat_ptr))Mblock;
}

//! get a new matrix that is a copy of some rows and columns of the original one,  
//! Mblock = M(imin:imax,jmin:jmax). The object Mblock must be created beforehand 
//! and the corresponding columns of M are copied into the value array    
//! of Mblock. M is not modified.
extern "C" void SUBR(sdMat_get_block)(TYPE(const_sdMat_ptr) vM,
                             TYPE(sdMat_ptr) vMblock,
                             int imin, int imax, int jmin, int jmax, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
PHIST_GHOST_TASK_BEGIN
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,M,vM,*ierr);
  CAST_PTR_FROM_VOID(ghost_densemat_t,Mblock,vMblock,*ierr);
 
  /*
  ghost_densemat_t* Mb_view=NULL;
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(vM,(TYPE(sdMat_ptr)*)&Mb_view,imin,imax,jmin,jmax,ierr),*ierr);
  Mb_view->fromVec(Mb_view,Mblock,0,0);
  Mb_view->destroy(Mb_view);
  */ 
  
  M->fromVec(Mblock,M,imin,jmin);
PHIST_GHOST_TASK_END
}

//! given a serial dense matrix Mblock, set M(imin:imax,jmin:jmax)=Mblock by 
//! copying the corresponding elements. Mblock is not modified.
extern "C" void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) vM, 
                             TYPE(const_sdMat_ptr) vMblock,
                             int imin, int imax, int jmin, int jmax, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
PHIST_GHOST_TASK_BEGIN
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,Mblock,vMblock,*ierr);

  ghost_densemat_t* Mb_view=NULL;
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(vM,(TYPE(sdMat_ptr)*)&Mb_view,imin,imax,jmin,jmax,ierr),*ierr);
  Mb_view->fromVec(Mb_view,Mblock,0,0);
  Mb_view->destroy(Mb_view);
PHIST_GHOST_TASK_END
}

//! \name destructors

//@{

//!
extern "C" void SUBR(crsMat_delete)(TYPE(crsMat_ptr) vA, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_sparsemat_t,A,vA,*ierr);
  A->destroy(A);
}

//!
extern "C" void SUBR(mvec_delete)(TYPE(mvec_ptr) vV, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  if (vV==NULL) return;
  CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*ierr);
  V->destroy(V);
}

//!
extern "C" void SUBR(sdMat_delete)(TYPE(sdMat_ptr) vM, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  if (vM==NULL) return;
  CAST_PTR_FROM_VOID(ghost_densemat_t,M,vM,*ierr);
  M->destroy(M);
}

//@}

//! \name Numerical functions
//!@{

//! put scalar value into all elements of a multi-vector
extern "C" void SUBR(mvec_put_value)(TYPE(mvec_ptr) vV, _ST_ value, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
PHIST_GHOST_TASK_BEGIN
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*ierr);
  PHIST_DEB("put value, V @ %p. V->traits.nrows=%" PRlidx "\n",V,V->traits.nrows);
  V->fromScalar(V,(void*)&value);
PHIST_GHOST_TASK_END
}

//! put scalar value into all elements of a multi-vector
extern "C" void SUBR(sdMat_put_value)(TYPE(sdMat_ptr) vV, _ST_ value, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
PHIST_GHOST_TASK_BEGIN
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*ierr);
  V->fromScalar(V,(void*)&value);
PHIST_GHOST_TASK_END
}

//! put random numbers into all elements of a multi-vector
extern "C" void SUBR(mvec_random)(TYPE(mvec_ptr) vV, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
PHIST_GHOST_TASK_BEGIN
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*ierr);
  V->fromRand(V);
PHIST_GHOST_TASK_END
}

extern "C" void SUBR(mvec_print)(TYPE(const_mvec_ptr) vV, int* ierr)
{
  *ierr = 0;
  ENTER_FCN(__FUNCTION__);
  CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*ierr);
  std::cout << "# local rows: "<<V->traits.nrows<<std::endl;
  std::cout << "# vectors:    "<<V->traits.ncols<<std::endl;
  char *str;
  V->string(V,&str);
  std::cout << str <<std::endl;
  free(str); str = NULL;
}

extern "C" void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) vM, int* ierr)
{
  *ierr=0;
  ENTER_FCN(__FUNCTION__);
  CAST_PTR_FROM_VOID(ghost_densemat_t,M,vM,*ierr);
  std::cout << "# rows: "<<M->traits.nrows<<std::endl;
  std::cout << "# cols: "<<M->traits.ncols<<std::endl;
  char *str;
  M->string(M,&str);
  std::cout << str <<std::endl;
  free(str); str = NULL;
}


//! put random numbers into all elements of a serial dense matrix
extern "C" void SUBR(sdMat_random)(TYPE(sdMat_ptr) vM, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
PHIST_GHOST_TASK_BEGIN
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,M,vM,*ierr);
  M->fromRand(M);
PHIST_GHOST_TASK_END
}

//! \name Numerical functions

//! compute the 2-norm) of each column of v                   
//! (vnrm[i] must be pre-allocated by caller)
  extern "C" void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) vV,
                            _MT_* vnrm, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"  
PHIST_GHOST_TASK_BEGIN
  int i,nv;
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*ierr);  
  nv=V->traits.ncols;
  _ST_ tmp[V->traits.ncols];
  ghost_dot(tmp,V,V);
  for (i=0;i<nv;i++) vnrm[i]=mt::sqrt(st::real(tmp[i]));
  return;
PHIST_GHOST_TASK_END
}

  //! normalize (in the 2-norm) each column of v and return ||v||_2
  //! for each vector i in vnrm[i] (must be pre-allocated by caller)
  extern "C" void SUBR(mvec_normalize)(TYPE(mvec_ptr) vV,
                            _MT_* vnrm, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"  
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*ierr);  
  // TODO - this call doesn't return the norm as we wish
  //V->normalize(V);    
  PHIST_CHK_IERR(SUBR(mvec_norm2)(vV,vnrm,ierr),*ierr);
  _ST_ inrm[V->traits.ncols];
  for (int i=0;i<V->traits.ncols;i++) inrm[i]=st::one()/vnrm[i];
  PHIST_CHK_IERR(SUBR(mvec_vscale)(vV,inrm,ierr),*ierr);
  return;
}

//! scale each column i of v and by scalar[i]
extern "C" void SUBR(mvec_scale)(TYPE(mvec_ptr) vV, 
                            _ST_ scalar, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
PHIST_GHOST_TASK_BEGIN
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*ierr);  
  V->scale(V,(void*)&scalar);
  return;
PHIST_GHOST_TASK_END
}

//! scale each column i of v and by scalar[i]
extern "C" void SUBR(mvec_vscale)(TYPE(mvec_ptr) vV, 
                            const _ST_* scalar, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"  
PHIST_GHOST_TASK_BEGIN
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*ierr);  
  V->vscale(V,(void*)scalar);
  return;
PHIST_GHOST_TASK_END
}

//! y=alpha*x+beta*y
extern "C" void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) vX,
                            _ST_ beta,  TYPE(mvec_ptr)       vY, 
                            int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
PHIST_GHOST_TASK_BEGIN
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,X,vX,*ierr);
  CAST_PTR_FROM_VOID(ghost_densemat_t,Y,vY,*ierr);
  if (beta==st::one())
  {
    Y->axpy(Y,X,(void*)&alpha);
  }
  else
  {
    Y->axpby(Y,X,(void*)&alpha,(void*)&beta);
  }
PHIST_GHOST_TASK_END
}

//! y[j]=alpha[j]*x[j]+beta[j]*y[j] for all columns j
extern "C" void SUBR(mvec_vadd_mvec)(const _ST_ *alpha, TYPE(const_mvec_ptr) vX,
                            _ST_ beta,  TYPE(mvec_ptr)       vY, 
                            int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
PHIST_GHOST_TASK_BEGIN
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,X,vX,*ierr);
  CAST_PTR_FROM_VOID(ghost_densemat_t,Y,vY,*ierr);
  if(beta == st::one())
  {
    Y->vaxpy(Y,X,(void*)alpha);
  }
  else
  {
    // ghost also expects a vector for beta, so construct one:
    int nvec = 0;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vY,&nvec,ierr),*ierr);
    _ST_ *b = (_ST_*)malloc(nvec*sizeof(_ST_));
    for(int i = 0; i < nvec; i++)
      b[i] = beta;

    Y->vaxpby(Y,X,(void*)alpha,(void*)b);
    free(b);
  }
PHIST_GHOST_TASK_END
}

//! B=alpha*A+beta*B
extern "C" void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA,
                            _ST_ beta,  TYPE(sdMat_ptr)       vB,
                            int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(alpha,vA,beta,vB,ierr), *ierr);
}

//! y=alpha*A*x+beta*y.
extern "C" void SUBR(crsMat_times_mvec)(_ST_ alpha, TYPE(const_crsMat_ptr) vA, TYPE(const_mvec_ptr) vx, 
_ST_ beta, TYPE(mvec_ptr) vy, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr=0;

#ifdef PHIST_TIMEMONITOR
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vx, &nvec, ierr), *ierr);
  for(int i = 0; i < nvec; i++)
    phist_totalMatVecCount();
#endif

  CAST_PTR_FROM_VOID(ghost_sparsemat_t,A,vA,*ierr);
  CAST_PTR_FROM_VOID(ghost_densemat_t,x,vx,*ierr);
  CAST_PTR_FROM_VOID(ghost_densemat_t,y,vy,*ierr);
  if (alpha==st::zero())
  {
    // no MVM needed
    if (beta==st::zero())
    {
      PHIST_CHK_IERR(SUBR(mvec_put_value)(vy,beta,ierr),*ierr);
    }
    else if (beta!=st::one())
    {
PHIST_GHOST_TASK_BEGIN
      y->scale(y,(void*)&beta);
PHIST_GHOST_TASK_END
    }
  }
  else
  {
PHIST_GHOST_TASK_BEGIN
    ghost_spmv_flags_t spMVM_opts=GHOST_SPMV_DEFAULT;
    // currently the vector mode is the only one working with MPI and multiple RHS
    spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_MODE_VECTOR);
    //void* old_scale = A->traits->scale;
    if (beta==st::one())
    {
      spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_AXPY);
    }
    else if (beta!=st::zero())
    {
      spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_AXPBY);
    }
    if (alpha!=st::one())
    {
      spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_SCALE);
      *ierr=ghost_spmv(y,A,x,&spMVM_opts,&alpha,&beta);
    }
    else
    {
      *ierr=ghost_spmv(y,A,x,&spMVM_opts,&beta);
    }
PHIST_GHOST_TASK_END
  }
}

//! y=alpha*A*x+beta*y.
extern "C" void SUBR(crsMatT_times_mvec)(_ST_ alpha, TYPE(const_crsMat_ptr) vA, TYPE(const_mvec_ptr) vx, 
_ST_ beta, TYPE(mvec_ptr) vy, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr=-99;
  return;
}
//! y[i]=alpha*(A*x[i]+shifts[i]*x[i]) + beta*y[i]
extern "C" void SUBR(crsMat_times_mvec_vadd_mvec)(_ST_ alpha, TYPE(const_crsMat_ptr) vA,
        const _ST_ shifts[], TYPE(const_mvec_ptr) vx, _ST_ beta, TYPE(mvec_ptr) vy, int* 
        ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr=0;

#ifdef PHIST_TIMEMONITOR
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vx, &nvec, ierr), *ierr);
  for(int i = 0; i < nvec; i++)
    phist_totalMatVecCount();
#endif

  CAST_PTR_FROM_VOID(ghost_sparsemat_t,A,vA,*ierr);
  CAST_PTR_FROM_VOID(ghost_densemat_t,x,vx,*ierr);
  CAST_PTR_FROM_VOID(ghost_densemat_t,y,vy,*ierr);
  if (alpha==st::zero())
  {
    PHIST_CHK_IERR(SUBR(crsMat_times_mvec)(alpha,vA,vx,beta,vy,ierr),*ierr);
  }
  else
  {
PHIST_GHOST_TASK_BEGIN
    int nvec;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vx, &nvec, ierr), *ierr);

    ghost_spmv_flags_t spMVM_opts=GHOST_SPMV_DEFAULT;
    // currently the vector mode is the only one working with MPI and multiple RHS
    spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_MODE_VECTOR);
    spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_VSHIFT);
    ST ghost_shifts[nvec];
    for (int i=0;i<nvec;i++) ghost_shifts[i]=-shifts[i];
    
    if (beta==st::one())
    {
      spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_AXPY);
      if (alpha!=st::one())
      {
        spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_SCALE);
        *ierr=ghost_spmv(y,A,x,&spMVM_opts,&alpha,ghost_shifts);
      }
      else
      {
        *ierr=ghost_spmv(y,A,x,&spMVM_opts,ghost_shifts);
      }
    }
    else if (beta!=st::zero())
    {
      spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_AXPBY);
      if (alpha!=st::one())
      {
        spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_SCALE);
        *ierr=ghost_spmv(y,A,x,&spMVM_opts,&alpha,&beta,ghost_shifts);
      }
      else
      {
        *ierr=ghost_spmv(y,A,x,&spMVM_opts,&beta,ghost_shifts);
      }
    }
    else
    {
      if (alpha!=st::one())
      {
        spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_SCALE);
        *ierr=ghost_spmv(y,A,x,&spMVM_opts,&alpha,ghost_shifts);
      }
      else
      {
        *ierr=ghost_spmv(y,A,x,&spMVM_opts,ghost_shifts);
      }
    }
PHIST_GHOST_TASK_END
  }
}

//! dot product of vectors v_i and w_i, i=1..numvecs
extern "C" void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) vV, TYPE(const_mvec_ptr) vW, _ST_* s, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
PHIST_GHOST_TASK_BEGIN
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*ierr);
  CAST_PTR_FROM_VOID(ghost_densemat_t,W,vW,*ierr);
  ghost_dot(s,V,W);
PHIST_GHOST_TASK_END
  }

//! dense tall skinny matrix-matrix product yielding a serial dense matrix
//! C=alpha*V'*W+beta*C. C is replicated on all MPI processes sharing V and W.
extern "C" void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) vV, TYPE(const_mvec_ptr) vW, _ST_ beta, TYPE(sdMat_ptr) vC, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
PHIST_GHOST_TASK_BEGIN
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*ierr);
  CAST_PTR_FROM_VOID(ghost_densemat_t,W,vW,*ierr);
  CAST_PTR_FROM_VOID(ghost_densemat_t,C,vC,*ierr);
#ifdef IS_COMPLEX
  char trans[]="C";
#else
  char trans[]="T";
#endif  

    lidx_t ncC = C->traits.ncols;

  PHIST_DEB("VtV=C, V %" PRlidx "x%" PRlidx ", \n"
            "       W %" PRlidx "x%" PRlidx ", \n"
            "       C %" PRlidx "x%" PRlidx "\n", 
  V->traits.nrows,V->traits.ncols,
  W->traits.nrows,W->traits.ncols,
  C->traits.nrows,C->traits.ncols);

  PHIST_CHK_GERR(ghost_gemm(C,V,trans,W,(char*)"N",(void*)&alpha,(void*)&beta,GHOST_GEMM_ALL_REDUCE),*ierr);
PHIST_GHOST_TASK_END
  }


//! n x m multi-vector times m x k dense matrix gives n x k multi-vector,
//! W=alpha*V*C + beta*W
extern "C" void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) vV,
                                       TYPE(const_sdMat_ptr) vC,
                           _ST_ beta,  TYPE(mvec_ptr) vW,
                                       int* ierr)
  {
#include "phist_std_typedefs.hpp"
    ENTER_FCN(__FUNCTION__);
PHIST_GHOST_TASK_BEGIN
    *ierr=0;
    CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*ierr);
    CAST_PTR_FROM_VOID(ghost_densemat_t,C,vC,*ierr);
    CAST_PTR_FROM_VOID(ghost_densemat_t,W,vW,*ierr);

    lidx_t nrV,nrW;
    int ncV, ncW, nrC, ncC;
    nrV=V->traits.nrows;  ncV=V->traits.ncols;
    nrW=W->traits.nrows;  ncW=V->traits.ncols;
    nrC=C->traits.nrows;  ncC=V->traits.ncols;

#ifdef TESTING
    PHIST_CHK_IERR(*ierr=nrV-nrW,*ierr);
    PHIST_CHK_IERR(*ierr=nrC-ncV,*ierr);
    PHIST_CHK_IERR(*ierr=ncC-ncW,*ierr);
    //PHIST_DEB("V'C with V %" PRlidx "x%d, C %dx%d and result %" PRlidx "x%d\n", nrV,ncV,nrC,ncC,nrW,ncW);
#endif
    // note: C is replicated, so this operation is a purely local one.
    PHIST_CHK_GERR(ghost_gemm(W,V,(char*)"N",C,(char*)"N",(void*)&alpha,(void*)&beta,GHOST_GEMM_NO_REDUCE),*ierr);
PHIST_GHOST_TASK_END
  }

//! C <- V*C
extern "C" void SUBR(mvec_times_sdMat_inplace)(TYPE(mvec_ptr) vV,
                                       TYPE(const_sdMat_ptr) vC,
                                       int* ierr)
  {
#include "phist_std_typedefs.hpp"
    ENTER_FCN(__FUNCTION__);
PHIST_GHOST_TASK_BEGIN
    *ierr=0;
    CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*ierr);
    CAST_PTR_FROM_VOID(ghost_densemat_t,C,vC,*ierr);

#ifdef TESTING
    int ncV, nrC, ncC;
    ncV=V->traits.ncols;
     nrC=C->traits.nrows;  ncC=V->traits.ncols;

    PHIST_CHK_IERR(*ierr=nrC-ncV,*ierr);
    PHIST_CHK_IERR(*ierr=nrC-ncC,*ierr);
    //PHIST_DEB("V'C with V %" PRlidx "x%d, C %dx%d and result %" PRlidx "x%d\n", nrV,ncV,nrC,ncC,nrW,ncW);
#endif
    // note: C is replicated, so this operation is a purely local one.
    ST alpha=st::one();
    ST beta=st::zero();
    // ghost internally picks the in-place variant of possible
    PHIST_CHK_GERR(ghost_gemm(V,V,(char*)"N",C,(char*)"N",(void*)&alpha,(void*)&beta,GHOST_GEMM_NO_REDUCE),*ierr);
PHIST_GHOST_TASK_END
  }
//! n x m serial dense matrix times m x k serial dense matrix gives n x k sdMat,
//! C=alpha*V*W + beta*C (serial XGEMM wrapper)
extern "C" void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vV,
                                         TYPE(const_sdMat_ptr) vW,
                              _ST_ beta, TYPE(sdMat_ptr) vC,
                                         int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
PHIST_GHOST_TASK_BEGIN
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*ierr);
  CAST_PTR_FROM_VOID(ghost_densemat_t,W,vW,*ierr);
  CAST_PTR_FROM_VOID(ghost_densemat_t,C,vC,*ierr);
  char trans[]="N";  
  PHIST_CHK_GERR(ghost_gemm(C,V,trans,W,trans,(void*)&alpha,(void*)&beta,GHOST_GEMM_NO_REDUCE),*ierr);
PHIST_GHOST_TASK_END
  }

//! n x m conj. transposed serial dense matrix times m x k serial dense matrix gives m x k sdMat,
//! C=alpha*V*W + beta*C (serial XGEMM wrapper)
extern "C" void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vV,
                                         TYPE(const_sdMat_ptr) vW,
                              _ST_ beta, TYPE(sdMat_ptr) vC,
                                         int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
PHIST_GHOST_TASK_BEGIN
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*ierr);
  CAST_PTR_FROM_VOID(ghost_densemat_t,W,vW,*ierr);
  CAST_PTR_FROM_VOID(ghost_densemat_t,C,vC,*ierr);
#ifdef IS_COMPLEX
  char trans[]="C";
#else
  char trans[]="T";
#endif  
  PHIST_CHK_GERR(ghost_gemm(C, V, trans,W, (char*)"N", (void*)&alpha, (void*)&beta, GHOST_GEMM_NO_REDUCE),*ierr);
PHIST_GHOST_TASK_END
  }


//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.   
//! Q is computed in place of V. If V does not have full rank, ierr>0   
//! indicates the dimension of the null-space of V. The first m-ierr    
//! columns of Q are an orthogonal basis of the column space of V, the  
//! remaining columns form a basis for the null space.  
extern "C" void SUBR(mvec_QR)(TYPE(mvec_ptr) vV, TYPE(sdMat_ptr) vR, int* ierr)
  {
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*ierr);
  CAST_PTR_FROM_VOID(ghost_densemat_t,R,vR,*ierr);

  int rank;
  MT rankTol=1000*mt::eps();
  int ncols=V->traits.ncols;
  if (ncols==1)
  {
    // we need a special treatment here because TSQR
    // uses a relative tolerance to determine rank deficiency,
    // so a single zero vector is not detected to be rank deficient.
    PHIST_DEB("mvec_QR: single-vector case\n");
    MT nrm;
    PHIST_CHK_IERR(SUBR(mvec_normalize)(vV,&nrm,ierr),*ierr);
    PHIST_DEB("single vector QR, R=%8.4e\n",nrm);
    rank=1;
    if (nrm<rankTol)
    {
      PHIST_DEB("zero vector detected\n");
      // randomize the vector
      PHIST_CHK_IERR(SUBR(mvec_random)(vV,ierr),*ierr);
      PHIST_CHK_IERR(SUBR(mvec_normalize)(vV,&nrm,ierr),*ierr);
      rank=0;// dimension of null space
    }
    PHIST_CHK_IERR(SUBR(sdMat_put_value)(R,(ST)nrm,ierr),*ierr);
    *ierr=1-rank;
    return;
  }// case ncols=1: normalize single vector

  PHIST_DEB("mvec_QR: multi-vector case\n");

  if (
  (V->traits.flags&GHOST_DENSEMAT_SCATTERED) ||
  (R->traits.flags&GHOST_DENSEMAT_SCATTERED))
  {
    PHIST_SOUT(PHIST_ERROR,"mvec_QR: cannot handle scattered vectors\n");
    *ierr=-1; // can't handle non-constant stride
    return;
  }//case vectors scattered: not implemented!

    //TSQR for row major storage not available yet - explicit memtranspose
    //     if either mvecs or sdMats are row-major. 
  bool transV=(V->traits.storage==GHOST_DENSEMAT_ROWMAJOR);
  bool transR=(R->traits.storage==GHOST_DENSEMAT_ROWMAJOR);
  
  if (transR||transV)
  {
    PHIST_DEB("we need to make the memory layout of V and/or R conform with TSQR\n");
    if (transV)
    {
      if (V->traits.flags&GHOST_DENSEMAT_VIEW) 
      {
        PHIST_SOUT(PHIST_ERROR,"mvec_QR: cannot handle scattered vectors\n");
        *ierr=-1; // can't handle view + memtranspose in ghost, yet
      }
      PHIST_DEB("memtranspose V\n");
      PHIST_CHK_GERR(V->memtranspose(V),*ierr);
    }
    if (transR)
    {
      PHIST_DEB("memtranspose R\n");
      PHIST_CHK_GERR(R->memtranspose(R),*ierr);
    }
  
    // do not change ierr after this call because
    // it may carry rank information
    int ierr_final;
    SUBR(mvec_QR)(V,R,&ierr_final);
    if (transV)
    {
      PHIST_DEB("memtranspose back V\n");
      PHIST_CHK_GERR(V->memtranspose(V),*ierr);
    }
    if (transR)
    {
      PHIST_DEB("memtranspose back R\n");
      PHIST_CHK_GERR(R->memtranspose(R),*ierr);
    }
    *ierr=ierr_final;
    return;
  }// need memtranspose of V or R
  
  // Here the actual TSQR call with col-major V and R begins...

  // wrapper class for ghost_densemat_t for calling Belos.
  // The wrapper does not own the vector so it doesn't destroy it.
  phist::GhostMV mv_V(V,false);
    
#ifdef TESTING
  int nrows = R->traits.nrows;
  ncols = R->traits.ncols;
    
  PHIST_CHK_IERR(*ierr=nrows-ncols,*ierr);
  PHIST_CHK_IERR(*ierr=nrows-(V->traits.ncols),*ierr);
#endif

#if defined(PHIST_HAVE_TEUCHOS)&&defined(PHIST_HAVE_KOKKOS)

  PHIST_DEB("do TSQR on col-major ghost data structures\n");
  PHIST_DEB("create Teuchos view of R\n");
  Teuchos::RCP<Traits<_ST_ >::Teuchos_sdMat_t> R_view;
  PHIST_CHK_IERR(R_view = Traits<_ST_ >::CreateTeuchosViewNonConst
        (Teuchos::rcp(R,false),ierr),*ierr);

  PHIST_DEB("create TSQR ortho manager\n");  
  Belos::TsqrOrthoManager<_ST_, phist::GhostMV> tsqr("phist/ghost");
  Teuchos::RCP<const Teuchos::ParameterList> valid_params = 
        tsqr.getValidParameters();
  // faster but numerically less robust settings:
  Teuchos::RCP<const Teuchos::ParameterList> fast_params = 
        tsqr.getFastParameters();
  Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp
        (new Teuchos::ParameterList(*valid_params));
  params->set("randomizeNullSpace",true);
  params->set("relativeRankTolerance",rankTol);
  PHIST_DEB("set TSQR parameters\n");
  tsqr.setParameterList(params);

  TRY_CATCH(rank = tsqr.normalize(mv_V,R_view),*ierr);
  PHIST_DEB("V has %d columns and rank %d\n",ncols,rank);
  *ierr = ncols-rank;// return positive number if rank not full.
#else
  *ierr=-99; // no Trilinos, no TSQR, no mvec_QR (right now)
#endif
  return;
  }


//!@}

//! mixed real/complex operation: split mvec into real and imag part.
//! if either reV or imV are NULL, it is not touched.
#ifdef IS_COMPLEX
# ifdef IS_DOUBLE
extern "C" void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, Dmvec_t* reV, Dmvec_t* imV, int *ierr)
{
  *ierr=-99;
}
# else
extern "C" void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, Smvec_t* reV, Smvec_t* imV, int *ierr)
{
  *ierr=-99;
}
# endif
#endif

void SUBR(crsMat_create_fromRowFunc)(TYPE(crsMat_ptr) *vA, const_comm_ptr_t vcomm,
        gidx_t nrows, gidx_t ncols, lidx_t maxnne, 
                int (*rowFunPtr)(ghost_gidx_t,ghost_lidx_t*,ghost_gidx_t*,void*), int *ierr)
{
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
PHIST_GHOST_TASK_BEGIN

  ghost_sparsemat_t* mat;
  ghost_context_t *ctx;

  ghost_sparsemat_traits_t *mtraits=new ghost_sparsemat_traits_t;
        *mtraits=(ghost_sparsemat_traits_t)GHOST_SPARSEMAT_TRAITS_INITIALIZER;
#ifdef PHIST_USE_SELL
        mtraits->format = GHOST_SPARSEMAT_SELL;
        ghost_sell_aux_t aux = GHOST_SELL_AUX_INITIALIZER;
        aux.C = PHIST_SELL_C;
        mtraits->sortScope = PHIST_SELL_SIGMA;
        mtraits->aux = &aux;
        mtraits->flags = (ghost_sparsemat_flags_t)(GHOST_SPARSEMAT_DEFAULT|GHOST_SPARSEMAT_PERMUTE);
#else
#warning "GHOST interface compiled to use the reference CRS format, will probably not yield optimal performance"
        mtraits->format = GHOST_SPARSEMAT_CRS;
        mtraits->flags = (ghost_sparsemat_flags_t)(GHOST_SPARSEMAT_DEFAULT);
#endif
        mtraits->datatype = st::ghost_dt;
  PHIST_CHK_GERR(ghost_context_create(&ctx,0,0,
        GHOST_CONTEXT_DEFAULT,NULL,GHOST_SPARSEMAT_SRC_FUNC,MPI_COMM_WORLD,1.0),*ierr);
  PHIST_CHK_GERR(ghost_sparsemat_create(&mat,ctx,mtraits,1),*ierr);                               

  ghost_sparsemat_src_rowfunc_t src = GHOST_SPARSEMAT_SRC_ROWFUNC_INITIALIZER;
  src.func = rowFunPtr;
  src.maxrowlen = maxnne;
      
  PHIST_CHK_GERR(mat->fromRowFunc(mat,&src),*ierr);
#if PHIST_OUTLEV >= PHIST_VERBOSE
  char *str;
  ghost_context_string(&str,ctx);
  PHIST_SOUT(PHIST_VERBOSE,"%s\n",str);
  free(str); str = NULL;
  ghost_sparsemat_string(&str,mat);
  PHIST_SOUT(PHIST_VERBOSE,"%s\n",str);
  free(str); str = NULL;
#endif
  *vA = (TYPE(crsMat_ptr))mat;
PHIST_GHOST_TASK_END

  return;
}

