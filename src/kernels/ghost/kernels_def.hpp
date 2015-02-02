#include "phist_config.h"
#include "ghost/config.h"
#ifdef PHIST_USE_SELL
#include "ghost/sell.h"
#endif

#ifdef GHOST_HAVE_SCOTCH_disabled
#define USE_SCOTCH
#endif

/* previously we started every kernel function as a task, but since this
   means that the threads need to be created every time, we now just start
   the main program as a task instead
*/   
#ifdef PHIST_GHOST_TASK_BEGIN
#undef PHIST_GHOST_TASK_BEGIN
#endif
#ifdef PHIST_GHOST_TASK_END
#undef PHIST_GHOST_TASK_END
#endif

#define PHIST_GHOST_TASK_BEGIN
#define PHIST_GHOST_TASK_END


// we implement only the double precision real type D
extern "C" void SUBR(type_avail)(int* iflag)
{
  *iflag=0;
}

// \name Matrix input from a file
//@{


//! read a matrix from a MatrixMarket (ASCII) file
extern "C" void SUBR(sparseMat_read_mm)(TYPE(sparseMat_ptr)* vA, const_comm_ptr_t vcomm,
const char* filename,int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
#include "phist_std_typedefs.hpp"
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  PHIST_CAST_PTR_FROM_VOID(const MPI_Comm,comm,vcomm,*iflag);
  if (filename==NULL)
  {
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
  ghost_sparsemat_t* mat;
  ghost_context_t *ctx;

  ghost_sparsemat_traits_t *mtraits=new ghost_sparsemat_traits_t;
        *mtraits=(ghost_sparsemat_traits_t)GHOST_SPARSEMAT_TRAITS_INITIALIZER;
  ghost_sparsemat_flags_t flags=GHOST_SPARSEMAT_DEFAULT;
#ifdef PHIST_USE_SELL
        mtraits->format = GHOST_SPARSEMAT_SELL;
        ghost_sell_aux_t aux = GHOST_SELL_AUX_INITIALIZER;
        aux.C = PHIST_SELL_C;
        mtraits->sortScope = PHIST_SELL_SIGMA;
        mtraits->aux = &aux;
        if (mtraits->sortScope > 1) {
            flags=(ghost_sparsemat_flags_t)(flags|GHOST_SPARSEMAT_PERMUTE);
        }
#else
#warning "GHOST interface compiled to use the reference CRS format, will probably not yield optimal performance"
        mtraits->format = GHOST_SPARSEMAT_CRS;
#endif
#ifdef USE_SCOTCH
        flags = (ghost_sparsemat_flags_t)(flags|GHOST_SPARSEMAT_SCOTCHIFY);
#endif
        mtraits->datatype = st::ghost_dt;
        mtraits->flags = flags;
        char* cfname=const_cast<char*>(filename);
// TODO - check ghost return codes everywhere like this
  PHIST_CHK_GERR(ghost_context_create(&ctx,0,0,
        GHOST_CONTEXT_DEFAULT,cfname,GHOST_SPARSEMAT_SRC_MM,*comm,1.0),*iflag);
  PHIST_CHK_GERR(ghost_sparsemat_create(&mat,ctx,mtraits,1),*iflag);                               
  PHIST_CHK_GERR(mat->fromMM(mat,cfname),*iflag);
#if PHIST_OUTLEV >= PHIST_VERBOSE
  char *str;
  ghost_context_string(&str,ctx);
  PHIST_SOUT(PHIST_VERBOSE,"%s\n",str);
  free(str); str = NULL;
  ghost_sparsemat_string(&str,mat);
  PHIST_SOUT(PHIST_VERBOSE,"%s\n",str);
  free(str); str = NULL;
#endif
  *vA = (TYPE(sparseMat_ptr))mat;
PHIST_GHOST_TASK_END
}

//! read a matrix from a Ghost CRS (binary) file.
extern "C" void SUBR(sparseMat_read_bin)(TYPE(sparseMat_ptr)* vA, const_comm_ptr_t vcomm,
const char* filename,int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
#include "phist_std_typedefs.hpp"
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  PHIST_CAST_PTR_FROM_VOID(const MPI_Comm,comm,vcomm,*iflag);
  if (filename==NULL)
  {
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
  ghost_sparsemat_t* mat;
  ghost_context_t *ctx;

  ghost_sparsemat_traits_t *mtraits=new ghost_sparsemat_traits_t;
        *mtraits=(ghost_sparsemat_traits_t)GHOST_SPARSEMAT_TRAITS_INITIALIZER;
  ghost_sparsemat_flags_t flags=GHOST_SPARSEMAT_DEFAULT;
#ifdef PHIST_USE_SELL
        mtraits->format = GHOST_SPARSEMAT_SELL;
        ghost_sell_aux_t aux = GHOST_SELL_AUX_INITIALIZER;
        aux.C = PHIST_SELL_C;
        mtraits->sortScope = PHIST_SELL_SIGMA;
        mtraits->aux = &aux;
        if (mtraits->sortScope > 1) {
            flags=(ghost_sparsemat_flags_t)(flags|GHOST_SPARSEMAT_PERMUTE);
        }
#else
#warning "GHOST interface compiled to use the reference CRS format, will probably not yield optimal performance"
        mtraits->format = GHOST_SPARSEMAT_CRS;
#endif
#ifdef USE_SCOTCH
        flags = (ghost_sparsemat_flags_t)(flags|GHOST_SPARSEMAT_SCOTCHIFY);
#endif
        mtraits->datatype = st::ghost_dt;
        mtraits->flags = flags;
        char* cfname=const_cast<char*>(filename);
// TODO - check ghost return codes everywhere like this
  PHIST_CHK_GERR(ghost_context_create(&ctx,0,0,
        GHOST_CONTEXT_DEFAULT,cfname,GHOST_SPARSEMAT_SRC_FILE,*comm,1.0),*iflag);
  PHIST_CHK_GERR(ghost_sparsemat_create(&mat,ctx,mtraits,1),*iflag);                               
  PHIST_CHK_GERR(mat->fromFile(mat,cfname),*iflag);
#if PHIST_OUTLEV >= PHIST_VERBOSE
  char *str;
  ghost_context_string(&str,ctx);
  PHIST_SOUT(PHIST_VERBOSE,"%s\n",str);
  free(str); str = NULL;
  ghost_sparsemat_string(&str,mat);
  PHIST_SOUT(PHIST_VERBOSE,"%s\n",str);
  free(str); str = NULL;
#endif
  *vA = (TYPE(sparseMat_ptr))mat;
PHIST_GHOST_TASK_END
}

//! read a matrix from a Harwell-Boeing (HB) file
extern "C" void SUBR(sparseMat_read_hb)(TYPE(sparseMat_ptr)* vA, const_comm_ptr_t vcomm,
const char* filename,int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_TOUCH(vA);
  PHIST_TOUCH(vcomm);
  PHIST_TOUCH(filename);
  *iflag = -99; // not implemented in ghost, use converter script to bin crs
}

//!@}

//! \name get information about the data distribution in a matrix (maps)

//!@{
//! get the row distribution of the matrix
extern "C" void SUBR(sparseMat_get_row_map)(TYPE(const_sparseMat_ptr) vA, const_map_ptr_t* vmap, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_sparsemat_t,A,vA,*iflag);
  ghost_map_t* map = new ghost_map_t;
  map->ctx = A->context;
  map->vtraits_template=phist_default_vtraits();
  *vmap = (const_map_ptr_t)map;
}

//! get column distribution of a matrix
//! we currently treat all maps as the same as we don't allow any fancy
//! operations using them anyway and ghost can handle both halo'd (colmap)
//! and standard (rowmap) vectors in the mvm.
extern "C" void SUBR(sparseMat_get_col_map)(TYPE(const_sparseMat_ptr) vA, const_map_ptr_t* vmap, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_sparsemat_t,A,vA,*iflag);
  ghost_map_t* map = new ghost_map_t;
  map->ctx = A->context;
  map->vtraits_template=phist_default_vtraits();
  *vmap = (const_map_ptr_t)map;
}

//! get the map for vectors x in y=A*x
//! we currently treat all maps as the same as we don't allow any fancy
//! operations using them anyway and ghost can handle both halo'd (colmap)
//! and standard (rowmap) vectors in the mvm.
extern "C" void SUBR(sparseMat_get_domain_map)(TYPE(const_sparseMat_ptr) vA, const_map_ptr_t* vmap, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  SUBR(sparseMat_get_col_map)(vA,vmap,iflag);
}

//! get the map for vectors y in y=A*x
//! we currently treat all maps as the same as we don't allow any fancy
//! operations using them anyway and ghost can handle both halo'd (colmap)
//! and standard (rowmap) vectors in the mvm.
extern "C" void SUBR(sparseMat_get_range_map)(TYPE(const_sparseMat_ptr) vA, const_map_ptr_t* vmap, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  SUBR(sparseMat_get_row_map)(vA,vmap,iflag);
}
//@}

//! \name constructors

//@{
//! create a block-vector. The entries are stored contiguously
//! at val in column major ordering.
extern "C" void SUBR(mvec_create)(TYPE(mvec_ptr)* vV, 
        const_map_ptr_t vmap, int nvec, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
#include "phist_std_typedefs.hpp"
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  PHIST_CAST_PTR_FROM_VOID(const ghost_map_t, map,vmap,*iflag);
  ghost_densemat_t* result;
  ghost_densemat_traits_t vtraits = map->vtraits_template;/*ghost_cloneVtraits(map->vtraits_template);*/
        vtraits.ncols=nvec;
        vtraits.ncolsorig=nvec;
        vtraits.ncolspadded=0;
        vtraits.datatype = st::ghost_dt;
        vtraits.flags = (ghost_densemat_flags_t)(vtraits.flags & ~GHOST_DENSEMAT_VIEW);
  PHIST_CHK_GERR(ghost_densemat_create(&result,map->ctx,vtraits),*iflag);
  ST zero = st::zero();
  // this allocates the vector and fills it with zeros
  PHIST_CHK_GERR(result->fromScalar(result,&zero),*iflag);
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
        int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);

  PHIST_CAST_PTR_FROM_VOID(const ghost_map_t, map,vmap,*iflag);
  ghost_densemat_t* result;
  ghost_densemat_traits_t vtraits = map->vtraits_template;/*ghost_cloneVtraits(map->vtraits_template);*/
        vtraits.ncols=nvec;
        vtraits.datatype = st::ghost_dt;

  PHIST_CHK_GERR(ghost_densemat_create(&result,map->ctx,vtraits),*iflag);

  if ((result->traits.nrows!=result->traits.nrowshalo)||(result->traits.nrowshalo!=lda))
  {
    PHIST_OUT(PHIST_ERROR,"viewing plain data as ghost_vec only works "
                          "for node-local or trivially parallel matrices!\n");
    PHIST_OUT(PHIST_ERROR,"nrows=%" PRlidx ", nrowshalo=%" PRlidx ", lda=%" PRlidx "\n",
        vtraits.nrows,vtraits.nrowshalo,lda);
    *iflag=-1;
    return;
  }

  PHIST_CHK_GERR(result->viewPlain(result,(void*)values,0,0,lda),*iflag);
  *vV=(TYPE(mvec_ptr))(result);
  return;
}


//! create a serial dense n x m matrix on all procs, with column major
//! ordering.
extern "C" void SUBR(sdMat_create)(TYPE(sdMat_ptr)* vM, int nrows, int ncols, 
        const_comm_ptr_t vcomm, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  MPI_Comm comm_self=MPI_COMM_SELF;
  MPI_Comm* comm=&comm_self;
  if (vcomm!=NULL) comm=(MPI_Comm*)vcomm;
#include "phist_std_typedefs.hpp"
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  ghost_densemat_t* result;
  ghost_densemat_traits_t dmtraits=GHOST_DENSEMAT_TRAITS_INITIALIZER;
        dmtraits.nrows=(ghost_lidx_t)nrows;
        dmtraits.ncols=(ghost_lidx_t)ncols;
        dmtraits.datatype=st::ghost_dt;
#ifdef PHIST_SDMATS_ROW_MAJOR
        dmtraits.storage=GHOST_DENSEMAT_ROWMAJOR;
#else
        dmtraits.storage=GHOST_DENSEMAT_COLMAJOR;
#endif
  // I think the sdMat should not have a context
  ghost_context_t* ctx=NULL;
  PHIST_CHK_GERR(ghost_context_create(&ctx,(ghost_gidx_t)nrows, 
        (ghost_gidx_t)ncols, GHOST_CONTEXT_REDUNDANT, 
        NULL, GHOST_SPARSEMAT_SRC_NONE, *comm, 1.0),*iflag);
  ghost_densemat_create(&result,ctx,dmtraits);
  ST zero = st::zero();
  result->fromScalar(result,&zero);
  *vM=(TYPE(sdMat_ptr))result;
PHIST_GHOST_TASK_END
}

void SUBR(sdMat_create_view)(TYPE(sdMat_ptr)* M, const_comm_ptr_t comm,
        _ST_* values, lidx_t lda, int nrows, int ncols,
        int* iflag)
{
  *iflag=-99;
}

//@}

//! retrieve local length of the vectors in V
extern "C" void SUBR(mvec_my_length)(TYPE(const_mvec_ptr) vV, lidx_t* len, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag = 0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_densemat_t,V,vV,*iflag);
  PHIST_CHK_IERR(*iflag=check_local_size(V->traits.nrows),*iflag);
  *len = V->traits.nrows;
}

//! retrieve the map of the vectors in V
extern "C" void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) vV, const_map_ptr_t* vmap, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_densemat_t,V,vV,*iflag);
  ghost_map_t* map = new ghost_map_t;
  map->ctx=V->context; 
  map->vtraits_template=V->traits;
  *vmap=(const_map_ptr_t)map;
}

//! retrieve the comm used for MPI communication in V
extern "C" void SUBR(mvec_get_comm)(TYPE(const_mvec_ptr) vV, const_comm_ptr_t* vcomm, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_densemat_t,V,vV,*iflag);

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
extern "C" void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) vV, int* nvec, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag = 0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_densemat_t,V,vV,*iflag);
  PHIST_CHK_IERR(*iflag=check_local_size(V->traits.ncols),*iflag);
  *nvec = (int)(V->traits.ncols);
}

//! get number of rows in local dense matrix
extern "C" void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) vM, int* nrows, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_densemat_t,M,vM,*iflag);
  PHIST_CHK_IERR(*iflag=check_local_size(M->traits.nrows),*iflag);
  *nrows = (int)(M->traits.nrows);
}
  
//! get number of cols in local dense matrix
extern "C" void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) vM, int* ncols, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_densemat_t,M,vM,*iflag);
  *ncols = (int)(M->traits.ncols);
}


extern "C" void SUBR(mvec_extract_view)(TYPE(mvec_ptr) vV, _ST_** val, lidx_t* lda, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"

  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V, vV, *iflag);
  if (V->traits.flags & GHOST_DENSEMAT_SCATTERED)
  {
    PHIST_OUT(PHIST_ERROR,"%s: cannot view data with non-constant stride using "
        "this function (file %s, line %d)\n", __FUNCTION__, __FILE__, __LINE__);
    *iflag=-1;
    return;
  }
  if (V->val==NULL)
  {
    PHIST_OUT(PHIST_ERROR,"%s, pointer is NULL\n",__FUNCTION__);
    *iflag=-2;
    return;
  }
  PHIST_CHK_GERR(ghost_densemat_valptr(V,(void**)val),*iflag);
  PHIST_CHK_IERR(*iflag=check_local_size(V->traits.nrowspadded),*iflag);

#ifdef PHIST_MVECS_ROW_MAJOR
  *lda = V->traits.ncolspadded;
#else
  *lda = V->traits.nrowspadded;
#endif
}

extern "C" void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) vM, _ST_** val, lidx_t* lda, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,M, vM, *iflag);

  if (M->traits.flags & GHOST_DENSEMAT_SCATTERED)
  {
    PHIST_OUT(PHIST_ERROR,"%s: cannot view data with non-constant stride using "
        "this function (file %s, line %d)\n", __FUNCTION__, __FILE__, __LINE__);
    *iflag=-1;
    return;
  }

  PHIST_CHK_GERR(ghost_densemat_valptr(M,(void**)val),*iflag);

  PHIST_CHK_IERR(*iflag=check_local_size(M->traits.nrowspadded),*iflag);

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

extern "C" void SUBR(mvec_to_device)(TYPE(mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V, vV, *iflag);
  PHIST_SOUT(PHIST_DEBUG,"ghost densemat upload\n"
                         "nrows=%" PRlidx ", ncols=%" PRlidx "\n"
                         "nrowshalo=%" PRlidx "\n"
                         "nrowspadded=%" PRlidx ", ncolspadded=%" PRlidx "\n",
                         V->traits.nrows, V->traits.ncols, 
                         V->traits.nrowshalo,
                         V->traits.nrowspadded, V->traits.ncolspadded);
  PHIST_SOUT(PHIST_DEBUG,"V flags: %d\n",(int)V->traits.flags);
  PHIST_CHK_GERR(V->upload(V),*iflag);
}

extern "C" void SUBR(mvec_from_device)(TYPE(mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V, vV, *iflag);
  PHIST_CHK_GERR(V->download(V),*iflag);
}

extern "C" void SUBR(sdMat_to_device)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,M, vM, *iflag);
  PHIST_CHK_GERR(M->upload(M),*iflag);
}

extern "C" void SUBR(sdMat_from_device)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,M, vM, *iflag);
  PHIST_CHK_GERR(M->download(M),*iflag);
}

extern "C" void SUBR(mvec_to_mvec)(TYPE(const_mvec_ptr) v_in, TYPE(mvec_ptr) v_out, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V_in,v_in,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V_out,v_out,*iflag);
  
  if (V_in->context!=V_out->context)
  {
    PHIST_SOUT(PHIST_WARNING,"function %s only implemented for simple permutation operations\n"
                             "where the result and input vectors have the same context and  \n"
                             "one of them may be permuted\n",__FUNCTION__);
    *iflag=PHIST_NOT_IMPLEMENTED;
    return;
  }
  
  bool resultPermuted=false;//V_out->flags && IS_PERMUTED;
  bool inputPermuted=false;//V_out->flags && IS_PERMUTED;
  
  // first copy the data
  PHIST_CHK_GERR(V_out->fromVec(V_out,V_in,0,0),*iflag);
  // check permutation state
  if (resultPermuted==inputPermuted) return;
  if (resultPermuted)
  {
    PHIST_CHK_GERR(V_out->permute(V_out,V_out->context->permutation,GHOST_PERMUTATION_ORIG2PERM),*iflag);
  }
  else
  {
    PHIST_CHK_GERR(V_out->permute(V_out,V_out->context->permutation,GHOST_PERMUTATION_PERM2ORIG),*iflag);
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
                             int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  ghost_densemat_t *Vblock;
  V->viewCols(V, &Vblock, jmax-jmin+1, jmin);

  if (*vVblock!=NULL)
  {
    PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,tmp,*vVblock,*iflag);
    //PHIST_DEB("destroying previous vector (view)\n");
    tmp->destroy(tmp);
  }
  PHIST_CHK_IERR(*iflag=((Vblock->traits.flags&GHOST_DENSEMAT_VIEW)-GHOST_DENSEMAT_VIEW),*iflag);
  *vVblock = (TYPE(mvec_ptr))Vblock;
}

//! get a new vector that is a copy of some columns of the original one,  
//! Vblock = V(:,jmin:jmax). The object Vblock must be created beforehand 
//! and the corresponding columns of V are copied into the value array    
//! of Vblock. V is not modified.
extern "C" void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) vV,
                             TYPE(mvec_ptr) vVblock,
                             int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,Vblock,vVblock,*iflag);
  *iflag=0;
#ifdef TESTING
// nonzero error code if #vectors in Vblock too small or large
  PHIST_CHK_IERR(*iflag=(jmax-jmin+1)-Vblock->traits.ncols,*iflag);
  PHIST_CHK_IERR(*iflag=((jmin<0)||(jmax>=V->traits.ncols))?PHIST_INVALID_INPUT:0,*iflag);
  // The phist kernel interface requires Vblock to be created by mvec_create, 
  // which calls fromScalar to allocate the block of memory. So we perform a few
  // safety checks here
  PHIST_CHK_IERR(*iflag=(Vblock->traits.nrows==V->traits.nrows)?0:PHIST_INVALID_INPUT,*iflag)
  // not sure what the ghost function fromVec actually supports, but I think 
  // this makes sense:
  PHIST_CHK_IERR(*iflag=(Vblock->traits.nrowspadded==V->traits.nrowspadded)?0:PHIST_INVALID_INPUT,*iflag)
#else
  PHIST_TOUCH(jmax);
#endif  
  Vblock->fromVec(Vblock,V,(ghost_lidx_t)0,(ghost_lidx_t)jmin);
PHIST_GHOST_TASK_END
}

//! given a multi-vector Vblock, set V(:,jmin:jmax)=Vblock by copying the corresponding
//! vectors. Vblock is not modified.
extern "C" void SUBR(mvec_set_block)(TYPE(mvec_ptr) vV,
                             TYPE(const_mvec_ptr) vVblock,
                             int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,Vblock,vVblock,*iflag);

#ifdef TESTING
int nv_v,nv_vb;
lidx_t nr_v,nr_vb;
PHIST_CHK_IERR(*iflag=V->elSize-Vblock->elSize,*iflag);
PHIST_CHK_IERR(*iflag=V->elSize-sizeof(_ST_),*iflag);
PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nv_v,iflag),*iflag);
PHIST_CHK_IERR(SUBR(mvec_num_vectors)(Vblock,&nv_vb,iflag),*iflag);
PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&nr_v,iflag),*iflag);
PHIST_CHK_IERR(SUBR(mvec_my_length)(Vblock,&nr_vb,iflag),*iflag);
  if ((nr_v!=nr_vb) || ((jmax-jmin+1)!=nv_vb) ||
      (jmin<0) || (jmax>nv_v))
      {
        PHIST_SOUT(PHIST_ERROR,"mvec_set_block: you are trying to set\n"
                               "V(%d:%" PRlidx ",%d:%d)=Vb(%d:%" PRlidx ",%d:%d)\n"
                               "(with V %" PRlidx "x%d)\n",
                               1,nr_v,jmin,jmax,1,nr_vb,1,nv_vb,nr_v,nv_v);
        *iflag=PHIST_INVALID_INPUT;
        return;
      }
#endif

  // as ghost uses memcpy, avoid passing in vectors that actually
  // view the same data (this may happen in subspacejada as it is
  // implemented right now)
  _ST_ *v_ptr=NULL,*vb_ptr=NULL;
  lidx_t ldv,ldvb;
  PHIST_CHK_IERR(SUBR(mvec_extract_view)((TYPE(mvec_ptr))vV,&v_ptr,&ldv,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_extract_view)((TYPE(mvec_ptr))vVblock,&vb_ptr,&ldvb,iflag),*iflag);
  
  // as far as I know ghost just calls memcpy(a,b), which
  // (I think) is ill-defined if the arrays overlap. If the
  // data is already in the correct location, return here.
  // otherwise, give a warning.

// row and column stride (rs and cs)
#ifdef PHIST_MVECS_ROW_MAJOR
  const int cs=1,rs=ldv;
#else
  const int cs=ldv,rs=1;
#endif  
  if (vb_ptr == v_ptr + cs*jmin)
  {
    PHIST_SOUT(PHIST_DEBUG,"mvec_set_block: data already in the correct location.\n");
    *iflag=0;
    return;
  }
  else if ( (vb_ptr                   <= v_ptr + cs*jmax) &&
            (vb_ptr+cs*(jmax-jmin)  >= v_ptr + cs*jmin) )
  {
    // note: this test is certainly not 
    // sufficient, but it catches the   
    // situation where Vblock is a view 
    // of V overlapping the target cols.
    PHIST_SOUT(PHIST_ERROR,"mvec_set_block: overlapping memory locations!\n");
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
  
  // create a view of the requested columns of V
  ghost_densemat_t *Vcols=NULL;
  V->viewCols(V,&Vcols,(ghost_lidx_t)(jmax-jmin+1),(ghost_lidx_t)jmin);
#ifdef TESTING
  _ST_* vc_ptr=NULL;
  lidx_t ldvc;
  PHIST_CHK_IERR(SUBR(mvec_extract_view)((TYPE(mvec_ptr))Vcols,&vc_ptr,&ldvc,iflag),*iflag);
  PHIST_CHK_IERR(*iflag=ldvc-ldv,*iflag);
  if (vc_ptr!=v_ptr+jmin*cs)
  {
    PHIST_SOUT(PHIST_ERROR,"tmp view incorrect in mvec_set_block (file %s, line %d)\n",
        __FILE__,__LINE__);
  }
#endif
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
                             int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,M,vM,*iflag);

#ifdef TESTING
  if (imin<0||jmin<0||imax>M->traits.nrows||jmax>M->traits.ncols)
  {
    PHIST_OUT(PHIST_ERROR,"%s: range out of bounds of matrix\n"
                          "requested range: (%d:%d,%d:%d), matrix dim: %dx%d\n",
                          __FUNCTION__,
                          imin,imax,jmin,jmax,M->traits.nrows,M->traits.ncols);
  }
#endif
  //TODO: we only view the host side of the vector here, this function should
  //      eventually be moved into ghost and the accelerator stuff added.

  // first just create a view of the corresponding columns
  ghost_densemat_t *Mblock;
  M->viewVec(M, &Mblock, imax-imin+1,imin,jmax-jmin+1, jmin);

  if (*vMblock!=NULL)
  {
    //PHIST_DEB("deleting previous object in %s\n",__FUNCTION__);
    PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,tmp,*vMblock,*iflag);
    tmp->destroy(tmp);
  }
  PHIST_CHK_IERR(*iflag=((Mblock->traits.flags&GHOST_DENSEMAT_VIEW)-GHOST_DENSEMAT_VIEW),*iflag);
  *vMblock = (TYPE(sdMat_ptr))Mblock;
}

//! get a new matrix that is a copy of some rows and columns of the original one,  
//! Mblock = M(imin:imax,jmin:jmax). The object Mblock must be created beforehand 
//! and the corresponding columns of M are copied into the value array    
//! of Mblock. M is not modified.
extern "C" void SUBR(sdMat_get_block)(TYPE(const_sdMat_ptr) vM,
                             TYPE(sdMat_ptr) vMblock,
                             int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,M,vM,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,Mblock,vMblock,*iflag);

#ifdef TESTING
  int nr=imax-imin+1;
  int nc=jmax-jmin+1;
  if (Mblock->traits.nrows!=nr || Mblock->traits.ncols!=nc)
  {
    PHIST_SOUT(PHIST_ERROR,"result block has wrong dimensions %dx%d, "
                           "requested range is (%d:%d,%d:%d)\n",
                           Mblock->traits.nrows,Mblock->traits.ncols,imin,imax,jmin,jmax);
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
  if (imin<0 || imax>M->traits.nrows ||
      jmin<0 || jmax>M->traits.ncols )
  {
    PHIST_SOUT(PHIST_ERROR,"requested range invalid. M is %dx%d, "
                           "requested range is (%d:%d,%d:%d)\n",
                           M->traits.nrows,M->traits.ncols,imin,imax,jmin,jmax);
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
#endif 
  _ST_ *m_ptr, *mb_ptr;
  lidx_t ldm, ldmb;
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)((void*)M,&m_ptr,&ldm,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)((void*)Mblock,&mb_ptr,&ldmb,iflag),*iflag);
#ifdef PHIST_SDMATS_ROW_MAJOR
// if we ever want that...
#error "row-major sdMats not implemented here"
#endif  
  if (mb_ptr==m_ptr+ldm*jmin+imin)
  {
    PHIST_SOUT(PHIST_DEBUG,"%s: data already in place.\n",__FUNCTION__);
    *iflag=0;
    return;
  }
  //TODO: check for overlapping data regions?
  Mblock->fromVec(Mblock,M,imin,jmin);
PHIST_GHOST_TASK_END
}

//! given a serial dense matrix Mblock, set M(imin:imax,jmin:jmax)=Mblock by 
//! copying the corresponding elements. Mblock is not modified.
extern "C" void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) vM, 
                             TYPE(const_sdMat_ptr) vMblock,
                             int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,Mblock,vMblock,*iflag);

  ghost_densemat_t* Mb_view=NULL;
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(vM,(TYPE(sdMat_ptr)*)&Mb_view,imin,imax,jmin,jmax,iflag),*iflag);
  Mb_view->fromVec(Mb_view,Mblock,0,0);
  Mb_view->destroy(Mb_view);
PHIST_GHOST_TASK_END
}

//! \name destructors

//@{

//!
extern "C" void SUBR(sparseMat_delete)(TYPE(sparseMat_ptr) vA, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  if (vA==NULL) return;
  PHIST_CAST_PTR_FROM_VOID(ghost_sparsemat_t,A,vA,*iflag);
  A->destroy(A);
}

//!
extern "C" void SUBR(mvec_delete)(TYPE(mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  if (vV==NULL) return;
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  V->destroy(V);
}

//!
extern "C" void SUBR(sdMat_delete)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  if (vM==NULL) return;
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,M,vM,*iflag);
  M->destroy(M);
}

//@}

//! \name Numerical functions
//!@{

//! put scalar value into all elements of a multi-vector
extern "C" void SUBR(mvec_put_value)(TYPE(mvec_ptr) vV, _ST_ value, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  PHIST_DEB("put value, V @ %p. V->traits.nrows=%" PRlidx "\n",V,V->traits.nrows);
  V->fromScalar(V,(void*)&value);
PHIST_GHOST_TASK_END
}

//! put scalar value into all elements of a multi-vector
extern "C" void SUBR(sdMat_put_value)(TYPE(sdMat_ptr) vV, _ST_ value, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  V->fromScalar(V,(void*)&value);
PHIST_GHOST_TASK_END
}

//! put random numbers into all elements of a multi-vector
extern "C" void SUBR(mvec_random)(TYPE(mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  V->fromRand(V);
PHIST_GHOST_TASK_END
}

extern "C" void SUBR(mvec_print)(TYPE(const_mvec_ptr) vV, int* iflag)
{
  *iflag = 0;
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  std::cout << "# local rows: "<<V->traits.nrows<<std::endl;
  std::cout << "# vectors:    "<<V->traits.ncols<<std::endl;
  char *str=NULL;
  V->string(V,&str);
  std::cout << str <<std::endl;
  free(str); str = NULL;
}

extern "C" void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) vM, int* iflag)
{
  *iflag=0;
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,M,vM,*iflag);
  std::cout << "# rows: "<<M->traits.nrows<<std::endl;
  std::cout << "# cols: "<<M->traits.ncols<<std::endl;
  char *str=NULL;
  M->string(M,&str);
  std::cout << str <<std::endl;
  free(str); str = NULL;
}


//! put random numbers into all elements of a serial dense matrix
extern "C" void SUBR(sdMat_random)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,M,vM,*iflag);
  M->fromRand(M);
  // use same values on all mpi processes if we have a communicator
  if( M->context && M->context->mpicomm )
    M->syncValues(M, 0);
PHIST_GHOST_TASK_END
}

//! \name Numerical functions

//! compute the 2-norm) of each column of v                   
//! (vnrm[i] must be pre-allocated by caller)
  extern "C" void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) vV,
                            _MT_* vnrm, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
#include "phist_std_typedefs.hpp"  
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  int i,nv;
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);  
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
                            _MT_* vnrm, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"  
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);  
  // TODO - this call doesn't return the norm as we wish
  //V->normalize(V);    
  PHIST_CHK_IERR(SUBR(mvec_norm2)(vV,vnrm,iflag),*iflag);
  _ST_ inrm[V->traits.ncols];
  for (int i=0;i<V->traits.ncols;i++) inrm[i]=st::one()/vnrm[i];
  PHIST_CHK_IERR(SUBR(mvec_vscale)(vV,inrm,iflag),*iflag);
  return;
}

//! scale each column i of v and by scalar[i]
extern "C" void SUBR(mvec_scale)(TYPE(mvec_ptr) vV, 
                            _ST_ scalar, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);  
  V->scale(V,(void*)&scalar);
  return;
PHIST_GHOST_TASK_END
}

//! scale each column i of v and by scalar[i]
extern "C" void SUBR(mvec_vscale)(TYPE(mvec_ptr) vV, 
                            const _ST_* scalar, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
#include "phist_std_typedefs.hpp"  
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);  
  V->vscale(V,(void*)scalar);
  return;
PHIST_GHOST_TASK_END
}

//! y=alpha*x+beta*y
extern "C" void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) vX,
                            _ST_ beta,  TYPE(mvec_ptr)       vY, 
                            int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
#include "phist_std_typedefs.hpp"
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,X,vX,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,Y,vY,*iflag);
  ST a=alpha, b=beta;
  if (alpha==st::zero())
  {
    if (beta!=st::one())
    {
      PHIST_DEB("scale output Y=beta*Y");
      PHIST_CHK_GERR(Y->scale(Y,(void*)&b),*iflag);
    }
  }
  else if (beta==st::one())
  {
    PHIST_DEB("axpy operation: Y=alpha*X+Y");
    PHIST_CHK_GERR(Y->axpy(Y,X,(void*)&a),*iflag);
  }
  else
  {
    PHIST_DEB("general case: Y=alpha*X+beta*Y");
    PHIST_CHK_GERR(Y->axpby(Y,X,(void*)&a,(void*)&b),*iflag);
  }
PHIST_GHOST_TASK_END
}

//! y[j]=alpha[j]*x[j]+beta[j]*y[j] for all columns j
extern "C" void SUBR(mvec_vadd_mvec)(const _ST_ *alpha, TYPE(const_mvec_ptr) vX,
                            _ST_ beta,  TYPE(mvec_ptr)       vY, 
                            int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
#include "phist_std_typedefs.hpp"
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,X,vX,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,Y,vY,*iflag);
  if(beta == st::one())
  {
    Y->vaxpy(Y,X,(void*)alpha);
  }
  else
  {
    // ghost also expects a vector for beta, so construct one:
    int nvec = 0;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vY,&nvec,iflag),*iflag);
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
                            int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(alpha,vA,beta,vB,iflag), *iflag);
}

//! y=alpha*A*x+beta*y.
extern "C" void SUBR(sparseMat_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, 
_ST_ beta, TYPE(mvec_ptr) vy, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=0;

#ifdef PHIST_TIMEMONITOR
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vx, &nvec, iflag), *iflag);
  for(int i = 0; i < nvec; i++)
    phist_totalMatVecCount();
#endif

  PHIST_CAST_PTR_FROM_VOID(ghost_sparsemat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,x,vx,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,y,vy,*iflag);
  if (alpha==st::zero())
  {
    // no MVM needed
    if (beta==st::zero())
    {
      PHIST_CHK_IERR(SUBR(mvec_put_value)(vy,beta,iflag),*iflag);
    }
    else if (beta!=st::one())
    {
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
      y->scale(y,(void*)&beta);
PHIST_GHOST_TASK_END
    }
  }
  else
  {
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
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
      *iflag=ghost_spmv(y,A,x,&spMVM_opts,&alpha,&beta);
    }
    else
    {
      *iflag=ghost_spmv(y,A,x,&spMVM_opts,&beta);
    }
PHIST_GHOST_TASK_END
  }
}

//! y=alpha*A*x+beta*y.
extern "C" void SUBR(sparseMatT_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, 
_ST_ beta, TYPE(mvec_ptr) vy, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=-99;
  return;
}
//! y[i]=alpha*(A*x[i]+shifts[i]*x[i]) + beta*y[i]
extern "C" void SUBR(sparseMat_times_mvec_vadd_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) vA,
        const _ST_ shifts[], TYPE(const_mvec_ptr) vx, _ST_ beta, TYPE(mvec_ptr) vy, int* 
        iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=0;

#ifdef PHIST_TIMEMONITOR
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vx, &nvec, iflag), *iflag);
  for(int i = 0; i < nvec; i++)
    phist_totalMatVecCount();
#endif

  PHIST_CAST_PTR_FROM_VOID(ghost_sparsemat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,x,vx,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,y,vy,*iflag);
  if (alpha==st::zero())
  {
    PHIST_CHK_IERR(SUBR(sparseMat_times_mvec)(alpha,vA,vx,beta,vy,iflag),*iflag);
  }
  else
  {
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
    int nvec;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vx, &nvec, iflag), *iflag);

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
        *iflag=ghost_spmv(y,A,x,&spMVM_opts,&alpha,ghost_shifts);
      }
      else
      {
        *iflag=ghost_spmv(y,A,x,&spMVM_opts,ghost_shifts);
      }
    }
    else if (beta!=st::zero())
    {
      spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_AXPBY);
      if (alpha!=st::one())
      {
        spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_SCALE);
        *iflag=ghost_spmv(y,A,x,&spMVM_opts,&alpha,&beta,ghost_shifts);
      }
      else
      {
        *iflag=ghost_spmv(y,A,x,&spMVM_opts,&beta,ghost_shifts);
      }
    }
    else
    {
      if (alpha!=st::one())
      {
        spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_SCALE);
        *iflag=ghost_spmv(y,A,x,&spMVM_opts,&alpha,ghost_shifts);
      }
      else
      {
        *iflag=ghost_spmv(y,A,x,&spMVM_opts,ghost_shifts);
      }
    }
PHIST_GHOST_TASK_END
  }
}

//! dot product of vectors v_i and w_i, i=1..numvecs
extern "C" void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) vV, TYPE(const_mvec_ptr) vW, _ST_* s, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,W,vW,*iflag);
  ghost_dot(s,V,W);
PHIST_GHOST_TASK_END
  }

//! dense tall skinny matrix-matrix product yielding a serial dense matrix
//! C=alpha*V'*W+beta*C. C is replicated on all MPI processes sharing V and W.
extern "C" void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) vV, TYPE(const_mvec_ptr) vW, _ST_ beta, TYPE(sdMat_ptr) vC, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,W,vW,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,C,vC,*iflag);
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

  PHIST_CHK_GERR(ghost_gemm(C,V,trans,W,(char*)"N",(void*)&alpha,(void*)&beta,GHOST_GEMM_ALL_REDUCE,GHOST_GEMM_DEFAULT),*iflag);
PHIST_GHOST_TASK_END
  }


//! n x m multi-vector times m x k dense matrix gives n x k multi-vector,
//! W=alpha*V*C + beta*W
extern "C" void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) vV,
                                       TYPE(const_sdMat_ptr) vC,
                           _ST_ beta,  TYPE(mvec_ptr) vW,
                                       int* iflag)
{
    PHIST_ENTER_FCN(__FUNCTION__);
    *iflag=0;
#include "phist_std_typedefs.hpp"
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
    PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
    PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,C,vC,*iflag);
    PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,W,vW,*iflag);

    lidx_t nrV,nrW;
    int ncV, ncW, nrC, ncC;
    nrV=V->traits.nrows;  ncV=V->traits.ncols;
    nrW=W->traits.nrows;  ncW=V->traits.ncols;
    nrC=C->traits.nrows;  ncC=V->traits.ncols;

#ifdef TESTING
    PHIST_CHK_IERR(*iflag=nrV-nrW,*iflag);
    PHIST_CHK_IERR(*iflag=nrC-ncV,*iflag);
    PHIST_CHK_IERR(*iflag=ncC-ncW,*iflag);
    //PHIST_DEB("V'C with V %" PRlidx "x%d, C %dx%d and result %" PRlidx "x%d\n", nrV,ncV,nrC,ncC,nrW,ncW);
#endif
    // note: C is replicated, so this operation is a purely local one.
    PHIST_CHK_GERR(ghost_gemm(W,V,(char*)"N",C,(char*)"N",(void*)&alpha,(void*)&beta,GHOST_GEMM_NO_REDUCE,GHOST_GEMM_DEFAULT),*iflag);
PHIST_GHOST_TASK_END
  }
#if 0
#include "../kernels_no_inplace_VC.cpp"
#else
//! C <- V*C
extern "C" void SUBR(mvec_times_sdMat_inplace)(TYPE(mvec_ptr) vV,
                                       TYPE(const_sdMat_ptr) vC,
                                       int* iflag)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
#include "phist_std_typedefs.hpp"
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
    PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
    PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,C,vC,*iflag);

#ifdef TESTING
    int ncV, nrC, ncC;
    ncV=V->traits.ncols;
     nrC=C->traits.nrows;  ncC=V->traits.ncols;

    PHIST_CHK_IERR(*iflag=nrC-ncV,*iflag);
    PHIST_CHK_IERR(*iflag=nrC-ncC,*iflag);
    //PHIST_DEB("V'C with V %" PRlidx "x%d, C %dx%d and result %" PRlidx "x%d\n", nrV,ncV,nrC,ncC,nrW,ncW);
#endif
    // note: C is replicated, so this operation is a purely local one.
    ST alpha=st::one();
    ST beta=st::zero();
    // ghost internally picks the in-place variant of possible
    PHIST_CHK_GERR(ghost_gemm(V,V,(char*)"N",C,(char*)"N",(void*)&alpha,(void*)&beta,GHOST_GEMM_NO_REDUCE,GHOST_GEMM_DEFAULT),*iflag);
PHIST_GHOST_TASK_END
  }
#endif
//! n x m serial dense matrix times m x k serial dense matrix gives n x k sdMat,
//! C=alpha*V*W + beta*C (serial XGEMM wrapper)
extern "C" void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vV,
                                         TYPE(const_sdMat_ptr) vW,
                              _ST_ beta, TYPE(sdMat_ptr) vC,
                                         int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,W,vW,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,C,vC,*iflag);
  char trans[]="N";  
  PHIST_CHK_GERR(ghost_gemm(C,V,trans,W,trans,(void*)&alpha,(void*)&beta,GHOST_GEMM_NO_REDUCE,GHOST_GEMM_DEFAULT),*iflag);
PHIST_GHOST_TASK_END
  }

//! n x m conj. transposed serial dense matrix times m x k serial dense matrix gives m x k sdMat,
//! C=alpha*V*W + beta*C (serial XGEMM wrapper)
extern "C" void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vV,
                                         TYPE(const_sdMat_ptr) vW,
                              _ST_ beta, TYPE(sdMat_ptr) vC,
                                         int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,W,vW,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,C,vC,*iflag);
#ifdef IS_COMPLEX
  char trans[]="C";
#else
  char trans[]="T";
#endif  
  PHIST_CHK_GERR(ghost_gemm(C, V, trans,W, (char*)"N", (void*)&alpha, (void*)&beta, GHOST_GEMM_NO_REDUCE,GHOST_GEMM_DEFAULT),*iflag);
PHIST_GHOST_TASK_END
  }


//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.   
//! Q is computed in place of V. If V does not have full rank, iflag>0   
//! indicates the dimension of the null-space of V. The first m-iflag    
//! columns of Q are an orthogonal basis of the column space of V, the  
//! remaining columns form a basis for the null space.  
extern "C" void SUBR(mvec_QR)(TYPE(mvec_ptr) vV, TYPE(sdMat_ptr) vR, int* iflag)
  {
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,R,vR,*iflag);

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
    PHIST_CHK_IERR(SUBR(mvec_normalize)(vV,&nrm,iflag),*iflag);
    PHIST_DEB("single vector QR, R=%8.4e\n",nrm);
    PHIST_CHK_IERR(SUBR(sdMat_put_value)(R,(ST)nrm,iflag),*iflag);
    rank=1;
    if (nrm<rankTol)
    {
      PHIST_DEB("zero vector detected\n");
      // randomize the vector
      PHIST_CHK_IERR(SUBR(mvec_random)(vV,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_normalize)(vV,&nrm,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_put_value)(R,st::zero(),iflag),*iflag);
      rank=0;// dimension of null space
      PHIST_CHK_IERR(SUBR(sdMat_put_value)(R,st::zero(),iflag),*iflag);
    }
    *iflag=1-rank;
    return;
  }// case ncols=1: normalize single vector

  PHIST_DEB("mvec_QR: multi-vector case\n");

#if defined(PHIST_HAVE_TEUCHOS)&&defined(PHIST_HAVE_KOKKOS)
  if (
  (V->traits.flags&GHOST_DENSEMAT_SCATTERED) ||
  (R->traits.flags&GHOST_DENSEMAT_SCATTERED))
  {
    PHIST_SOUT(PHIST_ERROR,"mvec_QR: cannot handle scattered vectors\n");
    *iflag=-1; // can't handle non-constant stride
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
      PHIST_CXX_TIMER("memtranspose for TSQR");
      PHIST_DEB("memtranspose V\n");
      PHIST_CHK_GERR(V->memtranspose(V),*iflag);
    }
    if (transR)
    {
      PHIST_CXX_TIMER("memtranspose for TSQR");
      PHIST_DEB("memtranspose R\n");
      PHIST_CHK_GERR(R->memtranspose(R),*iflag);
    }
  
    // do not change iflag after this call because
    // it may carry rank information
    int iflag_final;
    SUBR(mvec_QR)(V,R,&iflag_final);
    if (transV)
    {
      PHIST_CXX_TIMER("memtranspose for TSQR");
      PHIST_DEB("memtranspose back V\n");
      PHIST_CHK_GERR(V->memtranspose(V),*iflag);
    }
    if (transR)
    {
      PHIST_CXX_TIMER("memtranspose for TSQR");
      PHIST_DEB("memtranspose back R\n");
      PHIST_CHK_GERR(R->memtranspose(R),*iflag);
    }
    *iflag=iflag_final;
    return;
  }// need memtranspose of V or R
  
  // Here the actual TSQR call with col-major V and R begins...

  // wrapper class for ghost_densemat_t for calling Belos.
  // The wrapper does not own the vector so it doesn't destroy it.
  phist::GhostMV mv_V(V,false);
    
#ifdef TESTING
  int nrows = R->traits.nrows;
  ncols = R->traits.ncols;
    
  PHIST_CHK_IERR(*iflag=nrows-ncols,*iflag);
  PHIST_CHK_IERR(*iflag=nrows-(V->traits.ncols),*iflag);
#endif

  PHIST_DEB("do TSQR on col-major ghost data structures\n");
  PHIST_DEB("create Teuchos view of R\n");
  Teuchos::RCP<Traits<_ST_ >::Teuchos_sdMat_t> R_view;
  PHIST_CHK_IERR(R_view = Traits<_ST_ >::CreateTeuchosViewNonConst
        (Teuchos::rcp(R,false),iflag),*iflag);

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
#if PHIST_OUTLEV>=PHIST_VERBOSE
  static bool firstCall=true;
  if (firstCall)
  {
    std::stringstream ss;
    ss << *params;
    PHIST_SOUT(PHIST_VERBOSE,"TSQR parameters:\n%s\n",ss.str().c_str());
    firstCall=false;
  }
#endif

  PHIST_TRY_CATCH(rank = tsqr.normalize(mv_V,R_view),*iflag);
  PHIST_DEB("V has %d columns and rank %d\n",ncols,rank);
  *iflag = ncols-rank;// return positive number if rank not full.
#else
  *iflag=-99; // no Trilinos, no TSQR, no mvec_QR (right now)
#endif
  return;
  }


//!@}

//! mixed real/complex operation: split mvec into real and imag part.
//! if either reV or imV are NULL, it is not touched.
#ifdef IS_COMPLEX
# ifdef IS_DOUBLE
extern "C" void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, Dmvec_t* reV, Dmvec_t* imV, int *iflag)
{
  *iflag=-99;
}
# else
extern "C" void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, Smvec_t* reV, Smvec_t* imV, int *iflag)
{
  *iflag=-99;
}
# endif
#endif

void SUBR(sparseMat_create_fromRowFunc)(TYPE(sparseMat_ptr) *vA, const_comm_ptr_t vcomm,
        gidx_t nrows, gidx_t ncols, lidx_t maxnne, 
                int (*rowFunPtr)(ghost_gidx_t,ghost_lidx_t*,ghost_gidx_t*,void*), int *iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag = 0;
#include "phist_std_typedefs.hpp"
PHIST_GHOST_TASK_BEGIN
PHIST_GHOST_CHK_IN_TASK(__FUNCTION__, *iflag);

  ghost_sparsemat_t* mat = NULL;
  ghost_context_t *ctx = NULL;
  PHIST_CAST_PTR_FROM_VOID(const MPI_Comm, comm, vcomm, *iflag);

  ghost_sparsemat_traits_t *mtraits=new ghost_sparsemat_traits_t;
        *mtraits=(ghost_sparsemat_traits_t)GHOST_SPARSEMAT_TRAITS_INITIALIZER;
  ghost_sparsemat_flags_t flags=GHOST_SPARSEMAT_DEFAULT;
#ifdef PHIST_USE_SELL
        mtraits->format = GHOST_SPARSEMAT_SELL;
        ghost_sell_aux_t aux = GHOST_SELL_AUX_INITIALIZER;
        aux.C = PHIST_SELL_C;
        mtraits->sortScope = PHIST_SELL_SIGMA;
        mtraits->aux = &aux;
        if (mtraits->sortScope > 1) {
            flags=(ghost_sparsemat_flags_t)(flags|GHOST_SPARSEMAT_PERMUTE);
        }
#else
#warning "GHOST interface compiled to use the reference CRS format, will probably not yield optimal performance"
        mtraits->format = GHOST_SPARSEMAT_CRS;
#endif
#ifdef USE_SCOTCH
        flags = (ghost_sparsemat_flags_t)(flags|GHOST_SPARSEMAT_SCOTCHIFY);
#endif
        mtraits->datatype = st::ghost_dt;
        mtraits->flags = flags;
  PHIST_CHK_GERR(ghost_context_create(&ctx,nrows,ncols,
        GHOST_CONTEXT_DEFAULT,NULL,GHOST_SPARSEMAT_SRC_FUNC,*comm,1.0),*iflag);
  PHIST_CHK_GERR(ghost_sparsemat_create(&mat,ctx,mtraits,1),*iflag);                               

  ghost_sparsemat_src_rowfunc_t src = GHOST_SPARSEMAT_SRC_ROWFUNC_INITIALIZER;
  src.func = rowFunPtr;
  src.maxrowlen = maxnne;
      
  PHIST_CHK_GERR(mat->fromRowFunc(mat,&src),*iflag);
#if PHIST_OUTLEV >= PHIST_VERBOSE
  char *str;
  ghost_context_string(&str,ctx);
  PHIST_SOUT(PHIST_VERBOSE,"%s\n",str);
  free(str); str = NULL;
  ghost_sparsemat_string(&str,mat);
  PHIST_SOUT(PHIST_VERBOSE,"%s\n",str);
  free(str); str = NULL;
#endif
  *vA = (TYPE(sparseMat_ptr))mat;
PHIST_GHOST_TASK_END

  return;
}

