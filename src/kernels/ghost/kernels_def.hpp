#include "phist_config.h"
#include "ghost/config.h"
#include "ghost/sell.h"

#ifdef GHOST_HAVE_SCOTCH
#define USE_SCOTCH
#endif

#if defined(PHIST_HAVE_TEUCHOS)&&defined(PHIST_HAVE_KOKKOS)
template<>
Teuchos::RCP<node_t> ghost::TsqrAdaptor< _ST_ >::node_=Teuchos::null;
#endif

using namespace phist::ghost_internal;


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
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  bool repart = *iflag&PHIST_SPARSEMAT_REPARTITION;
  bool d2clr  = *iflag&PHIST_SPARSEMAT_DIST2_COLOR;
  int outlev = *iflag&PHIST_SPARSEMAT_QUIET ? PHIST_DEBUG : PHIST_INFO;

  int sellC, sellSigma;
  phist::ghost_internal::get_C_sigma(&sellC,&sellSigma,*iflag, *((MPI_Comm*)vcomm));
  PHIST_SOUT(PHIST_INFO, "Creating sparseMat with SELL-%d-%d format.\n", sellC, sellSigma);

  *iflag=0;

PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(const MPI_Comm,comm,vcomm,*iflag);
  if (filename==NULL)
  {
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
  ghost_sparsemat_t* mat;
  ghost_context_t *ctx;

  ghost_sparsemat_traits_t mtraits=(ghost_sparsemat_traits_t)GHOST_SPARSEMAT_TRAITS_INITIALIZER;
  ghost_sparsemat_flags_t flags=GHOST_SPARSEMAT_DEFAULT;

        mtraits.C = sellC;
        mtraits.sortScope = sellSigma;
        if (mtraits.sortScope > 1) {
            flags=(ghost_sparsemat_flags_t)(flags|GHOST_SPARSEMAT_PERMUTE);
        }
      if (repart)
      {
#ifdef USE_SCOTCH
          PHIST_SOUT(outlev, "Trying to repartition the matrix with SCOTCH\n");
          flags = (ghost_sparsemat_flags_t)(flags|GHOST_SPARSEMAT_SCOTCHIFY);
#else
          PHIST_SOUT(outlev, "SCOTCH not available, no matrix repartitioning\n");
#endif
      }
        mtraits.datatype = st::ghost_dt;
        mtraits.flags = flags;
        char* cfname=const_cast<char*>(filename);
// TODO - check ghost return codes everywhere like this
  PHIST_CHK_GERR(ghost_context_create(&ctx,0,0,
        GHOST_CONTEXT_DEFAULT,cfname,GHOST_SPARSEMAT_SRC_MM,*comm,repart ? 0. : 1.),*iflag);
  PHIST_CHK_GERR(ghost_sparsemat_create(&mat,ctx,&mtraits,1),*iflag);                               
  PHIST_CHK_GERR(mat->fromMM(mat,cfname),*iflag);
  char *str;
  ghost_context_string(&str,ctx);
  PHIST_SOUT(outlev,"%s\n",str);
  free(str); str = NULL;
  ghost_sparsemat_string(&str,mat);
  PHIST_SOUT(outlev,"%s\n",str);
  free(str); str = NULL;
//#endif
  *vA = (TYPE(sparseMat_ptr))mat;
PHIST_TASK_END(iflag);
}

//! read a matrix from a Ghost CRS (binary) file.
extern "C" void SUBR(sparseMat_read_bin)(TYPE(sparseMat_ptr)* vA, const_comm_ptr_t vcomm,
const char* filename,int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  bool repart = *iflag&PHIST_SPARSEMAT_REPARTITION;
  bool d2clr  = *iflag&PHIST_SPARSEMAT_DIST2_COLOR;
  int outlev = *iflag&PHIST_SPARSEMAT_QUIET ? PHIST_DEBUG : PHIST_INFO;

  int sellC, sellSigma;
  get_C_sigma(&sellC,&sellSigma,*iflag, *((MPI_Comm*)vcomm));
  PHIST_SOUT(outlev, "Creating sparseMat with SELL-%d-%d format.\n", sellC, sellSigma);

  *iflag=0;

PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(const MPI_Comm,comm,vcomm,*iflag);
  if (filename==NULL)
  {
    *iflag=PHIST_INVALID_INPUT;
    return;
  }

  ghost_sparsemat_t* mat;
  ghost_context_t *ctx;

  ghost_sparsemat_traits_t mtraits=(ghost_sparsemat_traits_t)GHOST_SPARSEMAT_TRAITS_INITIALIZER;
  ghost_sparsemat_flags_t flags=GHOST_SPARSEMAT_DEFAULT;

        mtraits.C = sellC;
        mtraits.sortScope = sellSigma;
        if (mtraits.sortScope > 1) {
            flags=(ghost_sparsemat_flags_t)(flags|GHOST_SPARSEMAT_PERMUTE);
        }

        if (repart)
        {
#ifdef USE_SCOTCH
          PHIST_SOUT(outlev, "Trying to repartition the matrix with SCOTCH\n");
          flags = (ghost_sparsemat_flags_t)(flags|GHOST_SPARSEMAT_SCOTCHIFY);
#else
          PHIST_SOUT(outlev, "SCOTCH not available, no matrix repartitioning\n");
#endif
        }
        mtraits.datatype = st::ghost_dt;
        mtraits.flags = flags;
        char* cfname=const_cast<char*>(filename);
// TODO - check ghost return codes everywhere like this
  PHIST_CHK_GERR(ghost_context_create(&ctx,0,0,
        GHOST_CONTEXT_DEFAULT,cfname,GHOST_SPARSEMAT_SRC_FILE,*comm,repart ? 0. : 1.),*iflag);
  PHIST_CHK_GERR(ghost_sparsemat_create(&mat,ctx,&mtraits,1),*iflag);                               
  PHIST_CHK_GERR(mat->fromFile(mat,cfname),*iflag);
//#if PHIST_OUTLEV >= PHIST_VERBOSE
  char *str;
  ghost_context_string(&str,ctx);
  PHIST_SOUT(outlev,"%s\n",str);
  free(str); str = NULL;
  ghost_sparsemat_string(&str,mat);
  PHIST_SOUT(outlev,"%s\n",str);
  free(str); str = NULL;
//#endif
  *vA = (TYPE(sparseMat_ptr))mat;
PHIST_TASK_END(iflag);
}

//! read a matrix from a Harwell-Boeing (HB) file
extern "C" void SUBR(sparseMat_read_hb)(TYPE(sparseMat_ptr)* vA, const_comm_ptr_t vcomm,
const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_TOUCH(vA);
  PHIST_TOUCH(vcomm);
  PHIST_TOUCH(filename);
  *iflag = PHIST_NOT_IMPLEMENTED; // not implemented in ghost, use converter script to bin crs
}

//!@}

//! \name get information about the data distribution in a matrix (maps)

//!@{
//! get the row distribution of the matrix
extern "C" void SUBR(sparseMat_get_row_map)(TYPE(const_sparseMat_ptr) vA, const_map_ptr_t* vmap, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_sparsemat_t,A,vA,*iflag);
  //!@TODO: cache this as it never gets deleted otherwise!
  ghost_map_t* map = mapGarbageCollector.new_map(vA);
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
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_sparsemat_t,A,vA,*iflag);
  //!@TODO: cache this as it never gets deleted otherwise!
  ghost_map_t* map = mapGarbageCollector.new_map(vA);
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
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  SUBR(sparseMat_get_col_map)(vA,vmap,iflag);
}

//! get the map for vectors y in y=A*x
//! we currently treat all maps as the same as we don't allow any fancy
//! operations using them anyway and ghost can handle both halo'd (colmap)
//! and standard (rowmap) vectors in the mvm.
extern "C" void SUBR(sparseMat_get_range_map)(TYPE(const_sparseMat_ptr) vA, const_map_ptr_t* vmap, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
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
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  bool replicate_cuda_mem=*iflag&PHIST_MVEC_REPLICATE_DEVICE_MEM;
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_MVEC_CREATE(vmap,nvec,iflag);
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(const ghost_map_t, map,vmap,*iflag);
  ghost_densemat_t* result;
  ghost_densemat_traits_t vtraits = map->vtraits_template;/*ghost_cloneVtraits(map->vtraits_template);*/
        vtraits.ncols=nvec;
        vtraits.ncolsorig=nvec;
        vtraits.ncolspadded=0;
        vtraits.datatype = st::ghost_dt;
        vtraits.flags = (ghost_densemat_flags_t)(vtraits.flags & ~GHOST_DENSEMAT_VIEW);

  // on CUDA nodes, allocate only device memory for mvecs
  ghost_type_t ghost_type;
  PHIST_CHK_GERR(ghost_type_get(&ghost_type),*iflag);
  if (ghost_type == GHOST_TYPE_CUDA) 
  {
    if (replicate_cuda_mem)
    {
      // allocate both host and device side (useful for testing)
      vtraits.location = (ghost_location_t)(GHOST_LOCATION_HOST|GHOST_LOCATION_DEVICE);
    }
    else
    {
      vtraits.location = GHOST_LOCATION_DEVICE;
    }
  } 
  else 
  {
    vtraits.location = GHOST_LOCATION_HOST;
  }


  PHIST_CHK_GERR(ghost_densemat_create(&result,map->ctx,vtraits),*iflag);
  ST zero = st::zero();
  // this allocates the vector and fills it with zeros
  PHIST_CHK_GERR(result->fromScalar(result,&zero),*iflag);
  PHIST_DEB("mvec nrows: %" PRlidx "\n",result->traits.nrows);
  *vV=(TYPE(mvec_ptr))(result);
PHIST_TASK_END(iflag);
}

//! create a block-vector as view of raw data. The map tells the object
//! how many rows it should 'see' in the data (at most lda, the leading
//! dimension of the 2D array values). CAVEAT: This function only works
//! if nrowshalo==nrowspadded in the map, which is in general only the case for
//! if there is only 1 MPI process or the matrix is trivially parallel.
extern "C" void SUBR(mvec_create_view)(TYPE(mvec_ptr)* vV, const_map_ptr_t vmap, 
        _ST_* values, lidx_t lda, int nvec,
        int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_CAST_PTR_FROM_VOID(const ghost_map_t, map,vmap,*iflag);
  ghost_densemat_t* result;
  ghost_densemat_traits_t vtraits = map->vtraits_template;/*ghost_cloneVtraits(map->vtraits_template);*/
        vtraits.flags|=GHOST_DENSEMAT_VIEW;
        vtraits.ncols=nvec;
        vtraits.datatype = st::ghost_dt;

  PHIST_CHK_GERR(ghost_densemat_create(&result,map->ctx,vtraits),*iflag);

#ifdef PHIST_MVECS_ROW_MAJOR
  if (result->traits.nrowshalo!=result->traits.nrowspadded+1)
  {
    PHIST_OUT(PHIST_ERROR,"viewing plain data as row-major ghost_densemat_t only works \n"
                          "for vectors without communciation buffers (for spMVM)\n");
    *iflag=-1;
    return;
  }
#else
  if ((result->traits.nrowshalopadded>lda+1))
  {
    PHIST_OUT(PHIST_ERROR,"viewing plain data as ghost_densemat_t only works \n"
                          "if the given lda can accomodate the required comm buffer of the vector!\n"
                          "nrows=%" PRlidx ", nrowshalopadded=%" PRlidx ", lda=%" PRlidx "\n",
        result->traits.nrows,result->traits.nrowshalopadded,lda);
    *iflag=-1;
    return;
  }
#endif
  PHIST_CHK_GERR(result->viewPlain(result,(void*)values,lda),*iflag);
  *vV=(TYPE(mvec_ptr))(result);
  return;
}


//! create a serial dense n x m matrix on all procs, with column major
//! ordering.
extern "C" void SUBR(sdMat_create)(TYPE(sdMat_ptr)* vM, int nrows, int ncols, 
        const_comm_ptr_t vcomm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_SMALL;
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN_SMALLDETERMINISTIC(ComputeTask)
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
  // on CUDA nodes, allocate both host and device memory for sdMats
  ghost_type_t ghost_type;
  PHIST_CHK_GERR(ghost_type_get(&ghost_type),*iflag);
  if (ghost_type == GHOST_TYPE_CUDA) 
  {
    dmtraits.location = (ghost_location_t)(GHOST_LOCATION_HOST|GHOST_LOCATION_DEVICE);
  } 
  else 
  {
    dmtraits.location = GHOST_LOCATION_HOST;
  }

  // I think the sdMat should not have a context
  ghost_densemat_create(&result,NULL,dmtraits);
  ST zero = st::zero();
  PHIST_CHK_GERR(result->fromScalar(result,&zero),*iflag);
  *vM=(TYPE(sdMat_ptr))result;
PHIST_TASK_END(iflag);
}

void SUBR(sdMat_create_view)(TYPE(sdMat_ptr)* M, const_comm_ptr_t comm,
        _ST_* values, lidx_t lda, int nrows, int ncols,
        int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

//@}

//! retrieve the map of the vectors in V
extern "C" void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) vV, const_map_ptr_t* vmap, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_densemat_t,V,vV,*iflag);
  //!@TODO: cache this as it never gets deleted otherwise!
  ghost_map_t* map = mapGarbageCollector.new_map(vV);
  map->ctx=V->context; 
  map->vtraits_template=V->traits;
  *vmap=(const_map_ptr_t)map;
}

//! retrieve number of vectors/columns in V
extern "C" void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) vV, int* nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag = 0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_densemat_t,V,vV,*iflag);
  PHIST_CHK_IERR(*iflag=check_local_size(V->traits.ncols),*iflag);
  *nvec = (int)(V->traits.ncols);
}

//! get number of rows in local dense matrix
extern "C" void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) vM, int* nrows, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_densemat_t,M,vM,*iflag);
  PHIST_CHK_IERR(*iflag=check_local_size(M->traits.nrows),*iflag);
  *nrows = (int)(M->traits.nrows);
}
  
//! get number of cols in local dense matrix
extern "C" void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) vM, int* ncols, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_densemat_t,M,vM,*iflag);
  *ncols = (int)(M->traits.ncols);
}


extern "C" void SUBR(mvec_extract_view)(TYPE(mvec_ptr) vV, _ST_** val, lidx_t* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
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
    if (V->traits.location == GHOST_LOCATION_DEVICE)
    {
      if (V->traits.flags & GHOST_DENSEMAT_NOT_RELOCATE) {      
        PHIST_OUT(PHIST_ERROR,"%s, host side of vector not allocated\n",__FUNCTION__);
        *iflag=PHIST_NOT_IMPLEMENTED;
        return;
      }
      V->download(V);
    }
    else
    {
      PHIST_OUT(PHIST_ERROR,"%s, pointer is NULL\n",__FUNCTION__);
      *iflag=-2;
      return;
    }
  }
  *val=(ST*)V->val;
  PHIST_CHK_IERR(*iflag=check_local_size(V->stride),*iflag);

  *lda = V->stride;
}

extern "C" void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) vM, _ST_** val, lidx_t* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,M, vM, *iflag);

  if (M->traits.flags & GHOST_DENSEMAT_SCATTERED)
  {
    PHIST_OUT(PHIST_ERROR,"%s: cannot view data with non-constant stride using "
        "this function (file %s, line %d)\n", __FUNCTION__, __FILE__, __LINE__);
    *iflag=-1;
    return;
  }

  *val=(_ST_*)M->val;

  PHIST_CHK_IERR(*iflag=check_local_size(M->stride),*iflag);

  *lda = M->stride;
}

extern "C" void SUBR(mvec_to_device)(TYPE(mvec_ptr) vV, int* iflag)
{
  *iflag=0;
#ifdef GHOST_HAVE_CUDA
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  //PHIST_PERFCHECK_VERIFY_TO_DEVICE(vV,iflag);
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
#endif
}

extern "C" void SUBR(mvec_from_device)(TYPE(mvec_ptr) vV, int* iflag)
{
  *iflag=0;
#ifdef GHOST_HAVE_CUDA
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  //PHIST_PERFCHECK_VERIFY_FROM_DEVICE(vV,iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V, vV, *iflag);
  ghost_type_t ghost_type;
  PHIST_CHK_GERR(ghost_type_get(&ghost_type),*iflag);
  if (ghost_type == GHOST_TYPE_CUDA) 
  {
    PHIST_CHK_GERR(V->download(V),*iflag);
  }
#endif
}

extern "C" void SUBR(sdMat_to_device)(TYPE(sdMat_ptr) vM, int* iflag)
{
  *iflag=0;
#ifdef GHOST_HAVE_CUDA
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,M, vM, *iflag);
  ghost_type_t ghost_type;
  PHIST_CHK_GERR(ghost_type_get(&ghost_type),*iflag);
  if (ghost_type == GHOST_TYPE_CUDA) 
  {
    PHIST_CHK_GERR(M->upload(M),*iflag);
  }
#endif
}

extern "C" void SUBR(sdMat_from_device)(TYPE(sdMat_ptr) vM, int* iflag)
{
  *iflag=0;
#ifdef GHOST_HAVE_CUDA
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,M, vM, *iflag);
  PHIST_CHK_GERR(M->download(M),*iflag);
#endif
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
//    PHIST_CHK_GERR(V_out->permute(V_out,GHOST_PERMUTATION_ORIG2PERM),*iflag);
  }
  else
  {
//    PHIST_CHK_GERR(V_out->permute(V_out,GHOST_PERMUTATION_PERM2ORIG),*iflag);
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
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  ghost_densemat_t *Vblock=(ghost_densemat_t*)(*vVblock);
  
  if (Vblock!=NULL)
  {
/*
    if ( &(V->val[0][0]) == &(Vblock->val[0][0]) )
    {
      // if the vector is already a view of some columns in the target
      // vector, update the view.
      PHIST_DEB("update existing view\n");
      int coffs_old = (Vblock->val[0] - Vblock->src->val[0])/Vblock->elSize;
      int coffs=jmin-coffs_old;
      Vblock->val[0] = Vblock->src->val[0]+jmin*Vblock->elSize;
      Vblock->traits.ncols = jmax-jmin+1;
    }
    else
    {
      // delete the vector pointed to by vVblock and create a new view
      PHIST_DEB("delete existing view\n");
    }
  */

    mapGarbageCollector.delete_maps(Vblock);
    Vblock->destroy(Vblock);
    Vblock=NULL;
  }
  PHIST_CHK_GERR(V->viewCols(V, &Vblock, jmax-jmin+1, jmin),*iflag);

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
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_MVEC_GET_BLOCK(vV,vVblock,jmin,jmax,iflag);
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
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
  //PHIST_CHK_IERR(*iflag=(Vblock->traits.nrowspadded==V->traits.nrowspadded)?0:PHIST_INVALID_INPUT,*iflag)
#else
  PHIST_TOUCH(jmax);
#endif  
  PHIST_CHK_GERR(Vblock->fromVec(Vblock,V,(ghost_lidx_t)0,(ghost_lidx_t)jmin),*iflag);
PHIST_TASK_END(iflag);
}

//! given a multi-vector Vblock, set V(:,jmin:jmax)=Vblock by copying the corresponding
//! vectors. Vblock is not modified.
extern "C" void SUBR(mvec_set_block)(TYPE(mvec_ptr) vV,
                             TYPE(const_mvec_ptr) vVblock,
                             int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_MVEC_SET_BLOCK(vV,vVblock,jmin,jmax,iflag);
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
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
  PHIST_CHK_GERR(Vcols->fromVec(Vcols,Vblock,0,0),*iflag);
  // delete the view
  Vcols->destroy(Vcols);
PHIST_TASK_END(iflag);
}

//! get a new sdMat that is a view of some rows and columns of the original one,
//! Mblock = M(imin:imax,jmin:jmax). The new object Vblock is created but does not
//! allocate memory for the vector entries, instead using the entries from V
//! directly.
extern "C" void SUBR(sdMat_view_block)(TYPE(mvec_ptr) vM, TYPE(mvec_ptr)* vMblock,
                             int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
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
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_SMALL;
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN_SMALLDETERMINISTIC(ComputeTask)
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
  PHIST_CHK_GERR(Mblock->fromVec(Mblock,M,imin,jmin),*iflag);
PHIST_TASK_END(iflag);
}

//! given a serial dense matrix Mblock, set M(imin:imax,jmin:jmax)=Mblock by 
//! copying the corresponding elements. Mblock is not modified.
extern "C" void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) vM, 
                             TYPE(const_sdMat_ptr) vMblock,
                             int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_SMALL;
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN_SMALLDETERMINISTIC(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,Mblock,vMblock,*iflag);

  ghost_densemat_t* Mb_view=NULL;
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(vM,(TYPE(sdMat_ptr)*)&Mb_view,imin,imax,jmin,jmax,iflag),*iflag);
  PHIST_CHK_GERR(Mb_view->fromVec(Mb_view,Mblock,0,0),*iflag);
  Mb_view->destroy(Mb_view);
PHIST_TASK_END(iflag);
}

//! \name destructors

//@{

//!
extern "C" void SUBR(sparseMat_delete)(TYPE(sparseMat_ptr) vA, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  if (vA==NULL) return;
  PHIST_CAST_PTR_FROM_VOID(ghost_sparsemat_t,A,vA,*iflag);

  mapGarbageCollector.delete_maps(vA);
  ghost_context_t *ctx = A->context;
  A->destroy(A);
  ghost_context_destroy(ctx);
}

//!
extern "C" void SUBR(mvec_delete)(TYPE(mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  if (vV==NULL) return;
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);

  mapGarbageCollector.delete_maps(vV);
  V->destroy(V);
}

//!
extern "C" void SUBR(sdMat_delete)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  if (vM==NULL) return;
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,M,vM,*iflag);
  ghost_context_t *ctx = NULL;
  if( !(M->traits.flags & GHOST_DENSEMAT_VIEW) )
    ctx = M->context;
  M->destroy(M);
  if( ctx != NULL )
    ghost_context_destroy(ctx);
}

//@}

//! \name Numerical functions
//!@{

//! put scalar value into all elements of a multi-vector
extern "C" void SUBR(mvec_put_value)(TYPE(mvec_ptr) vV, _ST_ value, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_MVEC_PUT_VALUE(vV,iflag);
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  PHIST_DEB("put value, V @ %p. V->traits.nrows=%" PRlidx "\n",V,V->traits.nrows);
  V->fromScalar(V,(void*)&value);
PHIST_TASK_END(iflag);
}

extern "C" void SUBR(mvec_put_func)(TYPE(mvec_ptr) vV,
        phist_mvec_elemFunc funPtr, void* last_arg, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_MVEC_PUT_VALUE(vV,iflag);
PHIST_TASK_DECLARE(ComputeTask)
  PHIST_TASK_BEGIN(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  PHIST_CHK_GERR(V->fromFunc(V,funPtr,last_arg),*iflag);
  PHIST_TASK_END(iflag);
}

//! put scalar value into all elements of a multi-vector
extern "C" void SUBR(sdMat_put_value)(TYPE(sdMat_ptr) vV, _ST_ value, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_SMALL;
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN_SMALLDETERMINISTIC(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  V->fromScalar(V,(void*)&value);
PHIST_TASK_END(iflag);
}

//! put scalar value into all elements of a multi-vector
extern "C" void SUBR(sdMat_identity)(TYPE(sdMat_ptr) V, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag = 0;
  PHIST_PERFCHECK_VERIFY_SMALL;

  _ST_ *V_raw = NULL;
  lidx_t lda;
  int m, n;
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(V, &V_raw, &lda, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(V, &m, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(V, &n, iflag), *iflag);
  for(int i = 0; i < m; i++)
    for(int j = 0; j < n; j++)
      V_raw[lda*i+j] = (i==j) ? st::one() : st::zero();
  PHIST_CHK_IERR(SUBR(sdMat_to_device)(V,iflag),*iflag);
}

#ifndef PHIST_BUILTIN_RNG
//! put random numbers into all elements of a multi-vector
extern "C" void SUBR(mvec_random)(TYPE(mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_MVEC_PUT_VALUE(vV,iflag);
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  V->fromRand(V);
PHIST_TASK_END(iflag);
}
#endif

extern "C" void SUBR(mvec_print)(TYPE(const_mvec_ptr) vV, int* iflag)
{
  *iflag = 0;
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  std::cout << "# local rows: "<<V->traits.nrows<<std::endl;
  std::cout << "# vectors:    "<<V->traits.ncols<<std::endl;
  std::cout << "# row major:  "<<(V->traits.storage & GHOST_DENSEMAT_ROWMAJOR)<<std::endl;
  std::cout << "# stride:     "<<V->stride<<std::endl;

  // if this is a GPU process, do not print the vector values.
  // the vec->string function will download the vector elements,
  // which allocates memory on the CPU. Also, the download changes
  // the semantic of the program. If the user wants to print vector
  // elements, we should add an input flag like PHIST_FORCE.
  if (V->traits.location == GHOST_LOCATION_HOST) 
  {
    char *str=NULL;
    V->string(V,&str);
    std::cout << str <<std::endl;
    free(str); str = NULL;
  }
  else
  {
    // cuda process
    std::cout << "(not printing vector entries on GPU nodes)\n";
  }
}

extern "C" void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) vM, int* iflag)
{
  *iflag=0;
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,M,vM,*iflag);
  std::cout << "# rows:       "<<M->traits.nrows<<std::endl;
  std::cout << "# cols:       "<<M->traits.ncols<<std::endl;
  std::cout << "# row major:  "<<(M->traits.storage & GHOST_DENSEMAT_ROWMAJOR)<<std::endl;
  std::cout << "# stride:     "<<M->stride<<std::endl;
  char *str=NULL;
  M->string(M,&str);
  std::cout << str <<std::endl;
  free(str); str = NULL;
}

#ifndef PHIST_BUILTIN_RNG
//! put random numbers into all elements of a serial dense matrix
extern "C" void SUBR(sdMat_random)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_SMALL
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN_SMALLDETERMINISTIC(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,M,vM,*iflag);
  M->fromRand(M);
PHIST_TASK_END(iflag);
}
#endif

//! put random numbers into all elements of a serial dense matrix
extern "C" void SUBR(sdMat_sync_values)(TYPE(sdMat_ptr) vM, const_comm_ptr_t vcomm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN_SMALLDETERMINISTIC(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,M,vM,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_mpi_comm_t,comm,vcomm,*iflag);
  M->syncValues(M, *comm, 0);
PHIST_TASK_END(iflag);
}

//! \name Numerical functions

//! compute the 2-norm) of each column of v                   
//! (vnrm[i] must be pre-allocated by caller)
  extern "C" void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) vV,
                            _MT_* vnrm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
#include "phist_std_typedefs.hpp" 
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);  
  {
    _ST_ tmp[V->traits.ncols];
    // call mvec_dot_mvec, so communication can be done asynchronously
    //ghost_dot(tmp,V,V);
    PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(vV,vV,tmp,iflag),*iflag);
    for(int i = 0; i < V->traits.ncols; i++)
      vnrm[i] = mt::sqrt(st::real(tmp[i]));
  }
}

//! normalize (in the 2-norm) each column of v and return ||v||_2
//! for each vector i in vnrm[i] (must be pre-allocated by caller)
extern "C" void SUBR(mvec_normalize)(TYPE(mvec_ptr) vV,
                            _MT_* vnrm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"  
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);  
  // TODO - this call doesn't return the norm as we wish
  //V->normalize(V);    
  PHIST_CHK_IERR(SUBR(mvec_norm2)(vV,vnrm,iflag),*iflag);
  _ST_ inrm[V->traits.ncols];
  for (int i=0;i<V->traits.ncols;i++) inrm[i]=st::one()/vnrm[i];
  PHIST_CHK_IERR(SUBR(mvec_vscale)(vV,inrm,iflag),*iflag);
}

//! scale each column i of v and by scalar[i]
extern "C" void SUBR(mvec_scale)(TYPE(mvec_ptr) vV, 
                            _ST_ scalar, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_MVEC_SCALE(vV,iflag);
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);  
  V->scale(V,(void*)&scalar);
PHIST_TASK_END(iflag);
}

//! scale each column i of v and by scalar[i]
extern "C" void SUBR(mvec_vscale)(TYPE(mvec_ptr) vV, 
                            const _ST_* scalar, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"  
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_MVEC_SCALE(vV,iflag);
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);  
  V->vscale(V,(void*)scalar);
PHIST_TASK_END(iflag);
}

//! y=alpha*x+beta*y
extern "C" void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) vX,
                            _ST_ beta,  TYPE(mvec_ptr)       vY, 
                            int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_MVEC_ADD_MVEC(alpha,vX,beta,vY,iflag);
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,X,vX,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,Y,vY,*iflag);
  ST a=alpha, b=beta;
  if (alpha==st::one() && beta==st::zero())
  {
    // copy operation
      PHIST_DEB("copy Y=X");
    PHIST_CHK_GERR(Y->fromVec(Y,X,0,0),*iflag);
  }
  else if (alpha==st::zero())
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
PHIST_TASK_END(iflag);
}

//! y[j]=alpha[j]*x[j]+beta[j]*y[j] for all columns j
extern "C" void SUBR(mvec_vadd_mvec)(const _ST_ *alpha, TYPE(const_mvec_ptr) vX,
                            _ST_ beta,  TYPE(mvec_ptr)       vY, 
                            int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_MVEC_VADD_MVEC(alpha,vX,beta,vY,iflag);
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
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
PHIST_TASK_END(iflag);
}

//! B=alpha*A+beta*B
extern "C" void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA,
                            _ST_ beta,  TYPE(sdMat_ptr)       vB,
                            int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(alpha,vA,beta,vB,iflag), *iflag);
}

//! B=alpha*A+beta*B
extern "C" void SUBR(sdMatT_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA,
                            _ST_ beta,  TYPE(sdMat_ptr)       vB,
                            int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  // simple workaround
  TYPE(sdMat_ptr) I = NULL;
  int m = 0;
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(vA,&m,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&I,m,m,NULL,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_identity)(I,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMatT_times_sdMat)(alpha,vA,I,beta,vB,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(I,iflag),*iflag);
}

//! spMVM communication
extern "C" void SUBR(sparseMat_times_mvec_communicate)(TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=0;

  PHIST_CAST_PTR_FROM_VOID(ghost_sparsemat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,x,vx,*iflag);

#ifdef GHOST_HAVE_MPI
    ghost_densemat_halo_comm_t comm = GHOST_DENSEMAT_HALO_COMM_INITIALIZER;
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
    PHIST_CHK_GERR(x->halocommInit(x,&comm),*iflag);
PHIST_TASK_END(iflag)
PHIST_TASK_POST_STEP(iflag)
    PHIST_CHK_GERR(x->halocommStart(x,&comm),*iflag);
    PHIST_CHK_GERR(x->halocommFinalize(x,&comm),*iflag);
#else
PHIST_TASK_POST_STEP(iflag)
#endif
}

//! y=alpha*A*x+beta*y.
extern "C" void SUBR(sparseMat_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, 
_ST_ beta, TYPE(mvec_ptr) vy, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sparseMat_times_mvec_fused_dot_norm2)(alpha,vA,vx,beta,vy,NULL,NULL,iflag),*iflag);
}

//! y=alpha*A*x+beta*y, ||y||
extern "C" void SUBR(sparseMat_times_mvec_fused_norm2)(_ST_ alpha, TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, 
_ST_ beta, TYPE(mvec_ptr) vy, _MT_* ynrm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sparseMat_times_mvec_fused_dot_norm2)(alpha,vA,vx,beta,vy,NULL,ynrm,iflag),*iflag);
}

//! y=alpha*A*x+beta*y, y'x
extern "C" void SUBR(sparseMat_times_mvec_fused_dot)(_ST_ alpha, TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, 
_ST_ beta, TYPE(mvec_ptr) vy, _ST_* yx, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sparseMat_times_mvec_fused_dot_norm2)(alpha,vA,vx,beta,vy,yx,NULL,iflag),*iflag);
}

extern "C" void SUBR(sparseMat_times_mvec_fused_dot_norm2)(_ST_ alpha, TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, 
_ST_ beta, TYPE(mvec_ptr) vy, _ST_* ydotx, _MT_* ynrm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  ghost_spmv_flags_t spMVM_opts=GHOST_SPMV_DEFAULT;
  // ghost spmvm mode
  if( *iflag & PHIST_SPMVM_ONLY_LOCAL )
    spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_MODE_NOMPI);
  else if( *iflag & PHIST_SPMVM_VECTOR )
    spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_MODE_VECTOR);
  else if( *iflag & PHIST_SPMVM_OVERLAP )
    spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_MODE_OVERLAP);
  else if( *iflag & PHIST_SPMVM_TASK )
    spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_MODE_TASK);
  else
  {
    // defaults to vector
    spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_MODE_VECTOR);
  }
  *iflag=0;

  PHIST_COUNT_MATVECS(vx);

  PHIST_CAST_PTR_FROM_VOID(ghost_sparsemat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,x,vx,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,y,vy,*iflag);

  int nvec = 0;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vy,&nvec,iflag),*iflag);

  if (alpha==st::zero())
  {
    // no MVM needed
    if (beta==st::zero())
    {
      PHIST_CHK_IERR(SUBR(mvec_put_value)(vy,beta,iflag),*iflag);
      if( ydotx != NULL )
        for(int i = 0; i < nvec; i++)
          ydotx[i] = st::zero();
      if( ynrm != NULL )
        for(int i = 0; i < nvec; i++)
          ynrm[i] = mt::zero();
    }
    else
    {
      if (beta!=st::one())
      {
        PHIST_CHK_IERR(SUBR(mvec_scale)(vy,beta,iflag),*iflag);
      }
      if( ydotx != NULL )
      {
        PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(vy,vx,ydotx,iflag),*iflag);
      }
      if( ynrm != NULL )
      {
        PHIST_CHK_IERR(SUBR(mvec_norm2)(vy,ynrm,iflag),*iflag);
      }
    }
  }
  else
  {
  PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)

    // gather args
    int narg = 0;
    void* args[4] = {NULL,NULL,NULL,NULL};

    if (alpha!=st::one())
    {
      spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_SCALE);
      args[narg++] = &alpha;
    }

    if (beta==st::one())
    {
      spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_AXPY);
    }
    else if (beta!=st::zero())
    {
      spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_AXPBY);
      args[narg++] = &beta;
    }

    if( ydotx != NULL )
    {
      spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_DOT_XY);
    }
    if( ynrm != NULL )
    {
      spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_DOT_YY);
    }

    std::vector<_ST_> dotBuff;
    if( ydotx != NULL || ynrm != NULL )
    {
      dotBuff.resize(3*nvec);
      args[narg++] = &dotBuff[0];
    }


    // call ghosts spMV
    PHIST_CHK_GERR(ghost_spmv(y,A,x,&spMVM_opts,args[0],args[1],args[2],args[3]),*iflag);

    if( ynrm != NULL )
      for(int i = 0; i < nvec; i++)
        ynrm[i] = mt::sqrt(st::real(dotBuff[i]));

    if( ydotx != NULL )
      for(int i = 0; i < nvec; i++)
        ydotx[i] = st::conj(dotBuff[nvec+i]);
    
PHIST_TASK_END(iflag);
  }
}

//! y=alpha*A*x+beta*y.
extern "C" void SUBR(sparseMatT_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, 
_ST_ beta, TYPE(mvec_ptr) vy, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
}
//! y[i]=alpha*(A*x[i]+shifts[i]*x[i]) + beta*y[i]
extern "C" void SUBR(sparseMat_times_mvec_vadd_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) vA,
        const _ST_ shifts[], TYPE(const_mvec_ptr) vx, _ST_ beta, TYPE(mvec_ptr) vy, int* 
        iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  ghost_spmv_flags_t spMVM_opts=GHOST_SPMV_DEFAULT;
  // ghost spmvm mode
  if( *iflag & PHIST_SPMVM_ONLY_LOCAL )
    spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_MODE_NOMPI);
  else if( *iflag & PHIST_SPMVM_VECTOR )
    spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_MODE_VECTOR);
  else if( *iflag & PHIST_SPMVM_OVERLAP )
    spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_MODE_OVERLAP);
  else if( *iflag & PHIST_SPMVM_TASK )
    spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_MODE_TASK);
  else
  {
    // defaults to vector
    spMVM_opts = (ghost_spmv_flags_t)((int)spMVM_opts | (int)GHOST_SPMV_MODE_VECTOR);
  }
  *iflag=0;

  PHIST_COUNT_MATVECS(vx);

  PHIST_CAST_PTR_FROM_VOID(ghost_sparsemat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,x,vx,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,y,vy,*iflag);
  if (alpha==st::zero())
  {
    PHIST_CHK_IERR(SUBR(sparseMat_times_mvec)(alpha,vA,vx,beta,vy,iflag),*iflag);
  }
  else
  {
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
    int nvec;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vx, &nvec, iflag), *iflag);

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
PHIST_TASK_END(iflag);
  }
}

//! dot product of vectors v_i and w_i, i=1..numvecs
extern "C" void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) vV, TYPE(const_mvec_ptr) vW, _ST_* s, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"  
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_MVEC_DOT_MVEC(vV,vW,iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,W,vW,*iflag);
  // NOTE: calculate local dot by hand and do the reduction by hand
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  //ghost_dot(s,V,W);
  PHIST_CHK_GERR(ghost_localdot(s,V,W),*iflag);
PHIST_TASK_END(iflag);

PHIST_TASK_POST_STEP(iflag);

#ifdef GHOST_HAVE_MPI
  if (V->context) {
    ghost_mpi_op_t sumOp;
    ghost_mpi_datatype_t mpiDt;
    ghost_mpi_op_sum(&sumOp,V->traits.datatype);
    ghost_mpi_datatype(&mpiDt,V->traits.datatype);
    PHIST_CHK_IERR(*iflag = MPI_Allreduce(MPI_IN_PLACE, s, V->traits.ncols, mpiDt, sumOp, V->context->mpicomm), *iflag);
  }
#endif
}

//! dense tall skinny matrix-matrix product yielding a serial dense matrix
//! C=alpha*V'*W+beta*C. C is replicated on all MPI processes sharing V and W.
extern "C" void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) vV, TYPE(const_mvec_ptr) vW, _ST_ beta, TYPE(sdMat_ptr) vC, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"  
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_MVECT_TIMES_MVEC(vV,vW,iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,W,vW,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,C,vC,*iflag);
#ifdef IS_COMPLEX
  char trans[]="C";
#else
  char trans[]="T";
#endif  
  _ST_ mybeta = st::zero();
  int rank = 0;
  PHIST_CHK_IERR(*iflag = MPI_Comm_rank(V->context->mpicomm,&rank),*iflag);
  if( rank == 0 )
    mybeta = beta;

    lidx_t ncC = C->traits.ncols;

  PHIST_DEB("VtV=C, V %" PRlidx "x%" PRlidx ", \n"
            "       W %" PRlidx "x%" PRlidx ", \n"
            "       C %" PRlidx "x%" PRlidx "\n", 
  V->traits.nrows,V->traits.ncols,
  W->traits.nrows,W->traits.ncols,
  C->traits.nrows,C->traits.ncols);

  // NOTE: we call the allreduction by hand afterwards to allow asynchronuous communication!
  /*
  PHIST_CHK_GERR(ghost_gemm(C,V,trans,W,(char*)"N",(void*)&alpha,(void*)&beta,GHOST_GEMM_ALL_REDUCE,GHOST_GEMM_DEFAULT),*iflag);
  */
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  ghost_error_t gemm_err = ghost_gemm(C,V,trans,W,(char*)"N",(void*)&alpha,(void*)&mybeta,GHOST_GEMM_NO_REDUCE,GHOST_GEMM_DEFAULT);
  if( gemm_err == GHOST_ERR_NOT_IMPLEMENTED )
  {
    // copy result
    ghost_densemat_t* Ccopy=NULL;
    ghost_densemat_traits_t vtraits = C->traits;
    vtraits.storage=GHOST_DENSEMAT_ROWMAJOR;  
    vtraits.flags = (ghost_densemat_flags_t)((int)vtraits.flags & ~(int)GHOST_DENSEMAT_VIEW);
    vtraits.ncolsorig=vtraits.ncols;
    vtraits.nrowsorig=vtraits.nrows;
    ghost_densemat_create(&Ccopy,C->context,vtraits);

    // this allocates the memory for the vector, copies and memTransposes the data
    PHIST_CHK_GERR(Ccopy->fromVec(Ccopy,C,0,0),*iflag);

    PHIST_CHK_GERR(gemm_err = ghost_gemm(Ccopy,V,trans,W,(char*)"N",(void*)&alpha,(void*)&mybeta,GHOST_GEMM_NO_REDUCE,GHOST_GEMM_DEFAULT),*iflag);

    // memtranspose data
    PHIST_CHK_GERR(C->fromVec(C,Ccopy,0,0),*iflag);
    Ccopy->destroy(Ccopy);
  }
  PHIST_CHK_GERR(gemm_err,*iflag);
PHIST_TASK_END(iflag);

PHIST_TASK_POST_STEP(iflag);

  PHIST_CHK_GERR(C->reduce(C,V->context->mpicomm,GHOST_ALLREDUCE),*iflag);
}


//! n x m multi-vector times m x k dense matrix gives n x k multi-vector,
//! W=alpha*V*C + beta*W
extern "C" void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) vV,
                                       TYPE(const_sdMat_ptr) vC,
                           _ST_ beta,  TYPE(mvec_ptr) vW,
                                       int* iflag)
{
    PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
    *iflag=0;
    PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT(alpha,vV,beta,vW,iflag);
  PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
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
PHIST_TASK_END(iflag);
  }
#ifndef PHIST_MVECS_ROW_MAJOR
#warning "Disabling GHOST tsmm_inplace for col-major mvecs for now..."
#include "../common/kernels_no_inplace_VC.cpp"
#else
//! C <- V*C
extern "C" void SUBR(mvec_times_sdMat_inplace)(TYPE(mvec_ptr) vV,
                                       TYPE(const_sdMat_ptr) vC,
                                       int* iflag)
  {
    PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
    *iflag=0;
    PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT_INPLACE(vV,vC,iflag);
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
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
PHIST_TASK_END(iflag);
  }
#endif
//! n x m serial dense matrix times m x k serial dense matrix gives n x k sdMat,
//! C=alpha*V*W + beta*C (serial XGEMM wrapper)
extern "C" void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vV,
                                         TYPE(const_sdMat_ptr) vW,
                              _ST_ beta, TYPE(sdMat_ptr) vC,
                                         int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_SMALL;
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN_SMALLDETERMINISTIC(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,W,vW,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,C,vC,*iflag);
  char trans[]="N";  
  PHIST_CHK_GERR(ghost_gemm(C,V,trans,W,trans,(void*)&alpha,(void*)&beta,GHOST_GEMM_NO_REDUCE,GHOST_GEMM_DEFAULT),*iflag);
PHIST_TASK_END(iflag);
  }

//! n x m conj. transposed serial dense matrix times m x k serial dense matrix gives m x k sdMat,
//! C=alpha*V*W + beta*C (serial XGEMM wrapper)
extern "C" void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vV,
                                         TYPE(const_sdMat_ptr) vW,
                              _ST_ beta, TYPE(sdMat_ptr) vC,
                                         int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_SMALL;
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN_SMALLDETERMINISTIC(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,W,vW,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,C,vC,*iflag);
#ifdef IS_COMPLEX
  char trans[]="C";
#else
  char trans[]="T";
#endif  
  PHIST_CHK_GERR(ghost_gemm(C, V, trans,W, (char*)"N", (void*)&alpha, (void*)&beta, GHOST_GEMM_NO_REDUCE,GHOST_GEMM_DEFAULT),*iflag);
PHIST_TASK_END(iflag);
  }

//! n x m serial dense matrix times k x m conj. transposed serial dense matrix gives m x k sdMat,
//! C=alpha*V*W + beta*C (serial XGEMM wrapper)
extern "C" void SUBR(sdMat_times_sdMatT)(_ST_ alpha, TYPE(const_sdMat_ptr) vV,
                                         TYPE(const_sdMat_ptr) vW,
                              _ST_ beta, TYPE(sdMat_ptr) vC,
                                         int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_SMALL;
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN_SMALLDETERMINISTIC(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,W,vW,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat_t,C,vC,*iflag);
#ifdef IS_COMPLEX
  char trans[]="C";
#else
  char trans[]="T";
#endif  
  PHIST_CHK_GERR(ghost_gemm(C, V, (char*)"N", W, trans, (void*)&alpha, (void*)&beta, GHOST_GEMM_NO_REDUCE,GHOST_GEMM_DEFAULT),*iflag);
PHIST_TASK_END(iflag);
  }


//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.   
//! Q is computed in place of V. If V does not have full rank, iflag>0   
//! indicates the dimension of the null-space of V. The first m-iflag    
//! columns of Q are an orthogonal basis of the column space of V, the  
//! remaining columns form a basis for the null space.  
extern "C" void SUBR(mvec_QR)(TYPE(mvec_ptr) vV, TYPE(sdMat_ptr) vR, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
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
    }
    *iflag=1-rank;
    return;
  }// case ncols=1: normalize single vector

  PHIST_DEB("mvec_QR: multi-vector case\n");

#if defined(PHIST_HAVE_TEUCHOS)&&defined(PHIST_HAVE_KOKKOS)

#ifdef GHOST_HAVE_CUDA
  static int any_cuda=-1;
  if (any_cuda==-1)
  {
    ghost_type_t ghost_type;
    PHIST_CHK_GERR(ghost_type_get(&ghost_type),*iflag);
    int is_cuda=(ghost_type==GHOST_TYPE_CUDA)?1:0;
    MPI_Allreduce(&is_cuda,&any_cuda,1,MPI_INT,MPI_SUM,V->context->mpicomm);
  }
  if (any_cuda)
  {
    *iflag=PHIST_NOT_IMPLEMENTED;
    return;
  }
#endif

  //TSQR for row major storage not available yet
  bool transV=(V->traits.storage==GHOST_DENSEMAT_ROWMAJOR);
  bool transR=(R->traits.storage==GHOST_DENSEMAT_ROWMAJOR);

  if( transR )
  {
    PHIST_CHK_IERR(*iflag=PHIST_NOT_IMPLEMENTED,*iflag);
  }
    
  // Here the actual TSQR call with col-major V and R begins...
  
  // TSQR does not do in-place QR, normalize copies the vector. So instead
  // we copy it ourelves (and memtranspose if in row-major order) and then
  // call the underlying normalizeOutOfPlace function directly.
  ghost_densemat_t* Vcopy=NULL, *Qcopy=NULL;
  ghost_densemat_traits_t vtraits = V->traits;
    vtraits.storage=GHOST_DENSEMAT_COLMAJOR;  
    vtraits.flags = (ghost_densemat_flags_t)((int)vtraits.flags & ~(int)GHOST_DENSEMAT_VIEW);
    //vtraits.flags = (ghost_densemat_flags_t)((int)vtraits.flags & ~(int)GHOST_DENSEMAT_PERMUTED); // melven: do we need this?
    vtraits.ncolsorig=vtraits.ncols;
    vtraits.nrowsorig=vtraits.nrows;
  ghost_densemat_create(&Vcopy,V->context,vtraits);
  ghost_densemat_create(&Qcopy,V->context,vtraits);
      
  {
    PHIST_ENTER_FCN("TSQR_memtranspose");
    // this allocates the memory for the vector, copies and memTransposes the data
    PHIST_CHK_GERR(Vcopy->fromVec(Vcopy,V,0,0),*iflag);
    // this allocates the result vector
    ST zero = st::zero();
    PHIST_CHK_GERR(Qcopy->fromScalar(Qcopy,&zero),*iflag);
  }

  // wrapper class for ghost_densemat_t for calling Belos.
  // The wrapper does not own the vector so it doesn't destroy it.
  phist::GhostMV mv_V(Vcopy,false);
  phist::GhostMV mv_Q(Qcopy,false);
    
  int nrows = R->traits.nrows;
  ncols = R->traits.ncols;
    
  PHIST_CHK_IERR(*iflag=nrows-ncols,*iflag);
  PHIST_CHK_IERR(*iflag=nrows-(V->traits.ncols),*iflag);

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

  {
    PHIST_ENTER_FCN("TSQR");
    PHIST_TRY_CATCH(rank = tsqr.normalizeOutOfPlace(mv_V,mv_Q,R_view),*iflag);
  }
  PHIST_DEB("V has %d columns and rank %d\n",ncols,rank);
  {
    PHIST_ENTER_FCN("TSQR_memtranspose");
    // copy (and memTranspose back if necessary)
    PHIST_CHK_GERR(V->fromVec(V,Qcopy,0,0),*iflag);
  }
  Vcopy->destroy(Vcopy);
  Qcopy->destroy(Qcopy);
  *iflag = ncols-rank;// return positive number if rank not full.
#else
  *iflag=PHIST_NOT_IMPLEMENTED; // no Trilinos, no TSQR, no mvec_QR (right now)
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
  *iflag=PHIST_NOT_IMPLEMENTED;
}
# else
extern "C" void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, Smvec_t* reV, Smvec_t* imV, int *iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}
# endif
#endif

void SUBR(sparseMat_create_fromRowFunc)(TYPE(sparseMat_ptr) *vA, const_comm_ptr_t vcomm,
        gidx_t nrows, gidx_t ncols, lidx_t maxnne,
                phist_sparseMat_rowFunc rowFunPtr, void* last_arg, int *iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  bool repart = *iflag&PHIST_SPARSEMAT_REPARTITION;
  bool d2clr  = *iflag&PHIST_SPARSEMAT_DIST2_COLOR;
  int outlev = *iflag&PHIST_SPARSEMAT_QUIET ? PHIST_DEBUG : PHIST_INFO;

  int sellC, sellSigma;
  get_C_sigma(&sellC,&sellSigma,*iflag, *((MPI_Comm*)vcomm));
  PHIST_SOUT(outlev, "Creating sparseMat with SELL-%d-%d format.\n", sellC, sellSigma);

  *iflag=0;
  
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  ghost_sparsemat_t* mat = NULL;
  ghost_context_t *ctx = NULL;
  PHIST_CAST_PTR_FROM_VOID(const MPI_Comm, comm, vcomm, *iflag);

  ghost_sparsemat_traits_t mtraits=(ghost_sparsemat_traits_t)GHOST_SPARSEMAT_TRAITS_INITIALIZER;
  ghost_sparsemat_flags_t flags=GHOST_SPARSEMAT_DEFAULT;

        mtraits.C = sellC;
        mtraits.sortScope = sellSigma;
        if (mtraits.sortScope > 1) {
            flags=(ghost_sparsemat_flags_t)(flags|GHOST_SPARSEMAT_PERMUTE);
        }

        if (repart)
        {
#ifdef USE_SCOTCH
          PHIST_SOUT(outlev, "Trying to repartition the matrix with SCOTCH\n");
          flags = (ghost_sparsemat_flags_t)(flags|GHOST_SPARSEMAT_SCOTCHIFY);
#else
          PHIST_SOUT(outlev, "SCOTCH not available, no matrix repartitioning\n");
#endif
        }
        else
        {
          PHIST_SOUT(outlev, "No matrix repartitioning requested\n");
        }
        mtraits.datatype = st::ghost_dt;
        mtraits.flags = flags;
  PHIST_CHK_GERR(ghost_context_create(&ctx,nrows,ncols,
        GHOST_CONTEXT_DEFAULT,NULL,GHOST_SPARSEMAT_SRC_FUNC,*comm,repart ? 0. : 1.),*iflag);
  PHIST_CHK_GERR(ghost_sparsemat_create(&mat,ctx,&mtraits,1),*iflag);                               

  ghost_sparsemat_src_rowfunc_t src = GHOST_SPARSEMAT_SRC_ROWFUNC_INITIALIZER;
  src.func = rowFunPtr;
  src.maxrowlen = maxnne;
      
  PHIST_CHK_GERR(mat->fromRowFunc(mat,&src),*iflag);
//#if PHIST_OUTLEV >= PHIST_VERBOSE
  char *str;
  ghost_context_string(&str,ctx);
  PHIST_SOUT(outlev,"%s\n",str);
  free(str); str = NULL;
  ghost_sparsemat_string(&str,mat);
  PHIST_SOUT(outlev,"%s\n",str);
  free(str); str = NULL;
//#endif
  *vA = (TYPE(sparseMat_ptr))mat;
PHIST_TASK_END(iflag);

  return;
}

