#include "phist_config.h"
#include <ghost/config.h>
#include <ghost/util.h>

#ifdef TEST_MVEC_MAPS_SAME
#undef TEST_MVEC_MAPS_SAME
#endif

#ifdef PHIST_TESTING_____
#define TEST_MVEC_MAPS_SAME(_v1,_v2,_iflag) \
if (_v1!=NULL && _v2!=NULL) { \
phist_const_map_ptr map1, map2; \
PHIST_CHK_IERR(SUBR(mvec_get_map)(_v1,&map1,_iflag),*_iflag); \
PHIST_CHK_IERR(SUBR(mvec_get_map)(_v2,&map2,_iflag),*_iflag); \
PHIST_CHK_IERR(phist_maps_compatible(map1,map2,_iflag),*_iflag); \
}
#else
#define TEST_MVEC_MAPS_SAME(_v1,_v2,_iflag)
#endif

/* helper macro to temporarily set the densemat's location to HOST and store the
   original value of the flag.
 */
#ifndef TMP_SET_DENSEMAT_LOCATION
#define TMP_SET_DENSEMAT_LOCATION(_void_ptr,_densemat_ptr,_orig_location) \
ghost_densemat* _densemat_ptr = (ghost_densemat*)(_void_ptr); \
ghost_location _orig_location = _densemat_ptr->traits.location; \
if (*iflag&PHIST_SDMAT_RUN_ON_HOST && false) _densemat_ptr->traits.location=GHOST_LOCATION_HOST;
#endif

#if defined(PHIST_HAVE_TEUCHOS)&&defined(PHIST_HAVE_KOKKOS)
template<>
Teuchos::RCP<node_type> ghost::TsqrAdaptor< _ST_ >::node_=Teuchos::null;
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
extern "C" void SUBR(sparseMat_read_mm)(TYPE(sparseMat_ptr)* vA, phist_const_comm_ptr vcomm,
const char* filename,int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  int iflag_in=*iflag;
  int outlev = *iflag&PHIST_SPARSEMAT_QUIET ? PHIST_DEBUG : PHIST_INFO;

  int sellC, sellSigma;
  phist::ghost_internal::get_C_sigma(&sellC,&sellSigma,*iflag, *((MPI_Comm*)vcomm));
  PHIST_SOUT(PHIST_VERBOSE, "Creating sparseMat with SELL-%d-%d format.\n", sellC, sellSigma);

  *iflag=0;

PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(const MPI_Comm,comm,vcomm,*iflag);
  if (filename==NULL)
  {
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
  ghost_sparsemat* mat;

  ghost_sparsemat_traits mtraits=(ghost_sparsemat_traits)GHOST_SPARSEMAT_TRAITS_INITIALIZER;
  ghost_sparsemat_flags flags=GHOST_SPARSEMAT_DEFAULT;

        mtraits.C = sellC;
        mtraits.sortScope = sellSigma;
        if (mtraits.sortScope > 1) 
        {
            flags=(ghost_sparsemat_flags)(flags|GHOST_SPARSEMAT_PERMUTE);
        }
        flags = (ghost_sparsemat_flags)(flags|get_perm_flag(iflag_in,outlev));

        mtraits.datatype = st::ghost_dt;
        mtraits.flags = flags;
        char* cfname=const_cast<char*>(filename);

//  PHIST_CHK_IERR(phist::ghost_internal::context_create(&ctx,0,0,
//        GHOST_CONTEXT_DEFAULT,cfname,GHOST_SPARSEMAT_SRC_MM,*comm,get_proc_weight(),iflag),*iflag);
  PHIST_CHK_GERR(ghost_sparsemat_create(&mat,NULL,&mtraits,1),*iflag);                               
  PHIST_CHK_GERR(ghost_sparsemat_init_mm(mat,cfname,*comm,get_proc_weight()),*iflag);
  char *str;
  ghost_context_string(&str,mat->context);
  PHIST_SOUT(outlev,"%s\n",str);
  free(str); str = NULL;
  ghost_sparsemat_info_string(&str,mat);
  PHIST_SOUT(outlev,"%s\n",str);
  free(str); str = NULL;
//#endif
  *vA = (TYPE(sparseMat_ptr))mat;
PHIST_TASK_END(iflag);
}

//! read a matrix from a Ghost CRS (binary) file.
extern "C" void SUBR(sparseMat_read_bin)(TYPE(sparseMat_ptr)* vA, phist_const_comm_ptr vcomm,
const char* filename,int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  int iflag_in=*iflag;
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

  ghost_sparsemat* mat;

  ghost_sparsemat_traits mtraits=(ghost_sparsemat_traits)GHOST_SPARSEMAT_TRAITS_INITIALIZER;
  ghost_sparsemat_flags flags=GHOST_SPARSEMAT_DEFAULT;

        mtraits.C = sellC;
        mtraits.sortScope = sellSigma;
        if (mtraits.sortScope > 1) {
            flags=(ghost_sparsemat_flags)(flags|GHOST_SPARSEMAT_PERMUTE);
        }

        flags = (ghost_sparsemat_flags)(flags|get_perm_flag(iflag_in,outlev));
        mtraits.datatype = st::ghost_dt;
        mtraits.flags = flags;
        char* cfname=const_cast<char*>(filename);

  PHIST_CHK_GERR(ghost_sparsemat_create(&mat,NULL,&mtraits,1),*iflag);                               
  PHIST_CHK_GERR(ghost_sparsemat_init_bin(mat,cfname,*comm,get_proc_weight()),*iflag);
//#if PHIST_OUTLEV >= PHIST_VERBOSE
  char *str;
  ghost_context_string(&str,mat->context);
  PHIST_SOUT(outlev,"%s\n",str);
  free(str); str = NULL;
  ghost_sparsemat_info_string(&str,mat);
  PHIST_SOUT(outlev,"%s\n",str);
  free(str); str = NULL;
//#endif
  *vA = (TYPE(sparseMat_ptr))mat;
PHIST_TASK_END(iflag);
}

//! read a matrix from a Harwell-Boeing (HB) file
extern "C" void SUBR(sparseMat_read_hb)(TYPE(sparseMat_ptr)* vA, phist_const_comm_ptr vcomm,
const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_TOUCH(vA);
  PHIST_TOUCH(vcomm);
  PHIST_TOUCH(filename);
  *iflag = PHIST_NOT_IMPLEMENTED; // not implemented in ghost, use converter script to bin crs
}

extern "C" void SUBR(sparseMat_read_mm_with_context)(TYPE(sparseMat_ptr)* A, phist_const_context_ptr ctx,
        const char* filename,int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
  return;
}

extern "C" void SUBR(sparseMat_read_bin_with_context)(TYPE(sparseMat_ptr)* A, phist_const_context_ptr ctx,
        const char* filename,int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_hb_with_context)(TYPE(sparseMat_ptr)* A, phist_const_ctx_ptr ctx,
        const char* filename,int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

//!@}

//! \name get information about the data distribution in a matrix (maps)

//!@{

//! get the complete context of the matrix
extern "C" void SUBR(sparseMat_get_context)(TYPE(const_sparseMat_ptr) vA, phist_const_context_ptr* vctx, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_sparsemat,A,vA,*iflag);
  *vctx = (phist_const_context_ptr)(A->context);
}

//! get the row distribution of the matrix
extern "C" void SUBR(sparseMat_get_row_map)(TYPE(const_sparseMat_ptr) vA, phist_const_map_ptr* vmap, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_sparsemat,A,vA,*iflag);
  *vmap = (phist_const_map_ptr)(A->context->row_map);
}

//! get column distribution of a matrix
//! we currently treat all maps as the same as we don't allow any fancy
//! operations using them anyway and ghost can handle both halo'd (colmap)
//! and standard (rowmap) vectors in the mvm.
extern "C" void SUBR(sparseMat_get_col_map)(TYPE(const_sparseMat_ptr) vA, phist_const_map_ptr* vmap, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_sparsemat,A,vA,*iflag);
  *vmap = (phist_const_map_ptr)(A->context->col_map);
}

//! get the map for vectors x in y=A*x
//! we currently treat all maps as the same as we don't allow any fancy
//! operations using them anyway and ghost can handle both halo'd (colmap)
//! and standard (rowmap) vectors in the mvm.
extern "C" void SUBR(sparseMat_get_domain_map)(TYPE(const_sparseMat_ptr) vA, phist_const_map_ptr* vmap, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const ghost_sparsemat,A,vA,*iflag);
  *vmap = (phist_const_map_ptr)(A->context->col_map);
}

//! get the map for vectors y in y=A*x
//! we currently treat all maps as the same as we don't allow any fancy
//! operations using them anyway and ghost can handle both halo'd (colmap)
//! and standard (rowmap) vectors in the mvm.
extern "C" void SUBR(sparseMat_get_range_map)(TYPE(const_sparseMat_ptr) vA, phist_const_map_ptr* vmap, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const ghost_sparsemat,A,vA,*iflag);
  *vmap = (phist_const_map_ptr)(A->context->col_map);
}
//@}

//! \name constructors

//@{
//! create a block-vector. The entries are stored contiguously
//! at val in column major ordering.
extern "C" void SUBR(mvec_create)(TYPE(mvec_ptr)* vV, 
        phist_const_map_ptr vmap, int nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  bool replicate_cuda_mem=*iflag&PHIST_MVEC_REPLICATE_DEVICE_MEM;
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_MVEC_CREATE(vmap,nvec,iflag);
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(const ghost_map, map,vmap,*iflag);
  ghost_densemat* result;
  ghost_densemat_traits vtraits = phist::ghost_internal::default_vtraits();
        vtraits.ncols=nvec;
        vtraits.ncolspadded=0;
        vtraits.datatype = st::ghost_dt;
        vtraits.flags = (ghost_densemat_flags)(vtraits.flags & ~GHOST_DENSEMAT_VIEW);

  // on CUDA nodes, allocate only device memory for mvecs
  ghost_type ghost_type;
  PHIST_CHK_GERR(ghost_type_get(&ghost_type),*iflag);
  if (ghost_type == GHOST_TYPE_CUDA) 
  {
    if (replicate_cuda_mem)
    {
      // allocate both host and device side (useful for testing)
      vtraits.location = (ghost_location)(GHOST_LOCATION_HOST|GHOST_LOCATION_DEVICE);
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


  PHIST_CHK_GERR(ghost_densemat_create(&result,(ghost_map*)map,vtraits),*iflag);
  ST zero = st::zero();
  // this allocates the vector and fills it with zeros
  PHIST_CHK_GERR(ghost_densemat_init_val(result,&zero),*iflag);
  PHIST_DEB("mvec nrows: %" PRlidx "\n",result->map->dim);
  *vV=(TYPE(mvec_ptr))(result);
PHIST_TASK_END(iflag);
}

//! create a block-vector as view of raw data. The map tells the object
//! how many rows it should 'see' in the data (at most lda, the leading
//! dimension of the 2D array values). CAVEAT: This function only works
//! if nrowshalo==nrowspadded in the map, which is in general only the case for
//! if there is only 1 MPI process or the matrix is trivially parallel.
extern "C" void SUBR(mvec_create_view)(TYPE(mvec_ptr)* vV, phist_const_map_ptr vmap, 
        _ST_* values, phist_lidx lda, int nvec,
        int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_CAST_PTR_FROM_VOID(const ghost_map, map,vmap,*iflag);
  ghost_densemat* result;
  ghost_densemat_traits vtraits = phist::ghost_internal::default_vtraits();
        vtraits.flags|=GHOST_DENSEMAT_VIEW;
        vtraits.ncols=nvec;
        vtraits.datatype = st::ghost_dt;

  PHIST_CHK_GERR(ghost_densemat_create(&result,(ghost_map*)map,vtraits),*iflag);

#ifdef PHIST_MVECS_ROW_MAJOR
  if (result->map->nhalo)
  {
    PHIST_OUT(PHIST_ERROR,"viewing plain data as row-major ghost_densemat only works \n"
                          "for vectors without communciation buffers (for spMVM)\n");
    *iflag=-1;
    return;
  }
#else
  if ((result->map->nhalo>lda))
  {
    PHIST_OUT(PHIST_ERROR,"viewing plain data as ghost_densemat only works \n"
                          "if the given lda can accomodate the required comm buffer of the vector!\n"
                          "nrows=%" PRlidx ", nrowshalopadded=%" PRlidx ", lda=%" PRlidx "\n",
        result->map->dim,result->map->dimhalopadded,lda);
    *iflag=-1;
    return;
  }
#endif
  PHIST_CHK_GERR(ghost_densemat_view_plain(result,(void*)values,lda),*iflag);
  *vV=(TYPE(mvec_ptr))(result);
  return;
}


//! create a serial dense n x m matrix on all procs, with column major
//! ordering.
extern "C" void SUBR(sdMat_create)(TYPE(sdMat_ptr)* vM, int nrows, int ncols, 
        phist_const_comm_ptr vcomm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_SMALL;
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN_SMALLDETERMINISTIC(ComputeTask)
  ghost_densemat* result;
  
  ghost_densemat_traits dmtraits=GHOST_DENSEMAT_TRAITS_INITIALIZER;
        dmtraits.ncols=(ghost_lidx)ncols;
        dmtraits.datatype=st::ghost_dt;
#ifdef PHIST_SDMATS_ROW_MAJOR
        dmtraits.storage=GHOST_DENSEMAT_ROWMAJOR;
#else
        dmtraits.storage=GHOST_DENSEMAT_COLMAJOR;
#endif
  // on CUDA nodes, allocate both host and device memory for sdMats
  ghost_type ghost_type;
  PHIST_CHK_GERR(ghost_type_get(&ghost_type),*iflag);
  if (ghost_type == GHOST_TYPE_CUDA) 
  {
    dmtraits.location = (ghost_location)(GHOST_LOCATION_HOST|GHOST_LOCATION_DEVICE);
  } 
  else 
  {
    dmtraits.location = GHOST_LOCATION_HOST;
  }

  // I think the sdMat should not have a context
  MPI_Comm comm = MPI_COMM_SELF;
  if (vcomm!=NULL) comm=*((MPI_Comm*)vcomm);
  ghost_map *map=ghost_map_create_light(nrows,comm);
  ghost_densemat_create(&result,map,dmtraits);
  // the map is created in this context, so we need to
  // release memory ownership here. The densemat has
  // taken ownership of the new map.
  map->ref_count--;
  ST zero = st::zero();
  PHIST_CHK_GERR(ghost_densemat_init_val(result,&zero),*iflag);
  *vM=(TYPE(sdMat_ptr))result;
PHIST_TASK_END(iflag);
}

void SUBR(sdMat_create_view)(TYPE(sdMat_ptr)* M, phist_const_comm_ptr comm,
        _ST_* values, phist_lidx lda, int nrows, int ncols,
        int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

//@}

//! retrieve the map of the vectors in V
extern "C" void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) vV, phist_const_map_ptr* vmap, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_densemat,V,vV,*iflag);
  PHIST_CHK_IERR(*iflag= (V->map==NULL)? PHIST_BAD_CAST:0, *iflag);
  *vmap=(phist_const_map_ptr)(V->map);
}

//! retrieve number of vectors/columns in V
extern "C" void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) vV, int* nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag = 0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_densemat,V,vV,*iflag);
  PHIST_CHK_IERR(*iflag=check_local_size(V->traits.ncols),*iflag);
  *nvec = (int)(V->traits.ncols);
}

//! get number of rows in local dense matrix
extern "C" void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) vM, int* nrows, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_densemat,M,vM,*iflag);
  PHIST_CHK_IERR(*iflag=check_local_size(M->map->dim),*iflag);
  *nrows = (int)(M->map->dim);
}
  
//! get number of cols in local dense matrix
extern "C" void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) vM, int* ncols, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(const ghost_densemat,M,vM,*iflag);
  *ncols = (int)(M->traits.ncols);
}


extern "C" void SUBR(mvec_extract_view)(TYPE(mvec_ptr) vV, _ST_** val, phist_lidx* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"

  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V, vV, *iflag);
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
      ghost_densemat_download(V);
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

extern "C" void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) vM, _ST_** val, phist_lidx* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,M, vM, *iflag);

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
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V, vV, *iflag);
/*
  PHIST_SOUT(PHIST_DEBUG,"ghost densemat upload\n"
                         "nrows=%" PRlidx ", ncols=%" PRlidx "\n"
                         "nrowshalo=%" PRlidx "\n"
                         "nrowspadded=%" PRlidx ", ncolspadded=%" PRlidx "\n",
                         V->traits.nrows, V->traits.ncols, 
                         V->traits.nrowshalo,
                         V->traits.nrowspadded, V->traits.ncolspadded);
  PHIST_SOUT(PHIST_DEBUG,"V flags: %d\n",(int)V->traits.flags);
*/
  PHIST_CHK_GERR(ghost_densemat_upload(V),*iflag);
#endif
}

extern "C" void SUBR(mvec_from_device)(TYPE(mvec_ptr) vV, int* iflag)
{
  *iflag=0;
#ifdef GHOST_HAVE_CUDA
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  //PHIST_PERFCHECK_VERIFY_FROM_DEVICE(vV,iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V, vV, *iflag);
  ghost_type ghost_type;
  PHIST_CHK_GERR(ghost_type_get(&ghost_type),*iflag);
  if (ghost_type == GHOST_TYPE_CUDA) 
  {
    PHIST_CHK_GERR(ghost_densemat_download(V),*iflag);
  }
#endif
}

extern "C" void SUBR(sdMat_to_device)(TYPE(sdMat_ptr) vM, int* iflag)
{
  *iflag=0;
#ifdef GHOST_HAVE_CUDA
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,M, vM, *iflag);
  ghost_type ghost_type;
  PHIST_CHK_GERR(ghost_type_get(&ghost_type),*iflag);
  if (ghost_type == GHOST_TYPE_CUDA) 
  {
    PHIST_CHK_GERR(ghost_densemat_upload(M),*iflag);
  }
#endif
}

extern "C" void SUBR(sdMat_from_device)(TYPE(sdMat_ptr) vM, int* iflag)
{
  *iflag=0;
#ifdef GHOST_HAVE_CUDA
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,M, vM, *iflag);
  PHIST_CHK_GERR(ghost_densemat_download(M),*iflag);
#endif
}

extern "C" void SUBR(mvec_to_mvec)(TYPE(const_mvec_ptr) v_in, TYPE(mvec_ptr) v_out, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V_in,v_in,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V_out,v_out,*iflag);

#ifdef TESTING
  // ask the map guru if the operation makes sense at all. We only do this in TESTING mode
  // because maps_compatible may do global reductions etc. (although it tries to avoid this)
  phist_const_map_ptr map_in, map_out;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(v_in,&map_in,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_get_map)(v_out,&map_out,iflag),*iflag);
  PHIST_CHK_NEG_IERR(phist_maps_compatible(map_in,map_out,iflag),*iflag);
#endif

  // set some convenient pointers
  ghost_map *map_in=V_in->map;
  ghost_map *map_out=V_out->map;
  
  ghost_densemat_flags flags_in=V_in->traits.flags;
  ghost_densemat_flags flags_out=V_out->traits.flags;
  
  ghost_lidx *lperm_in=V_in->map->loc_perm;
  ghost_lidx *lperm_out=V_out->map->loc_perm;

  ghost_gidx *gperm_in=V_in->map->glb_perm;
  ghost_gidx *gperm_out=V_out->map->glb_perm;
 
  // it seems like GHOST does not set these flags at all...
  // TODO: check this!!!
  bool outputPermuted=flags_out&GHOST_DENSEMAT_PERMUTED;
  bool  inputPermuted=flags_in &GHOST_DENSEMAT_PERMUTED;

  outputPermuted=(lperm_out!=NULL || gperm_out!=NULL);
  inputPermuted= (lperm_in!=NULL || gperm_in!=NULL);
  
  bool same_lperm = (lperm_in==lperm_out);
  bool same_gperm = (gperm_in==gperm_out);
  
  // if both are permuted with the same permutation, just copy
  bool no_perm_needed = (outputPermuted==inputPermuted && 
                         same_lperm && same_gperm);
  
  int me,nrank;
  ghost_rank(&me, V_in->map->mpicomm);
  ghost_nrank(&nrank,V_in->map->mpicomm);
  
  // check if the two maps have the same number of local elements on each process
  bool compatible_maps=(map_in==map_out);
  if (!compatible_maps)
  {
    compatible_maps=(map_in->mpicomm==map_out->mpicomm);
    for (int i=0; i<nrank;i++)
    {
      PHIST_SOUT(PHIST_DEBUG,"lnrows PE%d, v_in=%d, v_out=%d\n",i,map_in->ldim[i],map_out->ldim[i]); 
      compatible_maps &= (map_in->ldim[i]==map_out->ldim[i]);
    }
  }
  
  if (compatible_maps==false)
  {
    // I think ghost_densemat_init_densemat only works if both have the same map and does not call any
    // permute function. So what we'll do is temporarily replace the map of the output vec-
    // tor, copy the data by ghost_densemat_init_densemat, and then perform the permutations. This only   
    // works if the number of local elements is the same in both maps. We could also work  
    // with a temporary vector, but I think GHOST currently does not allow Zoltan/Scotch to    
    // change the number of local elements.
    PHIST_SOUT(PHIST_ERROR,"phist/ghost: mvec_to_mvec only implemented if both vectors have\n"
                           "the same local lengths on all MPI ranks (and the same MPI communicator)\n");
    *iflag=-1;
    return;
  }

// macro to perform a ghost call on a vector with a different map
#define GHOST_FUNC_MAP(_vec,_map,_permflag,_permmethod,_lperm,_gperm,_fnc,...) \
{\
  ghost_error _gerr=GHOST_SUCCESS; \
  ghost_map *_orig_map=(_vec)->map; \
  ghost_densemat_flags _orig_flags = (_vec)->traits.flags; \
  ghost_maptype _orig_permuted = _vec->map->type; \
  ghost_lidx *_orig_lperm=(_vec)->map->loc_perm; \
  ghost_gidx *_orig_gperm=(_vec)->map->glb_perm; \
  (_vec)->map=(_map);\
  (_vec)->traits.flags = (ghost_densemat_flags)( (_permflag&GHOST_DENSEMAT_PERMUTED) ? GHOST_DENSEMAT_PERMUTED|(_vec)->traits.flags : (~GHOST_DENSEMAT_PERMUTED)&(_vec)->traits.flags ); \
  _gerr=_fnc(__VA_ARGS__);\
  (_vec)->map=_orig_map;\
  (_vec)->traits.flags = _orig_flags; \
  PHIST_CHK_GERR(_gerr,*iflag); \
}

  GHOST_FUNC_MAP(V_out,map_in,V_in->traits.flags,V_in->map->type,lperm_in,gperm_in,ghost_densemat_init_densemat,V_out,V_in,0,0);

  if (no_perm_needed) return;
  
  if (inputPermuted)
  {
      PHIST_SOUT(PHIST_DEBUG,"mvec_to_mvec: unpermute input vector\n");
      GHOST_FUNC_MAP(V_out,map_in,V_in->traits.flags,V_in->map->type,lperm_in,gperm_in,ghost_densemat_permute,V_out,GHOST_PERMUTATION_PERM2ORIG);
  }
  
  if (outputPermuted)
  {
    PHIST_SOUT(PHIST_DEBUG,"mvec_to_mvec: permute output vector\n");
    PHIST_CHK_GERR(ghost_densemat_permute(V_out,GHOST_PERMUTATION_ORIG2PERM),*iflag);
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
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);
  ghost_densemat *Vblock=(ghost_densemat*)(*vVblock);
  
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

    ghost_densemat_destroy(Vblock);
    Vblock=NULL;
  }
  PHIST_CHK_GERR(ghost_densemat_create_and_view_densemat_cols(&Vblock, V, jmax-jmin+1, jmin),*iflag);

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
  TEST_MVEC_MAPS_SAME(vV,vVblock,iflag)
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,Vblock,vVblock,*iflag);
  *iflag=0;
#ifdef PHIST_TESTING
// nonzero error code if #vectors in Vblock too small or large
  PHIST_CHK_IERR(*iflag=(jmax-jmin+1)-Vblock->traits.ncols,*iflag);
  PHIST_CHK_IERR(*iflag=((jmin<0)||(jmax>=V->traits.ncols))?PHIST_INVALID_INPUT:0,*iflag);
  // The phist kernel interface requires Vblock to be created by mvec_create, 
  // which calls ghost_densemat_init_val to allocate the block of memory. So we perform a few
  // safety checks here
  PHIST_CHK_IERR(*iflag=(Vblock->map->dim==V->map->dim)?0:PHIST_INVALID_INPUT,*iflag)
  // not sure what the ghost function ghost_densemat_init_densemat actually supports, but I think 
  // this makes sense:
  //PHIST_CHK_IERR(*iflag=(Vblock->traits.nrowspadded==V->traits.nrowspadded)?0:PHIST_INVALID_INPUT,*iflag)
#else
  PHIST_TOUCH(jmax);
#endif  
  PHIST_CHK_GERR(ghost_densemat_init_densemat(Vblock,V,(ghost_lidx)0,(ghost_lidx)jmin),*iflag);
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
  TEST_MVEC_MAPS_SAME(vV,vVblock,iflag)
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,Vblock,vVblock,*iflag);

#ifdef PHIST_TESTING
int nv_v,nv_vb;
phist_lidx nr_v,nr_vb;
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

  // create a view of the requested columns of V
  ghost_densemat *Vcols=NULL;
  ghost_densemat_create_and_view_densemat_cols(&Vcols,V,(ghost_lidx)(jmax-jmin+1),(ghost_lidx)jmin);

  // copy the data
  PHIST_CHK_GERR(ghost_densemat_init_densemat(Vcols,Vblock,0,0),*iflag);
  // delete the view
  ghost_densemat_destroy(Vcols);
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
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,M,vM,*iflag);

#ifdef PHIST_TESTING
  if (imin<0||jmin<0||imax>M->map->dim||jmax>M->traits.ncols)
  {
    PHIST_OUT(PHIST_ERROR,"%s: range out of bounds of matrix\n"
                          "requested range: (%d:%d,%d:%d), matrix dim: %dx%d\n",
                          __FUNCTION__,
                          imin,imax,jmin,jmax,M->map->dim,M->traits.ncols);
  }
#endif
  //TODO: we only view the host side of the vector here, this function should
  //      eventually be moved into ghost and the accelerator stuff added.

  // first just create a view of the corresponding columns
  ghost_densemat *Mblock;
  ghost_densemat_create_and_view_densemat(&Mblock, M, imax-imin+1,imin,jmax-jmin+1, jmin);

  if (*vMblock!=NULL)
  {
    //PHIST_DEB("deleting previous object in %s\n",__FUNCTION__);
    PHIST_CAST_PTR_FROM_VOID(ghost_densemat,tmp,*vMblock,*iflag);
    ghost_densemat_destroy(tmp);
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
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,M,vM,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,Mblock,vMblock,*iflag);

#ifdef PHIST_TESTING
  int nr=imax-imin+1;
  int nc=jmax-jmin+1;
  if (Mblock->map->dim!=nr || Mblock->traits.ncols!=nc)
  {
    PHIST_SOUT(PHIST_ERROR,"result block has wrong dimensions %dx%d, "
                           "requested range is (%d:%d,%d:%d)\n",
                           Mblock->map->dim,Mblock->traits.ncols,imin,imax,jmin,jmax);
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
  if (imin<0 || imax>M->map->dim ||
      jmin<0 || jmax>M->traits.ncols )
  {
    PHIST_SOUT(PHIST_ERROR,"requested range invalid. M is %dx%d, "
                           "requested range is (%d:%d,%d:%d)\n",
                           M->map->dim,M->traits.ncols,imin,imax,jmin,jmax);
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
#endif 
  _ST_ *m_ptr, *mb_ptr;
  phist_lidx ldm, ldmb;
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
  PHIST_CHK_GERR(ghost_densemat_init_densemat(Mblock,M,imin,jmin),*iflag);
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
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,Mblock,vMblock,*iflag);

  ghost_densemat* Mb_view=NULL;
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(vM,(TYPE(sdMat_ptr)*)&Mb_view,imin,imax,jmin,jmax,iflag),*iflag);
  PHIST_CHK_GERR(ghost_densemat_init_densemat(Mb_view,Mblock,0,0),*iflag);
  ghost_densemat_destroy(Mb_view);
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
  PHIST_CAST_PTR_FROM_VOID(ghost_sparsemat,A,vA,*iflag);

  
  // delete the matrix data but not the context
  ghost_sparsemat_destroy(A);
}

//!
extern "C" void SUBR(mvec_delete)(TYPE(mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  if (vV==NULL) return;
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);

  ghost_densemat_destroy(V);
}

//!
extern "C" void SUBR(sdMat_delete)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  if (vM==NULL) return;
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,M,vM,*iflag);
  ghost_densemat_destroy(M);
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
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);
//  PHIST_DEB("put value, V @ %p. V->traits.nrows=%" PRlidx "\n",V,V->traits.nrows);
  ghost_densemat_init_val(V,(void*)&value);
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
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);
  PHIST_CHK_GERR(ghost_densemat_init_func(V,funPtr,last_arg),*iflag);
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
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);
  ghost_densemat_init_val(V,(void*)&value);
PHIST_TASK_END(iflag);
}

//! put scalar value into all elements of a multi-vector
extern "C" void SUBR(sdMat_identity)(TYPE(sdMat_ptr) V, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  bool host_only = (*iflag&PHIST_SDMAT_RUN_ON_HOST);
  *iflag = 0;
  PHIST_PERFCHECK_VERIFY_SMALL;

  _ST_ *V_raw = NULL;
  phist_lidx lda;
  int m, n;
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(V, &V_raw, &lda, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(V, &m, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(V, &n, iflag), *iflag);
  for(int i = 0; i < m; i++)
    for(int j = 0; j < n; j++)
      V_raw[lda*i+j] = (i==j) ? st::one() : st::zero();
  if (!host_only)
  {
    PHIST_CHK_IERR(SUBR(sdMat_to_device)(V,iflag),*iflag);
  }
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
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);
  V->fromRand(V);
PHIST_TASK_END(iflag);
}
#endif

extern "C" void SUBR(mvec_print)(TYPE(const_mvec_ptr) vV, int* iflag)
{
  *iflag = 0;
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);
  std::cout << "# local rows: "<<V->map->dim<<std::endl;
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
    ghost_densemat_string(&str,V);
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
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,M,vM,*iflag);
  std::cout << "# rows:       "<<M->map->dim<<std::endl;
  std::cout << "# cols:       "<<M->traits.ncols<<std::endl;
  std::cout << "# row major:  "<<(M->traits.storage & GHOST_DENSEMAT_ROWMAJOR)<<std::endl;
  std::cout << "# stride:     "<<M->stride<<std::endl;
  // always print the host side of the sdMat, if we don't replace the location by HOST in the traits
  // temporarily, GHOST will download the memory and allocate the host side.
  {
    char *str=NULL;
    ghost_location locM=M->traits.location;
    M->traits.location=GHOST_LOCATION_HOST;
    ghost_densemat_string(&str,M);
    M->traits.location=locM;
    std::cout << str <<std::endl;
    free(str); str = NULL;
  }
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
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,M,vM,*iflag);
  M->fromRand(M);
PHIST_TASK_END(iflag);
}
#endif

//! put random numbers into all elements of a serial dense matrix
extern "C" void SUBR(sdMat_sync_values)(TYPE(sdMat_ptr) vM, phist_const_comm_ptr vcomm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN_SMALLDETERMINISTIC(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,M,vM,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_mpi_comm,comm,vcomm,*iflag);
  ghost_densemat_sync_vals(M, *comm, 0);
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
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);  
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
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);  
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
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);  
  ghost_scale(V,(void*)&scalar);
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
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);  
  ghost_vscale(V,(void*)scalar);
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
  TEST_MVEC_MAPS_SAME(vX,vY,iflag)
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,X,vX,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,Y,vY,*iflag);
  ST a=alpha, b=beta;
  if (alpha==st::one() && beta==st::zero())
  {
      PHIST_DEB("copy Y=X\n");
    PHIST_CHK_GERR(ghost_densemat_init_densemat(Y,X,0,0),*iflag);
  }
  else if (alpha==st::zero())
  {
    if (beta!=st::one())
    {
      PHIST_DEB("scale output Y=beta*Y\n");
      PHIST_CHK_GERR(ghost_scale(Y,(void*)&b),*iflag);
    }
  }
  else if (beta==st::one())
  {
    PHIST_DEB("axpy operation: Y=alpha*X+Y\n");
    PHIST_CHK_GERR(ghost_axpy(Y,X,(void*)&a),*iflag);
  }
  else
  {
    PHIST_DEB("general case: Y=alpha*X+beta*Y\n");
    PHIST_CHK_GERR(ghost_axpby(Y,X,(void*)&a,(void*)&b),*iflag);
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
  TEST_MVEC_MAPS_SAME(vX,vY,iflag)
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,X,vX,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,Y,vY,*iflag);
  if(beta == st::one())
  {
    ghost_vaxpy(Y,X,(void*)alpha);
  }
  else
  {
    // ghost also expects a vector for beta, so construct one:
    int nvec = 0;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vY,&nvec,iflag),*iflag);
    _ST_ b[nvec];
    for(int i = 0; i < nvec; i++)
      b[i] = beta;

    ghost_vaxpby(Y,X,(void*)alpha,(void*)b);
  }
PHIST_TASK_END(iflag);
}

//! B=alpha*A+beta*B
extern "C" void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA,
                            _ST_ beta,  TYPE(sdMat_ptr)       vB,
                            int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  _ST_ a=alpha, b=beta;
  PHIST_PERFCHECK_VERIFY_SMALL;

PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN_SMALLDETERMINISTIC(ComputeTask)

  // if the user specifies PHIST_SDMAT_RUN_ON_HOST, manually switch the location setting
  // of the densemats and reset them after the call. The GHOST kernel is the same for sdMats
  // and mvecs, so we can simply call mvec_add_mvec.
  TMP_SET_DENSEMAT_LOCATION(vA,A,locA);
  TMP_SET_DENSEMAT_LOCATION(vB,B,locB);
    PHIST_CHK_GERR(ghost_axpby(B,A,(void*)&a,(void*)&b),*iflag);
  A->traits.location=locA;
  B->traits.location=locB;
PHIST_TASK_END(iflag);
}

//! B=alpha*A+beta*B
extern "C" void SUBR(sdMatT_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA,
                            _ST_ beta,  TYPE(sdMat_ptr)       vB,
                            int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  TMP_SET_DENSEMAT_LOCATION(vA,A,locA);
  TMP_SET_DENSEMAT_LOCATION(vB,B,locB);
  *iflag=0;
  // simple workaround
  TYPE(sdMat_ptr) I = NULL;
  int m = 0;
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(vA,&m,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_create)(&I,m,m,&A->map->mpicomm,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_identity)(I,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMatT_times_sdMat)(alpha,vA,I,beta,vB,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(I,iflag),*iflag);

  A->traits.location=locA;
  B->traits.location=locB;
}

//! spMVM communication
extern "C" void SUBR(sparseMat_times_mvec_communicate)(TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=0;

  PHIST_CAST_PTR_FROM_VOID(ghost_sparsemat,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,x,vx,*iflag);

#ifdef GHOST_HAVE_MPI
    ghost_densemat_halo_comm comm = GHOST_DENSEMAT_HALO_COMM_INITIALIZER;
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
    PHIST_CHK_GERR(ghost_densemat_halocomm_init(x,A->context,&comm),*iflag);
PHIST_TASK_END(iflag)
PHIST_TASK_POST_STEP(iflag)
    PHIST_CHK_GERR(ghost_densemat_halocomm_start(x,A->context,&comm),*iflag);
    PHIST_CHK_GERR(ghost_densemat_halocomm_finalize(x,A->context,&comm),*iflag);
#else
PHIST_TASK_POST_STEP(iflag)
#endif
}

//! y=alpha*A*x+beta*y.
extern "C" void SUBR(sparseMat_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, 
_ST_ beta, TYPE(mvec_ptr) vy, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(fused_spmv_mvdot_mvadd)(alpha,vA,vx,beta,vy,st::zero(),st::zero(),NULL,NULL,NULL,iflag),*iflag);
}

extern "C" void SUBR(fused_spmv_mvdot)(_ST_ alpha, TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, 
_ST_ beta, TYPE(mvec_ptr) vy, _ST_* ydoty, _ST_* xdoty, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(fused_spmv_mvdot_mvadd)(alpha,vA,vx,beta,vy,st::zero(),st::zero(),NULL,ydoty,xdoty,iflag),*iflag);
}

// This is the central place where we call the GHOST sparse matrix-vector product with all its bells and whistles
extern "C" void SUBR(fused_spmv_mvdot_mvadd)(_ST_ alpha, TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, 
_ST_ beta, TYPE(mvec_ptr) vy, 
_ST_ gamma, _ST_ delta, TYPE(mvec_ptr) vz,
_ST_* ydoty, _ST_* xdoty, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"

#ifdef PHIST_TESTING
  // check if the input and output have the correct maps, that is they live in the same index space
  // as the columns and rows of A, respectively, and have the same distribution and permutation.
  {
    phist_const_map_ptr map_x, map_y, range_map_A, domain_map_A;
    PHIST_CHK_IERR(SUBR(mvec_get_map)(vx,&map_x,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_get_map)(vy,&map_y,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sparseMat_get_range_map)(vA,&range_map_A,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sparseMat_get_domain_map)(vA,&domain_map_A,iflag),*iflag);
    
    // x and y must be correctly partitioned and permuted at this point, so demand *iflag=0 here:
    PHIST_CHK_IERR(phist_maps_compatible(map_x, domain_map_A,iflag),*iflag);
    PHIST_CHK_IERR(phist_maps_compatible(map_y, range_map_A,iflag),*iflag);
  }
#endif

  ghost_spmv_opts spMVM_opts=GHOST_SPMV_OPTS_INITIALIZER;
  // ghost spmvm mode
  if( *iflag & PHIST_SPMVM_ONLY_LOCAL )
    spMVM_opts.flags = (ghost_spmv_flags)((int)spMVM_opts.flags | (int)GHOST_SPMV_MODE_NOCOMM);
  else if( *iflag & PHIST_SPMVM_OVERLAP )
    spMVM_opts.flags = (ghost_spmv_flags)((int)spMVM_opts.flags | (int)GHOST_SPMV_MODE_OVERLAP);
  else if( *iflag & PHIST_SPMVM_TASK )
    spMVM_opts.flags = (ghost_spmv_flags)((int)spMVM_opts.flags | (int)GHOST_SPMV_MODE_TASK);
  bool no_reduce=false;
  if (*iflag & PHIST_NO_GLOBAL_REDUCTION )
  {
    spMVM_opts.flags = (ghost_spmv_flags)((int)spMVM_opts.flags | (int)GHOST_SPMV_NOT_REDUCE);
    no_reduce=true;
  }
    
  *iflag=0;

  PHIST_COUNT_MATVECS(vx);

  PHIST_CAST_PTR_FROM_VOID(ghost_sparsemat,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,x,vx,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,y,vy,*iflag);

  int nvec = 0;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vy,&nvec,iflag),*iflag);


  // this is not checked by maps_compatible because you can still add or dot-product
  // vectors with different number of halo elements. I'm not quite sure about the
  // difference between nrowshalo, nrowspadded and nrowshalopadded, I would think
  // that after the actual vector elements there's a padding, followed by the halo
  // and possibly more padding. However, in practice they are all set to the same
  // value, it seems, so I just check nrowshalo-nrows > halo_elements here for now:
/*if (A->context->halo_elements<=x->traits.nrowshalo-x->traits.nrowspadded == false)
{
  PHIST_OUT(PHIST_WARNING,"The following compatibility test fails: nrows=%d\nhalo_elements=%d\nnrowshalo=%d\nnrowspadded=%d\nnrowshalopadded=%d\nmaxnrowshalo=%d\n",
        x->traits.nrows,
        A->context->halo_elements,
        x->traits.nrowshalo,
        x->traits.nrowspadded,
        x->traits.nrowshalopadded,
        x->traits.maxnrowshalo);  
}
 PHIST_CHK_IERR(*iflag=(A->context->halo_elements<=x->traits.nrowshalo-x->traits.nrowspadded)?0:PHIST_INVALID_INPUT,*iflag);
*/
  // TODO: What was this test supposed to do?
  
  if (alpha==st::zero())
  {
    // no MVM needed
    if (beta==st::zero())
    {
      PHIST_CHK_IERR(SUBR(mvec_put_value)(vy,beta,iflag),*iflag);
      if( xdoty != NULL )
        for(int i = 0; i < nvec; i++)
          xdoty[i] = st::zero();
      if( ydoty != NULL )
        for(int i = 0; i < nvec; i++)
          ydoty[i] = st::zero();
    }
    else
    {
      if (beta!=st::one())
      {
        PHIST_CHK_IERR(SUBR(mvec_scale)(vy,beta,iflag),*iflag);
      }
      if( xdoty != NULL )
      {
        if (no_reduce) *iflag=PHIST_NO_GLOBAL_REDUCTION;
        PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(vx,vy,xdoty,iflag),*iflag);
      }
      if( ydoty != NULL )
      {
        if (no_reduce) *iflag=PHIST_NO_GLOBAL_REDUCTION;
        PHIST_CHK_IERR(SUBR(mvec_dot_mvec)(vy,vy,ydoty,iflag),*iflag);
      }
    }
  }
  else
  {
  PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)


    if (alpha!=st::one())
    {
      spMVM_opts.flags = (ghost_spmv_flags)((int)spMVM_opts.flags | (int)GHOST_SPMV_SCALE);
    }

    if (beta==st::one())
    {
      spMVM_opts.flags = (ghost_spmv_flags)((int)spMVM_opts.flags | (int)GHOST_SPMV_AXPY);
    }
    else if (beta!=st::zero())
    {
      spMVM_opts.flags = (ghost_spmv_flags)((int)spMVM_opts.flags | (int)GHOST_SPMV_AXPBY);
    }

    if( xdoty != NULL )
    {
      spMVM_opts.flags = (ghost_spmv_flags)((int)spMVM_opts.flags | (int)GHOST_SPMV_DOT_XY);
    }
    if( ydoty != NULL )
    {
      spMVM_opts.flags = (ghost_spmv_flags)((int)spMVM_opts.flags | (int)GHOST_SPMV_DOT_YY);
    }
    if( vz != NULL )
    {
      spMVM_opts.flags = (ghost_spmv_flags)((int)spMVM_opts.flags | (int)GHOST_SPMV_CHAIN_AXPBY);
    }

    std::vector<_ST_> dotBuff;
    if( xdoty != NULL || ydoty != NULL )
    {
      dotBuff.resize(3*nvec);
    }

    spMVM_opts.alpha = &alpha;
    spMVM_opts.beta = &beta;
    spMVM_opts.dot = &dotBuff[0];
    
    // ghost computes z = delta*z + eta*(spmv result),
    // so our delta is their eta and our gamma is their delta
    spMVM_opts.eta=&delta;
    spMVM_opts.delta=&gamma;
    spMVM_opts.z=(ghost_densemat*)vz;
    
    // ghost messes with the maps of the output vector: if it finds that it is the col_map of the matrix,
    // it simply swaps the pointer to the row_map. So we backup the pointers and reset them after the call
    ghost_map *map_x=x->map, *map_y=y->map;

    // call ghosts spMV
    PHIST_CHK_GERR(ghost_spmv(y,A,x,spMVM_opts),*iflag);

    x->map=map_x;
    y->map=map_y;
    
    // copy the dot results if necessary
    if (ydoty!=NULL)
    {
      for (int i=0;i<nvec;i++) ydoty[i]=dotBuff[i];
    }
    if (xdoty!=NULL)
    {
      for (int i=0;i<nvec;i++) xdoty[i]=dotBuff[nvec+i];
    }
    
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
  ghost_spmv_opts spMVM_opts=GHOST_SPMV_OPTS_INITIALIZER;
  // ghost spmvm mode
  if( *iflag & PHIST_SPMVM_ONLY_LOCAL )
    spMVM_opts.flags = (ghost_spmv_flags)((int)spMVM_opts.flags | (int)GHOST_SPMV_MODE_NOCOMM);
  else if( *iflag & PHIST_SPMVM_OVERLAP )
    spMVM_opts.flags = (ghost_spmv_flags)((int)spMVM_opts.flags | (int)GHOST_SPMV_MODE_OVERLAP);
  else if( *iflag & PHIST_SPMVM_TASK )
    spMVM_opts.flags = (ghost_spmv_flags)((int)spMVM_opts.flags | (int)GHOST_SPMV_MODE_TASK);
  *iflag=0;

  PHIST_COUNT_MATVECS(vx);

  PHIST_CAST_PTR_FROM_VOID(ghost_sparsemat,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,x,vx,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,y,vy,*iflag);
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

    spMVM_opts.flags = (ghost_spmv_flags)((int)spMVM_opts.flags | (int)GHOST_SPMV_VSHIFT);
    ST ghost_shifts[nvec];
    for (int i=0;i<nvec;i++) ghost_shifts[i]=-shifts[i];
    spMVM_opts.alpha = &alpha;
    spMVM_opts.beta = &beta;
    spMVM_opts.gamma = ghost_shifts;
    
    if (beta==st::one())
    {
      spMVM_opts.flags = (ghost_spmv_flags)((int)spMVM_opts.flags | (int)GHOST_SPMV_AXPY);
    }
    else if (beta!=st::zero())
    {
      spMVM_opts.flags = (ghost_spmv_flags)((int)spMVM_opts.flags | (int)GHOST_SPMV_AXPBY);
    }
    if (alpha!=st::one())
    {
      spMVM_opts.flags = (ghost_spmv_flags)((int)spMVM_opts.flags | (int)GHOST_SPMV_SCALE);
    }

    *iflag=ghost_spmv(y,A,x,spMVM_opts);
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
  TEST_MVEC_MAPS_SAME(vV,vW,iflag)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,W,vW,*iflag);
  // NOTE: calculate local dot by hand and do the reduction by hand
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  //ghost_dot(s,V,W);
  PHIST_CHK_GERR(ghost_localdot(s,V,W),*iflag);
PHIST_TASK_END(iflag);

PHIST_TASK_POST_STEP(iflag);

#ifdef GHOST_HAVE_MPI
  ghost_mpi_op sumOp;
  ghost_mpi_datatype mpiDt;
  ghost_mpi_op_sum(&sumOp,V->traits.datatype);
  ghost_mpi_datatype_get(&mpiDt,V->traits.datatype);
  PHIST_CHK_IERR(*iflag = MPI_Allreduce(MPI_IN_PLACE, s, V->traits.ncols, mpiDt, sumOp, V->map->mpicomm), *iflag);
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
  TEST_MVEC_MAPS_SAME(vV,vW,iflag)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,W,vW,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,C,vC,*iflag);
#ifdef IS_COMPLEX
  char trans[]="C";
#else
  char trans[]="T";
#endif  
  _ST_ mybeta = st::zero();
  phist_const_comm_ptr vcomm=NULL;
  phist_const_map_ptr map=NULL;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(W,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_comm(map,&vcomm,iflag),*iflag);
  PHIST_CAST_PTR_FROM_VOID(const MPI_Comm,comm,vcomm,*iflag);
  int rank = 0;
  PHIST_CHK_IERR(*iflag = MPI_Comm_rank(*comm,&rank),*iflag);
  if( rank == 0 )
  {
    mybeta = beta;
  }
  phist_lidx ncC = C->traits.ncols;
/*
  PHIST_DEB("VtV=C, V %" PRlidx "x%" PRlidx ", \n"
            "       W %" PRlidx "x%" PRlidx ", \n"
            "       C %" PRlidx "x%" PRlidx "\n", 
  V->traits.nrows,V->traits.ncols,
  W->traits.nrows,W->traits.ncols,
  C->traits.nrows,C->traits.ncols);
*/

  // NOTE: we call the allreduction by hand afterwards to allow asynchronuous communication!
  /*
  PHIST_CHK_GERR(ghost_gemm(C,V,trans,W,(char*)"N",(void*)&alpha,(void*)&beta,GHOST_GEMM_ALL_REDUCE,GHOST_GEMM_DEFAULT),*iflag);
  */
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  ghost_error gemm_err = ghost_gemm(C,V,trans,W,(char*)"N",(void*)&alpha,(void*)&mybeta,GHOST_GEMM_NO_REDUCE,GHOST_GEMM_DEFAULT);
  if( gemm_err == GHOST_ERR_NOT_IMPLEMENTED )
  {
    // copy result
    ghost_densemat* Ccopy=NULL;
    ghost_densemat_traits vtraits = C->traits;
    vtraits.storage=GHOST_DENSEMAT_ROWMAJOR;  
    vtraits.flags = (ghost_densemat_flags)((int)vtraits.flags & ~(int)GHOST_DENSEMAT_VIEW);
    ghost_densemat_create(&Ccopy,C->map,vtraits);

    // this allocates the memory for the vector, copies and memTransposes the data
    PHIST_CHK_GERR(ghost_densemat_init_densemat(Ccopy,C,0,0),*iflag);

    PHIST_CHK_GERR(gemm_err = ghost_gemm(Ccopy,V,trans,W,(char*)"N",(void*)&alpha,(void*)&mybeta,GHOST_GEMM_NO_REDUCE,GHOST_GEMM_DEFAULT),*iflag);

    // memtranspose data
    PHIST_CHK_GERR(ghost_densemat_init_densemat(C,Ccopy,0,0),*iflag);
    ghost_densemat_destroy(Ccopy);
  }
  PHIST_CHK_GERR(gemm_err,*iflag);
PHIST_TASK_END(iflag);

PHIST_TASK_POST_STEP(iflag);

  PHIST_CHK_GERR(ghost_densemat_reduce(C,GHOST_ALLREDUCE),*iflag);
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
  TEST_MVEC_MAPS_SAME(vV,vW,iflag)
  PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
    PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);
    PHIST_CAST_PTR_FROM_VOID(ghost_densemat,C,vC,*iflag);
    PHIST_CAST_PTR_FROM_VOID(ghost_densemat,W,vW,*iflag);

    phist_lidx nrV,nrW;
    int ncV, ncW, nrC, ncC;
    nrV=V->map->dim;  ncV=V->traits.ncols;
    nrW=W->map->dim;  ncW=V->traits.ncols;
    nrC=C->map->dim;  ncC=V->traits.ncols;

#ifdef PHIST_TESTING
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
    PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);
    PHIST_CAST_PTR_FROM_VOID(ghost_densemat,C,vC,*iflag);

#ifdef PHIST_TESTING
    int ncV, nrC, ncC;
    ncV=V->traits.ncols;
     nrC=C->map->dim;  ncC=V->traits.ncols;

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
  TMP_SET_DENSEMAT_LOCATION(vV,V,locV);
  TMP_SET_DENSEMAT_LOCATION(vW,W,locW);
  TMP_SET_DENSEMAT_LOCATION(vC,C,locC);
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_SMALL;
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN_SMALLDETERMINISTIC(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,W,vW,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,C,vC,*iflag);
  char trans[]="N";  
  PHIST_CHK_GERR(ghost_gemm(C,V,trans,W,trans,(void*)&alpha,(void*)&beta,GHOST_GEMM_NO_REDUCE,GHOST_GEMM_DEFAULT),*iflag);
PHIST_TASK_END(iflag);
  V->traits.location=locV;
  W->traits.location=locW;
  C->traits.location=locC;
}

//! n x m conj. transposed serial dense matrix times m x k serial dense matrix gives m x k sdMat,
//! C=alpha*V*W + beta*C (serial XGEMM wrapper)
extern "C" void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vV,
                                         TYPE(const_sdMat_ptr) vW,
                              _ST_ beta, TYPE(sdMat_ptr) vC,
                                         int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  TMP_SET_DENSEMAT_LOCATION(vV,V,locV);
  TMP_SET_DENSEMAT_LOCATION(vW,W,locW);
  TMP_SET_DENSEMAT_LOCATION(vC,C,locC);
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_SMALL;
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN_SMALLDETERMINISTIC(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,W,vW,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,C,vC,*iflag);
#ifdef IS_COMPLEX
  char trans[]="C";
#else
  char trans[]="T";
#endif  
  PHIST_CHK_GERR(ghost_gemm(C, V, trans,W, (char*)"N", (void*)&alpha, (void*)&beta, GHOST_GEMM_NO_REDUCE,GHOST_GEMM_DEFAULT),*iflag);
PHIST_TASK_END(iflag);
  V->traits.location=locV;
  W->traits.location=locW;
  C->traits.location=locC;
}

//! n x m serial dense matrix times k x m conj. transposed serial dense matrix gives m x k sdMat,
//! C=alpha*V*W + beta*C (serial XGEMM wrapper)
extern "C" void SUBR(sdMat_times_sdMatT)(_ST_ alpha, TYPE(const_sdMat_ptr) vV,
                                         TYPE(const_sdMat_ptr) vW,
                              _ST_ beta, TYPE(sdMat_ptr) vC,
                                         int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  TMP_SET_DENSEMAT_LOCATION(vV,V,locV);
  TMP_SET_DENSEMAT_LOCATION(vW,W,locW);
  TMP_SET_DENSEMAT_LOCATION(vC,C,locC);
  *iflag=0;
  PHIST_PERFCHECK_VERIFY_SMALL;
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN_SMALLDETERMINISTIC(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,W,vW,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,C,vC,*iflag);
#ifdef IS_COMPLEX
  char trans[]="C";
#else
  char trans[]="T";
#endif  
  PHIST_CHK_GERR(ghost_gemm(C, V, (char*)"N", W, trans, (void*)&alpha, (void*)&beta, GHOST_GEMM_NO_REDUCE,GHOST_GEMM_DEFAULT),*iflag);
PHIST_TASK_END(iflag);
  V->traits.location=locV;
  W->traits.location=locW;
  C->traits.location=locC;
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
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,R,vR,*iflag);

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
    PHIST_CHK_IERR(SUBR(mvec_from_device)(R,iflag),*iflag);
    *iflag=1-rank;
    return;
  }// case ncols=1: normalize single vector

  PHIST_DEB("mvec_QR: multi-vector case\n");

#if defined(PHIST_HAVE_TEUCHOS)&&defined(PHIST_HAVE_KOKKOS)&&defined(BELOS_HAVE_TSQR)

#ifdef GHOST_HAVE_CUDA
  static int any_cuda=-1;
  if (any_cuda==-1)
  {
    ghost_type ghost_type;
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
  ghost_densemat* Vcopy=NULL, *Qcopy=NULL;
  ghost_densemat_traits vtraits = V->traits;
    vtraits.storage=GHOST_DENSEMAT_COLMAJOR;  
    vtraits.flags = (ghost_densemat_flags)((int)vtraits.flags & ~(int)GHOST_DENSEMAT_VIEW);
    //vtraits.flags = (ghost_densemat_flags)((int)vtraits.flags & ~(int)GHOST_DENSEMAT_PERMUTED); // melven: do we need this?
    vtraits.ncolsorig=vtraits.ncols;
    //vtraits.nrowsorig=vtraits.nrows;
  ghost_densemat_create(&Vcopy,V->map,vtraits);
  ghost_densemat_create(&Qcopy,V->map,vtraits);
      
  {
    PHIST_ENTER_FCN("TSQR_memtranspose");
    // this allocates the memory for the vector, copies and memTransposes the data
    PHIST_CHK_GERR(ghost_densemat_init_densemat(Vcopy,V,0,0),*iflag);
    // this allocates the result vector
    ST zero = st::zero();
    PHIST_CHK_GERR(ghost_densemat_init_val(Qcopy,&zero),*iflag);
  }

  // wrapper class for ghost_densemat for calling Belos.
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
    PHIST_CHK_GERR(ghost_densemat_init_densemat(V,Qcopy,0,0),*iflag);
  }
  ghost_densemat_destroy(Vcopy);
  ghost_densemat_destroy(Qcopy);
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
extern "C" void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, phist_Dmvec* reV, phist_Dmvec* imV, int *iflag)
# else
extern "C" void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, phist_Smvec* reV, phist_Smvec* imV, int *iflag)
#endif
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=0;
  int _jmin=0,_jmax;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&_jmax,iflag),*iflag);
  _jmax--;
  PHIST_PERFCHECK_VERIFY_MVEC_SET_BLOCK(V,V,_jmin,_jmax,iflag);
  TEST_MVEC_MAPS_SAME(V,reV,iflag)
  TEST_MVEC_MAPS_SAME(V,imV,iflag)
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,src,V,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,re,reV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,im,imV,*iflag);
  
  PHIST_CHK_GERR(ghost_densemat_init_complex(re,im,src),*iflag);
  PHIST_TASK_END(iflag);
}

# ifdef IS_DOUBLE
extern "C" void SUBR(mvec_combine)(TYPE(mvec_ptr) V, phist_Dconst_mvec_ptr reV, phist_Dconst_mvec_ptr imV, int *iflag)
# else
extern "C" void SUBR(mvec_combine)(TYPE(mvec_ptr) V, phist_Sconst_mvec_ptr reV, phist_Sconst_mvec_ptr imV, int *iflag)
# endif
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *iflag=0;
  int _jmin=0;
  int _jmax;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&_jmax,iflag),*iflag);
  _jmax--;
  PHIST_PERFCHECK_VERIFY_MVEC_SET_BLOCK(V,V,_jmin,_jmax,iflag);
  TEST_MVEC_MAPS_SAME(V,reV,iflag)
  TEST_MVEC_MAPS_SAME(V,imV,iflag)
PHIST_TASK_DECLARE(ComputeTask)
PHIST_TASK_BEGIN(ComputeTask)
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,vec,V,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,re,reV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,im,imV,*iflag);
  
  PHIST_CHK_GERR(ghost_densemat_init_real(vec,re,im),*iflag);
  PHIST_TASK_END(iflag);
}
#endif

//! create a sparse matrix from a row func and use a distribution prescribed by a given map
extern "C" void SUBR(sparseMat_create_fromRowFuncAndContext)(TYPE(sparseMat_ptr) *vA,
        phist_const_context_ptr vctx,
        phist_lidx maxnne,phist_sparseMat_rowFunc rowFunPtr,void* last_arg,
        int *iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  int iflag_in=*iflag;
  int outlev = *iflag&PHIST_SPARSEMAT_QUIET ? PHIST_DEBUG : PHIST_INFO;
  int own_map= *iflag&PHIST_SPARSEMAT_OWN_MAPS;

  PHIST_CAST_PTR_FROM_VOID(const ghost_context,ctx,vctx,*iflag);

  *iflag=0;
  PHIST_TASK_DECLARE(ComputeTask)
  PHIST_TASK_BEGIN(ComputeTask)

  // TODO: since we're getting a map, it may be that there is already another 
  // sparseMat that has this map. We *actually* want to share the context and 
  // traits (e.g. C and sigma) of that sparseMat, not it's row_map or col_map.
  // So this is only a hotfix, in the long term we should have a function     
  // sparseMat_clone_shape or something alike (note mvec_clone_shape, which I 
  // recently added).
  int sellC=PHIST_SELL_C,sellSigma=PHIST_SELL_SIGMA;
  int sellC_suggested, sellSigma_suggested;
  phist::ghost_internal::get_C_sigma(&sellC_suggested,&sellSigma_suggested,iflag_in,map->mpicomm);
  if (sellC<0) sellC=sellC_suggested;
  if (sellSigma<0) sellSigma=sellSigma_suggested;
  
  PHIST_SOUT(outlev, "Creating sparseMat with SELL-%d-%d format.\n", sellC, sellSigma);

  ghost_sparsemat_src_rowfunc src = GHOST_SPARSEMAT_SRC_ROWFUNC_INITIALIZER;
  src.func = rowFunPtr;
  src.maxrowlen = maxnne;
  src.arg=last_arg;
  src.gnrows = map->gdim;
  src.gncols = map->gdim;


/*
  // TODO: introduce ghost functions context_comm_initialized() and context_clone() etc.
  if (map->ctx->wishlist!=NULL)
  {
    // Clone the context without the communication data structures
    //       (but with the complete permutation info)
    ghost_lidx nrows=map->ctx->row_map->dim;
    PHIST_CHK_IERR(context_create(&ctx,map->ctx->row_map->gdim,map->ctx->col_map->gdim,
        map->ctx->flags,&src,GHOST_SPARSEMAT_SRC_FUNC,map->ctx->mpicomm,map->ctx->weight,iflag),*iflag);
    // share permutation info with the original context
    ctx->col_map->loc_perm=map->ctx->col_map->loc_perm;
    ctx->col_map->loc_perm_inv=map->ctx->col_map->loc_perm_inv;
    ctx->col_map->glb_perm=map->ctx->col_map->glb_perm;
    ctx->col_map->glb_perm_inv=map->ctx->col_map->glb_perm_inv;
  }
*/
  ghost_sparsemat* mat = NULL;

  ghost_sparsemat_traits mtraits=(ghost_sparsemat_traits)GHOST_SPARSEMAT_TRAITS_INITIALIZER;

  mtraits.C=sellC;
  mtraits.sortScope=sellSigma;

#ifdef FIX_ME
  if (!own_map && (mtraits.sortScope>1 || mtraits.flags&GHOST_SPARSEMAT_PERMUTE) )
  {
    // we have to disable sorting for now and reset the flags after the matrix has been created with the
    // given permutations.
    mtraits.flags = (ghost_sparsemat_flags)(mtraits.flags & ~GHOST_SPARSEMAT_PERMUTE);
    mtraits.sortScope=1;
  }
#endif
  PHIST_SOUT(outlev, "Creating sparseMat with SELL-%d-%d format.\n", mtraits.C, mtraits.sortScope);

  mtraits.datatype = st::ghost_dt;
  PHIST_CHK_GERR(ghost_sparsemat_create(&mat,ctx,&mtraits,1),*iflag);                               
  mtraits.flags=GHOST_SPARSEMAT_DEFAULT; // map->mtraits_template.flags;
  double weight=phist::ghost_internal::get_proc_weight();
  //TODO: I want to create the matrix with the given context, at the moment this will re-initialize
  //      the context!!!
  PHIST_CHK_GERR(ghost_sparsemat_init_rowfunc(mat,&src,map->mpicomm,weight),*iflag);
//#if PHIST_OUTLEV >= PHIST_VERBOSE
  char *str;
  ghost_context_string(&str,ctx);
  PHIST_SOUT(outlev,"%s\n",str);
  free(str); str = NULL;
  ghost_sparsemat_info_string(&str,mat);
  PHIST_SOUT(outlev,"%s\n",str);
  free(str); str = NULL;
//#endif
  *vA = (TYPE(sparseMat_ptr))mat;

PHIST_TASK_END(iflag);
  return;
}

extern "C" void SUBR(sparseMat_create_fromRowFunc)(TYPE(sparseMat_ptr) *vA, phist_const_comm_ptr vcomm,
                phist_gidx nrows, phist_gidx ncols, phist_lidx maxnne,
                phist_sparseMat_rowFunc rowFunPtr, void* last_arg, int *iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#ifdef FIX_ME
  int iflag_in=*iflag;  
  int outlev = *iflag&PHIST_SPARSEMAT_QUIET ? PHIST_DEBUG : PHIST_INFO;
  *iflag=0;
  phist_map_ptr map=NULL;
  PHIST_CHK_IERR(phist_map_create(&map,vcomm,nrows,iflag),*iflag);
  

  *iflag=iflag_in | PHIST_SPARSEMAT_OWN_MAPS;
  PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFuncAndMap)(vA,map,maxnne,rowFunPtr,last_arg,iflag),*iflag);
#else
  *iflag=PHIST_NOT_IMPLEMENTED;
#endif
}

extern "C" void SUBR(mvec_write_bin)(TYPE(const_mvec_ptr) vV, const char* filename, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);
  PHIST_CHK_GERR(ghost_densemat_to_file(V,(char*)filename,V->map->mpicomm),*iflag);
}

extern "C" void SUBR(mvec_read_bin)(TYPE(mvec_ptr) vV, const char* filename, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,V,vV,*iflag);
  PHIST_CHK_GERR(ghost_densemat_init_file(V,(char*)filename,V->map->mpicomm),*iflag);
}

extern "C" void SUBR(sdMat_write_bin)(TYPE(const_sdMat_ptr) vM, const char* filename, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,M,vM,*iflag);
  PHIST_CHK_GERR(ghost_densemat_to_file(M,(char*)filename,MPI_COMM_SELF),*iflag);
}

extern "C" void SUBR(sdMat_read_bin)(TYPE(sdMat_ptr) vM, const char* filename, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(ghost_densemat,M,vM,*iflag);
  PHIST_CHK_GERR(ghost_densemat_init_file(M,(char*)filename,MPI_COMM_SELF),*iflag);
}
