/*! \file builtin/kernels_def.hpp
 * included by builtin/kernels.cpp
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 * \author "Jonas Thies <Jonas.Thies@DLR.de>
 *
*/

#ifndef PHIST_HAVE_MPI  
#error "builtin kernels don't work without MPI"
#endif

// Declaration of Fortran implemented functions
extern "C" {
  void SUBR(crsMat_create_fromRowFunc_f)(TYPE(sparseMat_ptr)*,const_comm_ptr_t comm,gidx_t,gidx_t, 
      lidx_t, int (*)(ghost_gidx_t,ghost_lidx_t*,ghost_gidx_t*,void*), int*);
  void SUBR(crsMat_delete_f)(TYPE(sparseMat_ptr) A, int* iflag);
  void SUBR(crsMat_get_map_f)(TYPE(const_sparseMat_ptr),const_map_ptr_t*,int*);
  void SUBR(crsMat_read_mm_f)(void*A,const_comm_ptr_t comm, int fname_len, const char* fname, int* iflag);
  void SUBR(mvecT_times_mvec_f)(_ST_,TYPE(const_mvec_ptr),TYPE(const_mvec_ptr),_ST_,TYPE(sdMat_ptr),int*);
  void SUBR(mvecT_times_mvec_times_sdMat_inplace_f)(_ST_,TYPE(const_mvec_ptr),TYPE(const_mvec_ptr),TYPE(const_sdMat_ptr),_ST_,TYPE(sdMat_ptr),int*);
  void SUBR(mvec_QR_f)(TYPE(mvec_ptr),TYPE(sdMat_ptr),int*);
  void SUBR(mvec_add_mvec_f)(_ST_,TYPE(const_mvec_ptr),_ST_,TYPE(mvec_ptr),int*);
  void SUBR(mvec_create_f)(TYPE(mvec_ptr)*,const_map_ptr_t,lidx_t,int*);
  void SUBR(mvec_delete_f)(TYPE(mvec_ptr),int*);
  void SUBR(mvec_dot_mvec_f)(TYPE(const_mvec_ptr),TYPE(const_mvec_ptr),_ST_*,int*);
  void SUBR(mvec_extract_view_f)(TYPE(mvec_ptr),_ST_**,lidx_t*,int*);
  void SUBR(mvec_gather_mvecs_f)(TYPE(mvec_ptr),TYPE(const_mvec_ptr) W[], int, int*);
  void SUBR(mvec_get_block_f)(TYPE(const_mvec_ptr),TYPE(mvec_ptr),int,int,int*);
  void SUBR(mvec_get_map_f)(TYPE(const_mvec_ptr),const_map_ptr_t*,int*);
  void SUBR(mvec_my_length_f)(TYPE(const_mvec_ptr),lidx_t*,int*);
  void SUBR(mvec_norm2_f)(TYPE(const_mvec_ptr),_MT_*,int*);
  void SUBR(mvec_num_vectors_f)(TYPE(const_mvec_ptr),int*,int*);
  void SUBR(mvec_print_f)(TYPE(const_mvec_ptr),int*);
  void SUBR(mvec_put_value_f)(TYPE(mvec_ptr),_ST_,int*);
  void SUBR(mvec_put_func_f)(TYPE(mvec_ptr),int(*)(ghost_gidx_t,ghost_lidx_t,void*),int*);
  void SUBR(mvec_random_f)(TYPE(mvec_ptr),int*);
  void SUBR(mvec_scale_f)(TYPE(mvec_ptr),_ST_,int*);
  void SUBR(mvec_scatter_mvecs_f)(TYPE(const_mvec_ptr),TYPE(mvec_ptr) W[], int, int*);
  void SUBR(mvec_set_block_f)(TYPE(mvec_ptr),TYPE(const_mvec_ptr),int,int,int*);
  void SUBR(mvec_times_sdMat_f)(_ST_,TYPE(const_mvec_ptr),TYPE(const_sdMat_ptr),_ST_,TYPE(mvec_ptr),int*);
  void SUBR(mvec_times_sdMat_augmented_f)(_ST_,TYPE(const_mvec_ptr),TYPE(const_sdMat_ptr),_ST_,TYPE(mvec_ptr),TYPE(sdMat_ptr),int*);
  void SUBR(mvec_times_sdMat_add_mvec_times_sdMat_f)(TYPE(const_mvec_ptr),TYPE(const_sdMat_ptr),TYPE(mvec_ptr),TYPE(const_sdMat_ptr),int*);
  void SUBR(mvec_times_sdMat_inplace_f)(TYPE(mvec_ptr) V, TYPE(const_sdMat_ptr) M, int*);
  void SUBR(mvec_to_mvec_f)(TYPE(const_mvec_ptr),TYPE(mvec_ptr),int*);
  void SUBR(mvec_vadd_mvec_f)(const _ST_*,TYPE(const_mvec_ptr),_ST_,TYPE(mvec_ptr),int*);
  void SUBR(mvec_view_block_f)(TYPE(mvec_ptr),TYPE(mvec_ptr)*,int,int,int*);
  void SUBR(mvec_vscale_f)(TYPE(mvec_ptr),const _ST_*,int*);
  void SUBR(sdMatT_times_sdMat_f)(_ST_,TYPE(const_sdMat_ptr),TYPE(const_sdMat_ptr),_ST_,TYPE(sdMat_ptr),int*);
  void SUBR(sdMat_times_sdMatT_f)(_ST_,TYPE(const_sdMat_ptr),TYPE(const_sdMat_ptr),_ST_,TYPE(sdMat_ptr),int*);
  void SUBR(sdMat_add_sdMat_f)(_ST_,TYPE(const_sdMat_ptr),_ST_,TYPE(sdMat_ptr),int*);
  void SUBR(sdMatT_add_sdMat_f)(_ST_,TYPE(const_sdMat_ptr),_ST_,TYPE(sdMat_ptr),int*);
  void SUBR(sdMat_create_f)(TYPE(sdMat_ptr)*,int,int,const_comm_ptr_t,int*);
  void SUBR(sdMat_create_view_f)(TYPE(sdMat_ptr)*,const_comm_ptr_t,void*, lidx_t, int,int,int*);
  void SUBR(sdMat_delete_f)(TYPE(sdMat_ptr),int*);
  void SUBR(sdMat_extract_view_f)(TYPE(sdMat_ptr),_ST_**,lidx_t*,int*);
#ifdef PHIST_HIGH_PRECISION_KERNELS
  void SUBR(sdMat_extract_error_f)(TYPE(sdMat_ptr),_ST_**,int*);
#endif
  void SUBR(sdMat_get_block_f)(TYPE(const_mvec_ptr),TYPE(mvec_ptr),int,int,int,int,int*);
  void SUBR(sdMat_get_ncols_f)(TYPE(const_sdMat_ptr),int*,int*);
  void SUBR(sdMat_get_nrows_f)(TYPE(const_sdMat_ptr),int*,int*);
  void SUBR(sdMat_print_f)(TYPE(const_sdMat_ptr),int*);
  void SUBR(sdMat_put_value_f)(TYPE(sdMat_ptr),_ST_,int*);
  void SUBR(sdMat_random_f)(TYPE(sdMat_ptr),int*);
  void SUBR(sdMat_identity_f)(TYPE(sdMat_ptr),int*);
  void SUBR(sdMat_set_block_f)(TYPE(mvec_ptr),TYPE(const_mvec_ptr),int,int,int,int,int*);
  void SUBR(sdMat_times_sdMat_f)(_ST_,TYPE(const_sdMat_ptr),TYPE(const_sdMat_ptr),_ST_,TYPE(sdMat_ptr),int*);
  void SUBR(sdMat_view_block_f)(TYPE(mvec_ptr),TYPE(mvec_ptr)*,int,int,int,int,int*);
}

extern "C" void SUBR(type_avail)(int *iflag)
{
  *iflag=0;
}

// NOTE: We accept the flags PHIST_SPARSEMAT_REPART (if ParMETIS is available,
//       the matrix will be repartitioned), and SPARSEMAT_OPT_CARP (compute a 
//       local dist-2 coloring if ColPack is available and omp_get_num_threads
//       is larger than 1). The flag SPARSEMAT_DIST2_COLOR can be used for    
//       benchmarking the performance impact of local permutation according to
//       the coloring, or for forcing D2-coloring even with one thread. The   
//       behavior is determined by the two flags like this:
//       
//                 DIST2_COLOR        0                 1
//      OPT_CARP
//          0                   no coloring         local permutation + WARNING
//                                                  (breaks MPI and CARP, only for
//                                                  benchmarking impact on spMVM)
//                                                  
//          1              coloring without local   like MC-CARP_CG but also use
//                         permutation (should be   coloring kernel if only one 
//                         used for MC-CARP_CG),    OpenMP thread is available.
//                         no coloring if only 
//                         one OpenMP thread
//      
//      note that specifying OPT_CARP has no performance impact on the spMVM at all.
extern "C" void SUBR(sparseMat_read_mm)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t vcomm,
        const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  if (filename==NULL)
  {
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
  if (*iflag&PHIST_SPARSEMAT_OPT_CARP)
  {
    if ((*iflag&PHIST_SPARSEMAT_DIST2_COLOR)==0)
    {
      int nthreads=1;
#pragma omp parallel
#pragma omp master
      nthreads=omp_get_num_threads();
      if (nthreads==1) 
      {
        // disable coloring for optimal performance.
        PHIST_SOUT(PHIST_INFO,"NOTE: You indicated that the matrix will be used in CARP-CG,\n"
                            "      as it seems that there is only one OpenMP thread, I \n"
                            "      will not construct a coloring, so subsequent CARP sweeps will be sequential per MPI process.\n"
                            "      If you want to use the coloring kernel anyway, specify the flag \n"
                            "       PHIST_SPARSEMAT_OPT_CARP|PHIST_SPARSEMAT_DIST2_COLOR to enforce it.\n"
                        );
        *iflag&=~PHIST_SPARSEMAT_OPT_CARP;
      }
    }
    else
    {
      // force coloring even for one thread, but otherwise like
      // PHIST_OPT_CARP
      *iflag&=~PHIST_SPARSEMAT_DIST2_COLOR;
    }
  }

  PHIST_CHK_IERR(SUBR(crsMat_read_mm_f)(A,vcomm,strlen(filename),filename,iflag),*iflag);
}

extern "C" void SUBR(sparseMat_read_bin)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t vcomm,
const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_hb)(TYPE(sparseMat_ptr)* A, const_comm_ptr_t vcomm,
const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_get_row_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR( SUBR(crsMat_get_map_f) (A,map,iflag), *iflag);
}

extern "C" void SUBR(sparseMat_get_col_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR( SUBR(sparseMat_get_row_map) (A,map,iflag), *iflag);
}

extern "C" void SUBR(sparseMat_get_domain_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR( SUBR(sparseMat_get_row_map) (A,map,iflag), *iflag);
}

extern "C" void SUBR(sparseMat_get_range_map)(TYPE(const_sparseMat_ptr) A, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR( SUBR(sparseMat_get_row_map) (A,map,iflag), *iflag);
}

extern "C" void SUBR(mvec_create)(TYPE(mvec_ptr)* V, 
    const_map_ptr_t map, lidx_t nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_CREATE(map,nvec,iflag);
  PHIST_CHK_IERR( SUBR(mvec_create_f) (V,map,nvec,iflag), *iflag);
}

extern "C" void SUBR(mvec_create_view)(TYPE(mvec_ptr)* V, const_map_ptr_t map, 
    _ST_* values, lidx_t lda, int nvec,
    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_create)(TYPE(sdMat_ptr)* M, 
    int nrows, int ncols, const_comm_ptr_t comm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR(SUBR(sdMat_create_f)(M,nrows,ncols,comm,iflag),*iflag);
}

extern "C" void SUBR(sdMat_create_view)(TYPE(sdMat_ptr)* M, const_comm_ptr_t comm,
        _ST_* values, lidx_t lda, int nrows, int ncols,
        int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_create_view_f)(M,comm,(void*)values,lda,nrows,ncols,iflag),*iflag);
}
                  

extern "C" void SUBR(mvec_my_length)(TYPE(const_mvec_ptr) V, lidx_t* len, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_my_length_f)(V,len,iflag),*iflag);
}

extern "C" void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) V, const_map_ptr_t* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_get_map_f)(V,map,iflag),*iflag);
}

extern "C" void SUBR(mvec_get_comm)(TYPE(const_mvec_ptr) V, const_comm_ptr_t* comm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  const_map_ptr_t map;
  PHIST_CHK_IERR(SUBR(mvec_get_map)(V,&map,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_comm(map,comm,iflag),*iflag);
}

extern "C" void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) V, int* nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors_f)(V,nvec,iflag),*iflag);
}

extern "C" void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) M, int* nrows, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows_f)(M,nrows,iflag),*iflag);
}

extern "C" void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) M, int* ncols, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols_f)(M,ncols,iflag),*iflag);
}

extern "C" void SUBR(mvec_extract_view)(TYPE(mvec_ptr) V, _ST_** val, lidx_t* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_extract_view_f)(V,val,lda,iflag),*iflag);
}

extern "C" void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) V, _ST_** val, lidx_t* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_extract_view_f)(V,val,lda,iflag),*iflag);
}

#ifdef PHIST_HIGH_PRECISION_KERNELS
extern "C" void SUBR(sdMat_extract_error)(TYPE(sdMat_ptr) V, _ST_** err, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_extract_error_f)(V,err,iflag),*iflag);
}
#endif
extern "C" void SUBR(mvec_to_mvec)(TYPE(const_mvec_ptr) v_in, TYPE(mvec_ptr) v_out, int* iflag)
{
#ifndef PHIST_HAVE_PARMETIS
  // we don't need this when no reordering occurs!
  *iflag=PHIST_NOT_IMPLEMENTED;
#else
  PHIST_CHK_IERR(SUBR(mvec_to_mvec_f)(v_in, v_out, iflag),*iflag);
#endif
}

extern "C" void SUBR(mvec_view_block)(TYPE(mvec_ptr) V,
    TYPE(mvec_ptr)* Vblock,
    int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR(SUBR(mvec_view_block_f)(V,Vblock,jmin,jmax,iflag),*iflag);
}

extern "C" void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) V,
    TYPE(mvec_ptr) Vblock,
    int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_GET_BLOCK(V,Vblock,jmin,jmax,iflag);
  PHIST_CHK_IERR(SUBR(mvec_get_block_f)(V,Vblock,jmin,jmax,iflag),*iflag);
}

extern "C" void SUBR(mvec_set_block)(TYPE(mvec_ptr) V,
    TYPE(const_mvec_ptr) Vblock,
    int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_SET_BLOCK(V,Vblock,jmin,jmax,iflag);
  PHIST_CHK_IERR(SUBR(mvec_set_block_f)(V,Vblock,jmin,jmax,iflag),*iflag);
}

extern "C" void SUBR(sdMat_view_block)(TYPE(mvec_ptr) M, 
    TYPE(mvec_ptr)* Mblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR(SUBR(sdMat_view_block_f)(M,Mblock,imin,imax,jmin,jmax,iflag),*iflag);
}

extern "C" void SUBR(sdMat_get_block)(TYPE(const_mvec_ptr) M, 
    TYPE(mvec_ptr) Mblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR(SUBR(sdMat_get_block_f)(M,Mblock,imin,imax,jmin,jmax,iflag),*iflag);
}

extern "C" void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) M, 
    TYPE(const_sdMat_ptr) Mblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR(SUBR(sdMat_set_block_f)(M,Mblock,imin,imax,jmin,jmax,iflag),*iflag);
}

extern "C" void SUBR(sparseMat_delete)(TYPE(sparseMat_ptr) A, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR(SUBR(crsMat_delete_f)(A,iflag), *iflag);
}

extern "C" void SUBR(mvec_delete)(TYPE(mvec_ptr) V, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR( SUBR(mvec_delete_f) (V,iflag),*iflag);
}

extern "C" void SUBR(sdMat_delete)(TYPE(sdMat_ptr) M, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR(SUBR(sdMat_delete_f)(M,iflag),*iflag);
}

extern "C" void SUBR(mvec_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_PUT_VALUE(V,iflag);
  PHIST_CHK_IERR(SUBR(mvec_put_value_f)(V,value,iflag),*iflag);
}

extern "C" void SUBR(mvec_put_func)(TYPE(mvec_ptr) V,
        int (*funPtr)(ghost_gidx_t,ghost_lidx_t,void*), int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_PUT_VALUE(V,iflag);
  PHIST_CHK_IERR(SUBR(mvec_put_func_f)(V,funPtr,iflag),*iflag);
}

extern "C" void SUBR(sdMat_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR(SUBR(sdMat_put_value_f)(V,value,iflag),*iflag);
}

extern "C" void SUBR(mvec_random)(TYPE(mvec_ptr) V, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_PUT_VALUE(V,iflag);
  PHIST_CHK_IERR(SUBR(mvec_random_f)(V,iflag),*iflag);
}

extern "C" void SUBR(mvec_print)(TYPE(const_mvec_ptr) V, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_print_f)(V,iflag),*iflag);
}

extern "C" void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) M, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(sdMat_print_f)(M,iflag),*iflag);
}

extern "C" void SUBR(sdMat_random)(TYPE(sdMat_ptr) M, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR(SUBR(sdMat_random_f)(M,iflag),*iflag);
}

extern "C" void SUBR(sdMat_identity)(TYPE(sdMat_ptr) M, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR(SUBR(sdMat_identity_f)(M,iflag),*iflag);
}

extern "C" void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) V,
    _MT_* vnrm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_DOT_MVEC(V,V,iflag);
  int iflag0=*iflag;
  SUBR(mvec_norm2_f)(V,vnrm,iflag);
  if (*iflag==PHIST_NOT_IMPLEMENTED)
  {
    PHIST_SOUT(PHIST_WARNING,"WARNING: try to use slow fallback version of %s\n",__FUNCTION__);
    Dmvec_ptr_t vtmp=NULL;
    *iflag=0;
    int nvec;
    const_map_ptr_t map=NULL;
    PHIST_CHK_IERR(SUBR(mvec_get_map)(V,&map,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nvec,iflag),*iflag);
    int i=0, istep=4;
    bool realloc=true;
    while (i<nvec)
    {
      while (i+istep>nvec)
      {
        istep/=2;
        if (!realloc)
        {
          PHIST_CHK_IERR(SUBR(mvec_delete)(vtmp,iflag),*iflag);
        }
        realloc=true;
      }
      if (realloc)
      {
        PHIST_CHK_IERR(SUBR(mvec_create)(&vtmp,map,istep,iflag),*iflag);
        realloc=false;
      }
      PHIST_CHK_IERR(SUBR(mvec_get_block)(V,vtmp,i,i+istep-1,iflag),*iflag);
      *iflag=iflag0;
      PHIST_CHK_IERR(SUBR(mvec_norm2_f)(vtmp,vnrm+i,iflag),*iflag);
      i+=istep;
    }//while
    if (!realloc)
    {
      PHIST_CHK_IERR(SUBR(mvec_delete)(vtmp,iflag),*iflag);
    }
  }
}

extern "C" void SUBR(mvec_normalize)(TYPE(mvec_ptr) V,
    _MT_* vnrm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"  
  PHIST_CHK_IERR(SUBR(mvec_norm2_f)(V,vnrm,iflag),*iflag);
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nvec,iflag),*iflag);
  _ST_*scale = new _ST_[nvec];
  for(int i = 0; i < nvec; i++)
    scale[i] = st::one()/vnrm[i];
  PHIST_CHK_IERR(SUBR(mvec_vscale_f)(V,scale,iflag),*iflag);
  delete[] scale;
}

extern "C" void SUBR(mvec_scale)(TYPE(mvec_ptr) V, 
    _ST_ scalar, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_SCALE(V,iflag);
  PHIST_CHK_IERR(SUBR(mvec_scale_f)(V,scalar,iflag),*iflag);
}

extern "C" void SUBR(mvec_vscale)(TYPE(mvec_ptr) V, 
    const _ST_* scalar, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_SCALE(V,iflag);
  PHIST_CHK_IERR(SUBR(mvec_vscale_f)(V,scalar,iflag),*iflag);
}

extern "C" void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) X,
    _ST_ beta,  TYPE(mvec_ptr)       Y, 
    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_ADD_MVEC(alpha,X,beta,Y,iflag);
  PHIST_CHK_IERR(SUBR(mvec_add_mvec_f)(alpha,X,beta,Y,iflag),*iflag);
}

extern "C" void SUBR(mvec_vadd_mvec)(const _ST_ alpha[], TYPE(const_mvec_ptr) X,
    _ST_ beta,  TYPE(mvec_ptr)       Y, 
    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_VADD_MVEC(alpha,X,beta,Y,iflag);
  PHIST_CHK_IERR(SUBR(mvec_vadd_mvec_f)(alpha,X,beta,Y,iflag),*iflag);
}

extern "C" void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) A,
    _ST_ beta,  TYPE(sdMat_ptr)       B, 
    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR(SUBR(sdMat_add_sdMat_f)(alpha,A,beta,B,iflag),*iflag);
}

extern "C" void SUBR(sdMatT_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) A,
    _ST_ beta,  TYPE(sdMat_ptr)       B, 
    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR(SUBR(sdMatT_add_sdMat_f)(alpha,A,beta,B,iflag),*iflag);
}

extern "C" void SUBR(sparseMat_times_mvec_communicate)(TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag = 0;
}

extern "C" void SUBR(sparseMat_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
    TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

#ifdef PHIST_TIMEMONITOR
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(x, &nvec, iflag), *iflag);
  for(int i = 0; i < nvec; i++)
    phist_totalMatVecCount();
#endif

  void SUBR(crsMat_times_mvec_f)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
      TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag);
  PHIST_CHK_IERR(SUBR(crsMat_times_mvec_f)(alpha,A,x,beta,y,iflag),*iflag);
}

extern "C" void SUBR(sparseMatT_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
    TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

//! y[i]=alpha*(A*x[i]+shifts[i]*x[i]) + beta*y[i]
extern "C" void SUBR(sparseMat_times_mvec_vadd_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A,
        const _ST_ shifts[], TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

#ifdef PHIST_TIMEMONITOR
  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(x, &nvec, iflag), *iflag);
  for(int i = 0; i < nvec; i++)
    phist_totalMatVecCount();
#endif

  void SUBR(crsMat_times_mvec_vadd_mvec_f)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
      const _ST_ shifts[], TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag);
  PHIST_CHK_IERR(SUBR(crsMat_times_mvec_vadd_mvec_f)(alpha,A,shifts,x,beta,y,iflag),*iflag);
}

extern "C" void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) v, 
    TYPE(const_mvec_ptr) w, 
    _ST_* s, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_DOT_MVEC(v,w,iflag);
  int iflag0=*iflag;
  SUBR(mvec_dot_mvec_f)(v,w,s,iflag);
  if (*iflag==PHIST_NOT_IMPLEMENTED)
  {
    PHIST_SOUT(PHIST_WARNING,"WARNING: try to use slow fallback version of %s\n",__FUNCTION__);
    Dmvec_ptr_t vtmp=NULL,wtmp=NULL;
    *iflag=0;
    int nvec;
    const_map_ptr_t map=NULL;
    PHIST_CHK_IERR(SUBR(mvec_get_map)(v,&map,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(v,&nvec,iflag),*iflag);
    int i=0, istep=4;
    bool realloc=true;
    while (i<nvec)
    {
      while (i+istep>nvec)
      {
        istep/=2;
        if (!realloc)
        {
          PHIST_CHK_IERR(SUBR(mvec_delete)(vtmp,iflag),*iflag);
          PHIST_CHK_IERR(SUBR(mvec_delete)(wtmp,iflag),*iflag);
        }
        realloc=true;
      }
      if (realloc)
      {
        PHIST_CHK_IERR(SUBR(mvec_create)(&vtmp,map,istep,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(mvec_create)(&wtmp,map,istep,iflag),*iflag);
        realloc=false;
      }
      PHIST_CHK_IERR(SUBR(mvec_get_block)(v,vtmp,i,i+istep-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_get_block)(w,wtmp,i,i+istep-1,iflag),*iflag);
      *iflag=iflag0;
      PHIST_CHK_IERR(SUBR(mvec_dot_mvec_f)(vtmp,wtmp,s+i,iflag),*iflag);
      i+=istep;
    }//while
    if (!realloc)
    {
      PHIST_CHK_IERR(SUBR(mvec_delete)(vtmp,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_delete)(wtmp,iflag),*iflag);
    }
  }
}

extern "C" void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
    TYPE(const_sdMat_ptr) C, 
    _ST_ beta, TYPE(mvec_ptr) W, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT(alpha,V,beta,W,iflag);
  int iflag0=*iflag;
  SUBR(mvec_times_sdMat_f)(alpha,V,C,beta,W,iflag);
  if (*iflag==PHIST_NOT_IMPLEMENTED)
  {
    PHIST_SOUT(PHIST_WARNING,"WARNING: try to use slow fallback version of %s\n",__FUNCTION__);
    const_comm_ptr_t comm=NULL;
    const_map_ptr_t map=NULL;
    Dmvec_ptr_t wtmp=NULL;
    DsdMat_ptr_t ctmp=NULL;
    *iflag=0;
    int nvecv,nvecw;
    PHIST_CHK_IERR(SUBR(mvec_get_map)(W,&map,iflag),*iflag);
    PHIST_CHK_IERR(phist_map_get_comm(map,&comm,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nvecv,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&nvecw,iflag),*iflag);
    int i=0, istep=4;
    bool realloc=true;
    while (i<nvecw)
    {
      while (i+istep>nvecw)
      {
        istep/=2;
        if (!realloc)
        {
          PHIST_CHK_IERR(SUBR(mvec_delete)(wtmp,iflag),*iflag);
          PHIST_CHK_IERR(SUBR(sdMat_delete)(ctmp,iflag),*iflag);
        }
        realloc=true;
      }
      if (realloc)
      {
        PHIST_CHK_IERR(SUBR(mvec_create)(&wtmp,map,istep,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_create)(&ctmp,nvecv,istep,comm,iflag),*iflag);
        realloc=false;
      }
      PHIST_CHK_IERR(SUBR(mvec_get_block)(W,wtmp,i,i+istep-1,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_get_block)(C,ctmp,0,nvecv-1,i,i+istep-1,iflag),*iflag);
      *iflag=iflag0;
      PHIST_CHK_IERR(SUBR(mvec_times_sdMat_f)(alpha,V,ctmp,beta,wtmp,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_set_block)(W,wtmp,i,i+istep-1,iflag),*iflag);
      i+=istep;
    }//while
    if (!realloc)
    {
      PHIST_CHK_IERR(SUBR(mvec_delete)(wtmp,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_delete)(ctmp,iflag),*iflag);
    }
  }
}

extern "C" void SUBR(mvec_times_sdMat_augmented)(_ST_ alpha,  TYPE(const_mvec_ptr)  V, 
                                                              TYPE(const_sdMat_ptr) C, 
                                                 _ST_ beta,   TYPE(mvec_ptr)        W,
                                                              TYPE(sdMat_ptr)       D,
                                                 int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT(alpha,V,beta,W,iflag);
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat_augmented_f)(alpha,V,C,beta,W,D,iflag),*iflag);
}

extern "C" void SUBR(mvec_times_sdMat_add_mvec_times_sdMat)(TYPE(const_mvec_ptr) V, 
                                                            TYPE(const_sdMat_ptr) C,
                                                            TYPE(mvec_ptr) W, 
                                                            TYPE(const_sdMat_ptr) D,
                                                            int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT(st::one(),V,st::one(),W,iflag);
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat_add_mvec_times_sdMat_f)(V,C,W,D,iflag),*iflag);
}

extern "C" void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
    TYPE(const_sdMat_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR(SUBR(sdMat_times_sdMat_f)(alpha,V,W,beta,C,iflag),*iflag);
}

extern "C" void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
    TYPE(const_sdMat_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR(SUBR(sdMatT_times_sdMat_f)(alpha,V,W,beta,C,iflag),*iflag);
}

extern "C" void SUBR(sdMat_times_sdMatT)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
    TYPE(const_sdMat_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR(SUBR(sdMat_times_sdMatT_f)(alpha,V,W,beta,C,iflag),*iflag);
}

extern "C" void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
    TYPE(const_mvec_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVECT_TIMES_MVEC(V,W,iflag);
  int iflag0=*iflag;
  SUBR(mvecT_times_mvec_f)(alpha,V,W,beta,C,iflag);
  if (*iflag==PHIST_NOT_IMPLEMENTED)
  {
    PHIST_SOUT(PHIST_WARNING,"WARNING: try to use slow fallback version of %s\n",__FUNCTION__);
    const_comm_ptr_t comm=NULL;
    const_map_ptr_t map=NULL;
    Dmvec_ptr_t vtmp=NULL;
    DsdMat_ptr_t ctmp=NULL;
    *iflag=0;
    int nvecv,nvecw;
    PHIST_CHK_IERR(SUBR(mvec_get_map)(W,&map,iflag),*iflag);
    PHIST_CHK_IERR(phist_map_get_comm(map,&comm,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nvecv,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&nvecw,iflag),*iflag);
    int i=0, istep=4;
    bool realloc=true;
    while (i<nvecv)
    {
      while (i+istep>nvecv)
      {
        istep/=2;
        if (!realloc)
        {
          PHIST_CHK_IERR(SUBR(mvec_delete)(vtmp,iflag),*iflag);
          PHIST_CHK_IERR(SUBR(sdMat_delete)(ctmp,iflag),*iflag);
        }
        realloc=true;
      }
      if (realloc)
      {
        PHIST_CHK_IERR(SUBR(mvec_create)(&vtmp,map,istep,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_create)(&ctmp,istep,nvecw,comm,iflag),*iflag);
        realloc=false;
      }
      PHIST_CHK_IERR(SUBR(mvec_get_block)(V,vtmp,i,i+istep-1,iflag),*iflag);
      PHIST_SOUT(PHIST_DEBUG,"compute C(%d:%d,%d:%d)\n",i,i+istep-1,0,nvecw-1);
      PHIST_CHK_IERR(SUBR(sdMat_get_block)(C,ctmp,i,i+istep-1,0,nvecw-1,iflag),*iflag);
      *iflag=iflag0;
      PHIST_CHK_IERR(SUBR(mvecT_times_mvec_f)(alpha,vtmp,W,beta,ctmp,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_set_block)(C,ctmp,i,i+istep-1,0,nvecw-1,iflag),*iflag);
      i+=istep;
    }//while
    if (!realloc)
    {
      PHIST_CHK_IERR(SUBR(mvec_delete)(vtmp,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_delete)(ctmp,iflag),*iflag);
    }
  }
}


extern "C" void SUBR(mvecT_times_mvec_times_sdMat_inplace)(_ST_ alpha, TYPE(const_mvec_ptr)  V,
                                                                       TYPE(mvec_ptr)        W,
                                                                       TYPE(const_sdMat_ptr) C,
                                                           _ST_ beta,  TYPE(sdMat_ptr)       D,
                                                           int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVECT_TIMES_MVEC(V,W,iflag);
  PHIST_CHK_IERR(SUBR(mvecT_times_mvec_times_sdMat_inplace_f)(alpha,V,W,C,beta,D,iflag),*iflag);
}


extern "C" void SUBR(mvec_QR)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  // trilinos tsqr would be better but seems to expect column-major mvecs
#ifdef PHIST_HIGH_PRECISION_KERNELS
  int robust=(*iflag&PHIST_ROBUST_REDUCTIONS);
  if (robust)
  {
    // use Cholesky-QR
    int m,rank;
    PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,iflag),*iflag);
    *iflag=PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(st::one(),V,V,st::zero(),R,iflag),*iflag);
    int perm[m];
    PHIST_CHK_IERR(SUBR(sdMat_cholesky)(R,perm,&rank,iflag),*iflag);
    
    // construct inv(R)
    TYPE(sdMat_ptr) R_1=NULL;
    const_comm_ptr_t comm;
    PHIST_CHK_IERR(SUBR(mvec_get_comm)(V,&comm,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_create)(&R_1,m,m,comm,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_identity)(R_1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_backwardSubst_sdMat)(R,perm,rank,R_1,iflag),*iflag);
    *iflag=PHIST_ROBUST_REDUCTIONS;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(V,R_1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(sdMat_delete)(R_1,iflag),*iflag);
    *iflag=m-rank;
    return;
  }
#endif
  PHIST_CHK_NEG_IERR(SUBR(mvec_QR_f)(V,R,iflag),*iflag);
}

extern "C" void SUBR(mvec_gather_mvecs)(TYPE(mvec_ptr) V, TYPE(const_mvec_ptr) W[], int nblocks, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_MARK_AS_EXPERIMENTAL(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_gather_mvecs_f)(V,W,nblocks,iflag),*iflag);
}

extern "C" void SUBR(mvec_scatter_mvecs)(TYPE(const_mvec_ptr) V, TYPE(mvec_ptr) W[], int nblocks, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_MARK_AS_EXPERIMENTAL(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_scatter_mvecs_f)(V,W,nblocks,iflag),*iflag);
}

extern "C" void SUBR(mvec_times_sdMat_inplace)(TYPE(mvec_ptr) V, TYPE(const_sdMat_ptr) M, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT_INPLACE(V,M,iflag);
  int iflag0=*iflag;
  SUBR(mvec_times_sdMat_inplace_f)(V, M, iflag);
  if (*iflag==PHIST_NOT_IMPLEMENTED)
  {
    PHIST_SOUT(PHIST_WARNING,"kernel not found, try to use out-of-place variant instead (performance hazard!)\n");
    TYPE(mvec_ptr) Vtmp=NULL;
    const_map_ptr_t map=NULL;
    PHIST_CHK_IERR(SUBR(mvec_get_map)(V,&map,iflag),*iflag);
    int nvecw;
    PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(M,&nvecw,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_create)(&Vtmp,map,nvecw,iflag),*iflag);
    *iflag=iflag0;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(1.0,V, M, 0.0, Vtmp, iflag), *iflag);    
    PHIST_CHK_IERR(SUBR(mvec_set_block)(V,Vtmp,0,nvecw-1,iflag),*iflag);
    PHIST_CHK_IERR(SUBR(mvec_delete)(Vtmp,iflag),*iflag);
  }
}

// NOTE: see the description of sparseMat_read_mm on how we treat input flags for this function
extern "C" void SUBR(sparseMat_create_fromRowFunc)(TYPE(sparseMat_ptr) *A, const_comm_ptr_t vcomm,
        gidx_t nrows, gidx_t ncols, lidx_t maxnne, 
                int (*rowFunPtr)(ghost_gidx_t,ghost_lidx_t*,ghost_gidx_t*,void*), 
                int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
/*
std::cout << "iflag="<<*iflag<<std::endl;
std::cout << "iflag&OPT_CARP="<<(*iflag&PHIST_SPARSEMAT_OPT_CARP)<<std::endl;
std::cout << "iflag&DIST2_COLOR="<<(*iflag&PHIST_SPARSEMAT_DIST2_COLOR)<<std::endl;
*/
  if (*iflag&PHIST_SPARSEMAT_OPT_CARP)
  {
    if ((*iflag&PHIST_SPARSEMAT_DIST2_COLOR)==0)
    {
      int nthreads=1;
#pragma omp parallel
#pragma omp master
      nthreads=omp_get_num_threads();
      std::cout << "nthreads="<<nthreads<<std::endl;
      if (nthreads==1)
      {
        // disable coloring for optimal performance.
        PHIST_SOUT(PHIST_INFO,"NOTE: You indicated that the matrix will be used in CARP-CG,\n"
                            "      as it seems that there is only one OpenMP thread, I \n"
                            "      will not construct a coloring, so subsequent CARP sweeps will be sequential per MPI process.\n"
                            "      If you want to use the coloring kernel anyway, specify the flag \n"
                            "       PHIST_SPARSEMAT_OPT_CARP|PHIST_SPARSEMAT_DIST2_COLOR to enforce it.\n"
                        );
        *iflag&=~PHIST_SPARSEMAT_OPT_CARP;
      }
    }
    else
    {
      // force coloring even for one thread, but otherwise like
      // PHIST_OPT_CARP
      *iflag&=~PHIST_SPARSEMAT_DIST2_COLOR;
    }
  }

  PHIST_CHK_IERR(SUBR(crsMat_create_fromRowFunc_f)(A, vcomm, nrows, ncols, maxnne, 
        rowFunPtr, iflag), *iflag);
}

#include "../kernels_nogpu.c"


