/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
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
  void SUBR(crsMat_create_fromRowFunc_f)(TYPE(sparseMat_ptr)*,phist_const_comm_ptr,phist_gidx,phist_gidx, 
      phist_lidx, phist_sparseMat_rowFunc, phist_sparseMat_rowFuncConstructor, void*, int*);
// note: this function doesn't exist yet!
//  void SUBR(crsMat_create_fromRowFuncAndContext_f)(TYPE(sparseMat_ptr)*,phist_const_comm_ptr,phist_const_context_ptr, 
//      phist_lidx, phist_sparseMat_rowFunc,void*, int*);
  void SUBR(crsMat_delete_f)(TYPE(sparseMat_ptr) A, int* iflag);
  void SUBR(crsMat_local_nnz_f)(TYPE(const_sparseMat_ptr),int64_t*,int*);
  void SUBR(crsMat_global_nnz_f)(TYPE(const_sparseMat_ptr),int64_t*,int*);
  void SUBR(crsMat_get_map_f)(TYPE(const_sparseMat_ptr),phist_const_map_ptr*,int*);
  void SUBR(crsMat_read_mm_f)(void*A,phist_const_comm_ptr comm, int fname_len, const char* fname, int* iflag);
  void SUBR(crsMat_read_mm_with_map_f)(void*A,phist_const_map_ptr map, int fname_len, const char* fname, int* iflag);
  void SUBR(mvecT_times_mvec_f)(_ST_,TYPE(const_mvec_ptr),TYPE(const_mvec_ptr),_ST_,TYPE(sdMat_ptr),int*);
  void SUBR(mvecT_times_mvec_times_sdMat_inplace_f)(_ST_,TYPE(const_mvec_ptr),TYPE(const_mvec_ptr),TYPE(const_sdMat_ptr),_ST_,TYPE(sdMat_ptr),int*);
  void SUBR(mvec_QR_f)(TYPE(mvec_ptr),TYPE(sdMat_ptr),int*);
  void SUBR(mvec_add_mvec_f)(_ST_,TYPE(const_mvec_ptr),_ST_,TYPE(mvec_ptr),int*);
  void SUBR(mvec_times_mvec_elemwise_f)(_ST_,TYPE(const_mvec_ptr),TYPE(mvec_ptr),int*);
  void SUBR(mvec_create_f)(TYPE(mvec_ptr)*,phist_const_map_ptr,phist_lidx,int*);
  void SUBR(mvec_delete_f)(TYPE(mvec_ptr),int*);
  void SUBR(mvec_dot_mvec_f)(TYPE(const_mvec_ptr),TYPE(const_mvec_ptr),_ST_*,int*);
  void SUBR(mvec_extract_view_f)(TYPE(mvec_ptr),_ST_**,phist_lidx*,int*);
  void SUBR(mvec_gather_mvecs_f)(TYPE(mvec_ptr),TYPE(const_mvec_ptr) W[], int, int*);
  void SUBR(mvec_get_block_f)(TYPE(const_mvec_ptr),TYPE(mvec_ptr),int,int,int*);
  void SUBR(mvec_get_map_f)(TYPE(const_mvec_ptr),phist_const_map_ptr*,int*);
//  void SUBR(mvec_my_length_f)(TYPE(const_mvec_ptr),phist_lidx*,int*);
  void SUBR(mvec_norm2_f)(TYPE(const_mvec_ptr),_MT_*,int*);
  void SUBR(mvec_num_vectors_f)(TYPE(const_mvec_ptr),int*,int*);
  void SUBR(mvec_print_f)(TYPE(const_mvec_ptr),int*);
  void SUBR(mvec_put_value_f)(TYPE(mvec_ptr),_ST_,int*);
  void SUBR(mvec_put_func_f)(TYPE(mvec_ptr),phist_mvec_elemFunc,void*,int*);
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
  void SUBR(sdMat_create_f)(TYPE(sdMat_ptr)*,int,int,phist_const_comm_ptr,int*);
  void SUBR(sdMat_delete_f)(TYPE(sdMat_ptr),int*);
  void SUBR(sdMat_extract_view_f)(TYPE(sdMat_ptr),_ST_**,phist_lidx*,int*);
#ifdef PHIST_HIGH_PRECISION_KERNELS
  void SUBR(sdMat_extract_error_f)(TYPE(sdMat_ptr),_ST_**,int*);
#endif
  void SUBR(sdMat_get_block_f)(TYPE(const_sdMat_ptr),TYPE(sdMat_ptr),int,int,int,int,int*);
  void SUBR(sdMat_get_ncols_f)(TYPE(const_sdMat_ptr),int*,int*);
  void SUBR(sdMat_get_nrows_f)(TYPE(const_sdMat_ptr),int*,int*);
  void SUBR(sdMat_print_f)(TYPE(const_sdMat_ptr),int*);
  void SUBR(sdMat_put_value_f)(TYPE(sdMat_ptr),_ST_,int*);
  void SUBR(sdMat_random_f)(TYPE(sdMat_ptr),int*);
  void SUBR(sdMat_identity_f)(TYPE(sdMat_ptr),int*);
  void SUBR(sdMat_set_block_f)(TYPE(sdMat_ptr),TYPE(const_sdMat_ptr),int,int,int,int,int*);
  void SUBR(sdMat_times_sdMat_f)(_ST_,TYPE(const_sdMat_ptr),TYPE(const_sdMat_ptr),_ST_,TYPE(sdMat_ptr),int*);
  void SUBR(sdMat_view_block_f)(TYPE(sdMat_ptr),TYPE(sdMat_ptr)*,int,int,int,int,int*);
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
extern "C" void SUBR(sparseMat_read_mm)(TYPE(sparseMat_ptr)* A, phist_const_comm_ptr vcomm,
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
        if( !(*iflag & PHIST_SPARSEMAT_QUIET) )
        {
          PHIST_SOUT(PHIST_INFO,"NOTE: You indicated that the matrix will be used in CARP-CG,\n"
                              "      as it seems that there is only one OpenMP thread, I \n"
                              "      will not construct a coloring, so subsequent CARP sweeps will be sequential per MPI process.\n"
                              "      If you want to use the coloring kernel anyway, specify the flag \n"
                              "       PHIST_SPARSEMAT_OPT_CARP|PHIST_SPARSEMAT_DIST2_COLOR to enforce it.\n"
                          );
        }
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

//
extern "C" void SUBR(sparseMat_read_mm_with_context)(TYPE(sparseMat_ptr)* A, phist_const_context_ptr vctx,
        const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  if (filename==NULL)
  {
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
  PHIST_CAST_PTR_FROM_VOID(phist::internal::default_context,ctx,vctx,*iflag);
  PHIST_CHK_IERR(SUBR(crsMat_read_mm_with_map_f)(A,ctx->row_map,strlen(filename),filename,iflag),*iflag);
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

extern "C" void SUBR(sparseMat_read_bin_with_context)(TYPE(sparseMat_ptr)* A, phist_const_context_ptr ctx,
        const char* filename,int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_hb_with_context)(TYPE(sparseMat_ptr)* A, phist_const_context_ptr ctx,
        const char* filename,int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_local_nnz)(TYPE(const_sparseMat_ptr) A, int64_t* local_nnz, int* iflag)
{
  PHIST_CHK_IERR( SUBR(crsMat_local_nnz_f) (A,local_nnz,iflag), *iflag);
}

extern "C" void SUBR(sparseMat_global_nnz)(TYPE(const_sparseMat_ptr) A, int64_t* global_nnz, int* iflag)
{
  PHIST_CHK_IERR( SUBR(crsMat_global_nnz_f) (A,global_nnz,iflag), *iflag);
}


extern "C" void SUBR(sparseMat_get_row_map)(TYPE(const_sparseMat_ptr) A, phist_const_map_ptr* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR( SUBR(crsMat_get_map_f) (A,map,iflag), *iflag);
}

extern "C" void SUBR(sparseMat_get_col_map)(TYPE(const_sparseMat_ptr) A, phist_const_map_ptr* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR( SUBR(sparseMat_get_row_map) (A,map,iflag), *iflag);
}

extern "C" void SUBR(sparseMat_get_domain_map)(TYPE(const_sparseMat_ptr) A, phist_const_map_ptr* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR( SUBR(sparseMat_get_row_map) (A,map,iflag), *iflag);
}

extern "C" void SUBR(sparseMat_get_range_map)(TYPE(const_sparseMat_ptr) A, phist_const_map_ptr* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR( SUBR(sparseMat_get_row_map) (A,map,iflag), *iflag);
}

extern "C" void SUBR(mvec_create)(TYPE(mvec_ptr)* V, 
    phist_const_map_ptr map, phist_lidx nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_CREATE(map,nvec,iflag);
  PHIST_CHK_IERR( SUBR(mvec_create_f) (V,map,nvec,iflag), *iflag);
}

extern "C" void SUBR(sdMat_create)(TYPE(sdMat_ptr)* M, 
    int nrows, int ncols, phist_const_comm_ptr comm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR(SUBR(sdMat_create_f)(M,nrows,ncols,comm,iflag),*iflag);
}

extern "C" void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) V, phist_const_map_ptr* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_get_map_f)(V,map,iflag),*iflag);
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

extern "C" void SUBR(mvec_extract_view)(TYPE(mvec_ptr) V, _ST_** val, phist_lidx* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(mvec_extract_view_f)(V,val,lda,iflag),*iflag);
}

extern "C" void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) V, _ST_** val, phist_lidx* lda, int* iflag)
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
  // the only permutation we may apply is by parmetis, so copy the vectors if they are compatible.
  // Otherwise mvec_add_mvec will return an error code (I hope...)
  PHIST_CHK_IERR(SUBR(mvec_add_mvec)(1.0,v_in,0.0,v_out,iflag),*iflag);
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

extern "C" void SUBR(sdMat_view_block)(TYPE(sdMat_ptr) M, 
    TYPE(sdMat_ptr)* Mblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR(SUBR(sdMat_view_block_f)(M,Mblock,imin,imax,jmin,jmax,iflag),*iflag);
}

extern "C" void SUBR(sdMat_get_block)(TYPE(const_sdMat_ptr) M, 
    TYPE(sdMat_ptr) Mblock,
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
  // from common/default_context.h: if anyone obtained the context from this matrix
  //    it was created at that moment and is now deleted:
  phist::internal::delete_default_context(A);
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
        phist_mvec_elemFunc funPtr,void* last_arg, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_PUT_VALUE(V,iflag);
  PHIST_CHK_IERR(SUBR(mvec_put_func_f)(V,funPtr,last_arg,iflag),*iflag);
}

extern "C" void SUBR(sdMat_put_value)(TYPE(sdMat_ptr) V, _ST_ value, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR(SUBR(sdMat_put_value_f)(V,value,iflag),*iflag);
}

#ifndef PHIST_BUILTIN_RNG
// actually this uses  the same RNG from tools/phist_random.h, but if the user deactivates
// the option, the common function won't be compiled (see common/kernels_common_impl_def.hpp).
extern "C" void SUBR(mvec_random)(TYPE(mvec_ptr) V, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_PUT_VALUE(V,iflag);
  PHIST_CHK_IERR(SUBR(mvec_random_f)(V,iflag),*iflag);
}
#endif

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

#ifndef PHIST_BUILTIN_RNG
extern "C" void SUBR(sdMat_random)(TYPE(sdMat_ptr) M, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_PERFCHECK_VERIFY_SMALL;
  PHIST_CHK_IERR(SUBR(sdMat_random_f)(M,iflag),*iflag);
}
#endif

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
  if (*iflag!=PHIST_NOT_IMPLEMENTED)
  {
    PHIST_CHK_IERR(void(),*iflag);
  }
  else // *iflag==PHIST_NOT_IMPLEMENTED
  {
    static bool first_time=true;
    if (first_time) {PHIST_SOUT(PHIST_WARNING,"Warning: try to use slow fallback version of %s\n",__FUNCTION__); first_time=false;}
    phist_Dmvec_ptr vtmp=NULL;
    *iflag=0;
    int nvec;
    phist_const_map_ptr map=NULL;
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
  PHIST_CHK_IERR(SUBR(mvec_norm2)(V,vnrm,iflag),*iflag);
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

extern "C" void SUBR(mvec_times_mvec_elemwise)(_ST_ alpha, TYPE(const_mvec_ptr) X,
                                                  TYPE(mvec_ptr)       Y, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_TIMES_MVEC_ELEMWISE(alpha,X,Y,iflag);
  PHIST_CHK_IERR(SUBR(mvec_times_mvec_elemwise_f)(alpha,X,Y,iflag),*iflag);
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

  PHIST_COUNT_MATVECS(x)

  PHIST_PERFCHECK_VERIFY_SPMV(alpha,A,_ST_(0),x,beta,y,_ST_(0),_ST_(0),_ST_(0),iflag);

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

  PHIST_COUNT_MATVECS(x);

  PHIST_PERFCHECK_VERIFY_SPMV(alpha,A,_ST_(1),x,beta,y,0.0,0.0,0,iflag);

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
    static bool first_time=true;
    if (first_time) {PHIST_SOUT(PHIST_WARNING,"Warning: try to use slow fallback version of %s\n",__FUNCTION__); first_time=false;}
    phist_Dmvec_ptr vtmp=NULL,wtmp=NULL;
    *iflag=0;
    int nvec;
    phist_const_map_ptr map=NULL;
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
  if (*iflag!=PHIST_NOT_IMPLEMENTED)
  {
    PHIST_CHK_IERR(void(),*iflag);
  }
  else // *iflag==PHIST_NOT_IMPLEMENTED
  {
    static bool first_time=true;
    if (first_time) {PHIST_SOUT(PHIST_WARNING,"Warning: try to use slow fallback version of %s\n",__FUNCTION__); first_time=false;}
    phist_const_comm_ptr comm=NULL;
    phist_const_map_ptr map=NULL;
    phist_Dmvec_ptr wtmp=NULL;
    phist_DsdMat_ptr ctmp=NULL;
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
      PHIST_SOUT(PHIST_DEBUG,"compute from C(%d:%d,%d:%d)\n",0,nvecv-1,i,i+istep-1);
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

extern "C" void SUBR(fused_mvsd_mvTmv)(_ST_ alpha,  TYPE(const_mvec_ptr)  V, 
                                                              TYPE(const_sdMat_ptr) C, 
                                                 _ST_ beta,   TYPE(mvec_ptr)        W,
                                                              TYPE(sdMat_ptr)       D,
                                                 int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVEC_TIMES_SDMAT(alpha,V,beta,W,iflag);
#ifndef PHIST_HIGH_PRECISION_KERNELS
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat_augmented_f)(alpha,V,C,beta,W,D,iflag),*iflag);
#else
  int iflag_ = *iflag;
  SUBR(mvec_times_sdMat_augmented_f)(alpha,V,C,beta,W,D,iflag);
  if( *iflag == PHIST_NOT_IMPLEMENTED )
  {
    // call separatee kernels, not that this introduces intermediate roundoff that is avoided by the fused kernel
    static bool first_time=true;
    if (first_time) 
    {
      PHIST_SOUT(PHIST_WARNING,"Warning: using less accurate fallback implementtion of fused "
                                             "high precision kernel (file %s, line %d)\n",__FILE__,__LINE__);
      first_time=false;
    }                                             
    int iflag1 = iflag_;
    SUBR(mvec_times_sdMat)(alpha,V,C,beta,W,&iflag1);
    int iflag2 = iflag_;
    SUBR(mvecT_times_mvec)(st::one(),W,W,st::zero(),D,&iflag2);
    PHIST_CHK_IERR(*iflag = iflag1,*iflag);
    PHIST_CHK_IERR(*iflag = iflag2,*iflag);
  }
  else
  {
    PHIST_CHK_IERR(void(),*iflag);
  }
#endif
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
#ifndef PHIST_HIGH_PRECISION_KERNELS
  PHIST_CHK_IERR(SUBR(mvec_times_sdMat_add_mvec_times_sdMat_f)(V,C,W,D,iflag),*iflag);
#else
  int flags = *iflag;
  SUBR(mvec_times_sdMat_add_mvec_times_sdMat_f)(V,C,W,D,iflag);
  if( *iflag == PHIST_NOT_IMPLEMENTED )
  {
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat_inplace)(W,D,iflag),*iflag);
    *iflag = flags;
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(st::one(),V,C,st::one(),W,iflag),*iflag);
  }
  else
  {
    PHIST_CHK_IERR(void(),*iflag);
  }
#endif
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
  if (*iflag!=PHIST_NOT_IMPLEMENTED)
  {
    PHIST_CHK_IERR(void(),*iflag);
  }
  else // *iflag==PHIST_NOT_IMPLEMENTED
  {
    static bool first_time=true;
    if (first_time) {PHIST_SOUT(PHIST_WARNING,"Warning: try to use slow fallback version of %s\n",__FUNCTION__); first_time=false;}
    phist_const_comm_ptr comm=NULL;
    phist_const_map_ptr map=NULL;
    phist_Dmvec_ptr vtmp=NULL, wtmp=NULL;
    phist_DsdMat_ptr ctmp=NULL;
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
        istep/=2;
      int j=0, jstep=4;
      realloc=true;
      while (j<nvecw)
      {
        while (j+jstep>nvecw)
        {
          jstep/=2;
          realloc=true;
        }
        if (realloc)
        {
          PHIST_CHK_IERR(SUBR(mvec_delete)(vtmp,iflag),*iflag); vtmp = NULL;
          PHIST_CHK_IERR(SUBR(mvec_delete)(wtmp,iflag),*iflag); vtmp = NULL;
          PHIST_CHK_IERR(SUBR(sdMat_delete)(ctmp,iflag),*iflag); ctmp = NULL;

          PHIST_CHK_IERR(SUBR(mvec_create)(&vtmp,map,istep,iflag),*iflag);
          PHIST_CHK_IERR(SUBR(mvec_create)(&wtmp,map,jstep,iflag),*iflag);
          PHIST_CHK_IERR(SUBR(sdMat_create)(&ctmp,istep,jstep,comm,iflag),*iflag);
          realloc=false;
        }
        PHIST_CHK_IERR(SUBR(mvec_get_block)(V,vtmp,i,i+istep-1,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(mvec_get_block)(W,wtmp,j,j+jstep-1,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_get_block)(C,ctmp,i,i+istep-1,j,j+jstep-1,iflag),*iflag);
        *iflag=iflag0;
        PHIST_CHK_IERR(SUBR(mvecT_times_mvec_f)(alpha,vtmp,wtmp,beta,ctmp,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(sdMat_set_block)(C,ctmp,i,i+istep-1,j,j+jstep-1,iflag),*iflag);
        j+=jstep;
      }
      i+=istep;
    }//while
    if (!realloc)
    {
      PHIST_CHK_IERR(SUBR(mvec_delete)(vtmp,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(mvec_delete)(wtmp,iflag),*iflag);
      PHIST_CHK_IERR(SUBR(sdMat_delete)(ctmp,iflag),*iflag);
    }
  }
}


extern "C" void SUBR(fused_mvsdi_mvTmv)(_ST_ alpha, TYPE(const_mvec_ptr)  V,
                                                                       TYPE(mvec_ptr)        W,
                                                                       TYPE(const_sdMat_ptr) C,
                                                           _ST_ beta,  TYPE(sdMat_ptr)       D,
                                                           int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_PERFCHECK_VERIFY_MVECT_TIMES_MVEC_TIMES_SDMAT(V,W,iflag);
#ifndef PHIST_HIGH_PRECISION_KERNELS
  PHIST_CHK_IERR(SUBR(mvecT_times_mvec_times_sdMat_inplace_f)(alpha,V,W,C,beta,D,iflag),*iflag);
#else
  int iflag_ = *iflag;
  SUBR(mvecT_times_mvec_times_sdMat_inplace_f)(alpha,V,W,C,beta,D,iflag);
  if( *iflag == PHIST_NOT_IMPLEMENTED )
  {
    // call separatee kernels, not that this introduces intermediate roundoff that is avoided by the fused kernel
    static bool first_time=false;
    if (first_time) 
    {
      PHIST_SOUT(PHIST_WARNING,"Warning: using less accurate fallback implementtion of fused "
                                             "high precision kernel (file %s, line %d)\n",__FILE__,__LINE__);
      first_time=false;
    }                                             
    int iflag1 = iflag_;
    SUBR(mvec_times_sdMat_inplace)(W,C,&iflag1);
    int iflag2 = iflag_;
    SUBR(mvecT_times_mvec)(alpha,V,W,beta,D,&iflag2);
    PHIST_CHK_IERR(*iflag = iflag1,*iflag);
    PHIST_CHK_IERR(*iflag = iflag2,*iflag);
  }
  else
  {
    PHIST_CHK_IERR(void(),*iflag);
  }
#endif
}


extern "C" void SUBR(mvec_QR)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  // the builtin mvec_QR based on iterated MGS is not too nice from a performance point of view,
  // if orthog finds NOT_IMPLEMENTED it will use CholQR instead
  *iflag=PHIST_NOT_IMPLEMENTED;
//  PHIST_CHK_NEG_IERR(SUBR(mvec_QR_f)(V,R,iflag),*iflag);
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
    phist_const_map_ptr map=NULL;
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

extern "C" void SUBR(sparseMat_create_fromRowFuncAndContext)(TYPE(sparseMat_ptr) *vA, phist_const_context_ptr vctx,
        phist_lidx maxnne,phist_sparseMat_rowFunc rowFunPtr,
        void* last_arg,
        int *iflag)
{
  PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFuncWithConstructorAndContext)
        (vA, vctx, maxnne, rowFunPtr, nullptr, last_arg, iflag), *iflag);
}

extern "C" void SUBR(sparseMat_create_fromRowFuncWithConstructorAndContext)(TYPE(sparseMat_ptr) *vA, phist_const_context_ptr vctx,
        phist_lidx maxnne,phist_sparseMat_rowFunc rowFunPtr,
        phist_sparseMat_rowFuncConstructor rowFunConstructorPtr,
        void* last_arg,
        int *iflag)
{
  //TODO - we don't actually support this up to now. So as a fallback, check if the given map
  //       happens to be the default map anyway, and otherwise return -99 (not implemented)
  phist_gidx N;
  PHIST_CAST_PTR_FROM_VOID(phist::internal::default_context,ctx,vctx,*iflag);
  PHIST_CHK_IERR(*iflag=(ctx->range_map!=NULL &&
                          ctx->row_map!=NULL   &&
                          ctx->domain_map!=NULL)? 0:PHIST_NOT_IMPLEMENTED,*iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(ctx->range_map,&N,iflag),*iflag);
  phist_const_comm_ptr vcomm=NULL;
  PHIST_CHK_IERR(phist_map_get_comm(ctx->range_map,&vcomm,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(crsMat_create_fromRowFunc_f)(vA, vcomm, N, N, maxnne,
        rowFunPtr, rowFunConstructorPtr, last_arg, iflag), *iflag);
  phist_const_map_ptr new_map=NULL;
  PHIST_CHK_IERR(SUBR(sparseMat_get_range_map)(*vA,&new_map,iflag),*iflag);
  phist_maps_compatible(ctx->range_map,new_map,iflag);
  if (*iflag!=0) 
  {
    PHIST_CHK_IERR(SUBR(sparseMat_delete)(*vA,iflag),*iflag);
    *iflag=-99;
  }
  return;
}

// NOTE: see the description of sparseMat_read_mm on how we treat input flags for this function
extern "C" void SUBR(sparseMat_create_fromRowFunc)(TYPE(sparseMat_ptr) *vA, phist_const_comm_ptr vcomm,
        phist_gidx nrows, phist_gidx ncols, phist_lidx maxnne,
                phist_sparseMat_rowFunc rowFunPtr, void* last_arg,
                int *iflag)
{
  PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFuncWithConstructor)
        (vA, vcomm, nrows, ncols, maxnne, rowFunPtr, nullptr, last_arg, iflag), *iflag);
}

extern "C" void SUBR(sparseMat_create_fromRowFuncWithConstructor)(TYPE(sparseMat_ptr) *A, phist_const_comm_ptr vcomm,
        phist_gidx nrows, phist_gidx ncols, phist_lidx maxnne,
                phist_sparseMat_rowFunc rowFunPtr, 
                phist_sparseMat_rowFuncConstructor rowFunConstructorPtr, 
                void* last_arg, int *iflag)
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
        if( !(*iflag & PHIST_SPARSEMAT_QUIET) )
        {
          PHIST_SOUT(PHIST_INFO,"NOTE: You indicated that the matrix will be used in CARP-CG,\n"
                              "      as it seems that there is only one OpenMP thread, I \n"
                              "      will not construct a coloring, so subsequent CARP sweeps will be sequential per MPI process.\n"
                              "      If you want to use the coloring kernel anyway, specify the flag \n"
                              "       PHIST_SPARSEMAT_OPT_CARP|PHIST_SPARSEMAT_DIST2_COLOR to enforce it.\n"
                          );
        }
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
        rowFunPtr, rowFunConstructorPtr, last_arg, iflag), *iflag);
}



