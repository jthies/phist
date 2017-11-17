/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

#include <Kokkos_View.hpp>
#include <Tpetra_MultiVector_def.hpp>

extern "C" void SUBR(type_avail)(int *iflag)
{
  *iflag = PHIST_SUCCESS;
}

extern "C"
{

void SUBR(sparseMat_read_mm)(TYPE(sparseMat_ptr)* matrixPtr,
                              phist_const_comm_ptr vcomm,
                              const char* fileName, int* iflag)
{ 
    PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
    if (fileName == nullptr)
    {
        *iflag = PHIST_INVALID_INPUT;
        return;
    }
    Tpetra::MatrixMarket::Reader<Traits<_ST_>::sparseMat_t> reader{};

    Teuchos::RCP<Traits<_ST_>::sparseMat_t> resultMatrix{};
    PHIST_CAST_PTR_FROM_VOID(const comm_type, comm, vcomm, *iflag);

    auto phist_comm_ptr = Teuchos::rcp(comm, false);

    PHIST_TRY_CATCH
    (
        resultMatrix = reader.readSparseFile(std::string{fileName}, phist_comm_ptr),
        *iflag       
    );

    auto resultMatrixPtr = resultMatrix.release();
    //a  *matrixPtr = static_cast<TYPE(sparseMat_ptr)>(resultMatrixPtr);

    *matrixPtr = (TYPE(sparseMat_ptr))(resultMatrixPtr.get()); 
    *iflag = PHIST_NOT_IMPLEMENTED;
}

} // extern "C"

extern "C" void SUBR(sparseMat_read_bin)(TYPE(sparseMat_ptr)* A, phist_const_comm_ptr comm,
        const char* filename,int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_hb)(TYPE(sparseMat_ptr)* A, phist_const_comm_ptr comm,
        const char* filename,int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_mm_with_context)(TYPE(sparseMat_ptr)* A, phist_const_context_ptr ctx,
        const char* filename,int* iflag)
{
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

extern "C" void SUBR(sparseMat_create_fromRowFunc)(TYPE(sparseMat_ptr) *A, phist_const_comm_ptr comm,
        phist_gidx nrows, phist_gidx ncols, phist_lidx maxnne,
        phist_sparseMat_rowFunc rowFunPtr, void* last_arg, int *iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_create_fromRowFuncAndContext)(TYPE(sparseMat_ptr) *vA, phist_const_context_ptr ctx,
        phist_lidx maxnne,phist_sparseMat_rowFunc rowFunPtr,void* last_arg,
        int *iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}


extern "C" void SUBR(sparseMat_get_row_map)(TYPE(const_sparseMat_ptr) A, phist_const_map_ptr* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t, mat, A, *iflag);
  *map = (phist_const_map_ptr)(mat->getRowMap().get());
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_get_col_map)(TYPE(const_sparseMat_ptr) A, phist_const_map_ptr* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t, mat, A, *iflag);
  *map = (phist_const_map_ptr)(mat->getColMap().get());
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_get_domain_map)(TYPE(const_sparseMat_ptr) A, phist_const_map_ptr* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t, mat, A, *iflag);
  *map = (phist_const_map_ptr)(mat->getDomainMap().get());
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_get_range_map)(TYPE(const_sparseMat_ptr) A, phist_const_map_ptr* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t, mat, A, *iflag);
  *map = (phist_const_map_ptr)(mat->getRangeMap().get());
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_create)(TYPE(mvec_ptr)* vec, 
    phist_const_map_ptr map, phist_lidx nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const map_type, inMap, map, *iflag);
  
  auto map_ptr = Teuchos::rcp(inMap, false);
  auto result = new Traits<_ST_>::mvec_t(map_ptr, nvec);
  *vec = (TYPE(mvec_ptr))(result);

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_create)(TYPE(sdMat_ptr)* mat, 
    int nrows, int ncols, phist_const_comm_ptr comm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const comm_type, localComm, comm, *iflag);

  auto phistCommPtr = localComm == nullptr ? 
        Teuchos::DefaultComm<int>::getDefaultSerialComm(Teuchos::null)
      :
        Teuchos::rcp(localComm, false);

  auto localMap = Teuchos::rcp(new map_type(nrows, 0, phistCommPtr, Tpetra::LocallyReplicated));
  auto matrix = new Traits<_ST_>::mvec_t(localMap, ncols);

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) V, phist_const_map_ptr* map, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, mVec, V, *iflag);
  *map = (phist_const_map_ptr)(mVec->getMap().get());
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) V, int* nvec, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, mVec, V, *iflag);
  *nvec = mVec->getNumVectors();
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) M, int* nrows, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t, sdMat, M, *iflag);
  *nrows = sdMat->getLocalLength();
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) M, int* ncols, int* iflag)
{
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t, sdMat, M, *iflag);
  *ncols = sdMat->getNumVectors();
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_to_mvec)(TYPE(const_mvec_ptr) v_in, TYPE(mvec_ptr) v_out, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}


extern "C" void SUBR(mvec_extract_view)(TYPE(mvec_ptr) V, _ST_** val, phist_lidx* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mVec, V, *iflag);

  auto val_ptr = mVec->getLocalView<Kokkos::DefaultHostExecutionSpace>();
  *val = (_ST_*)(val_ptr.ptr_on_device());

  *lda = mVec->getStride();

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) V, _ST_** val, phist_lidx* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t, sdMat, V, *iflag);

  auto val_ptr = sdMat->getLocalView<Kokkos::DefaultHostExecutionSpace>();
  *val = (_ST_*)(val_ptr.ptr_on_device());

  *lda = sdMat->getStride();

  *iflag = PHIST_SUCCESS;
}

#ifdef PHIST_HIGH_PRECISION_KERNELS
extern "C" void SUBR(sdMat_extract_error)(TYPE(sdMat_ptr) M, _ST_** MC_raw, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}
#endif

extern "C" void SUBR(mvec_view_block)(TYPE(mvec_ptr) V,
    TYPE(mvec_ptr)* Vblock,
    int jmin, int jmax, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) V,
    TYPE(mvec_ptr) Vblock,
    int jmin, int jmax, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_set_block)(TYPE(mvec_ptr) V,
    TYPE(const_mvec_ptr) Vblock,
    int jmin, int jmax, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_view_block)(TYPE(mvec_ptr) M, 
    TYPE(mvec_ptr)* Mblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_get_block)(TYPE(const_mvec_ptr) M, 
    TYPE(mvec_ptr) Mblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) M, 
    TYPE(const_sdMat_ptr) Mblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_delete)(TYPE(sparseMat_ptr) A, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_delete)(TYPE(mvec_ptr) V, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_delete)(TYPE(sdMat_ptr) M, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_put_func)(TYPE(mvec_ptr) V,
        phist_mvec_elemFunc funPtr, void* last_arg, int *iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}


extern "C" void SUBR(sdMat_put_value)(TYPE(mvec_ptr) V, _ST_ value, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

#ifndef PHIST_BUILTIN_RNG
extern "C" void SUBR(mvec_random)(TYPE(mvec_ptr) V, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}
extern "C" void SUBR(sdMat_random)(TYPE(sdMat_ptr) M, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}
#endif
extern "C" void SUBR(mvec_print)(TYPE(const_mvec_ptr) vec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, mvec, vec, *iflag);

  mvec->print(std::cout);

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) mat, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t, sdMat, mat, *iflag);

  sdMat->print(std::cout);

  *iflag = PHIST_SUCCESS;
}

// TODO: Maybe utilize a Kokkos::View and nested parallel_for
//       or has phist some parallelism on its own?
//       Depends on the size if it is worth to use parallelism,
//       so we need to benchmark it to see what is best.
extern "C" void SUBR(sdMat_identity)(TYPE(sdMat_ptr) mat, int* iflag)
{
  #include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  _ST_* raw_values = nullptr; 

  phist_lidx lda;
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(mat, &raw_values, &lda, iflag), *iflag);

  int numRows;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(mat, &numRows, iflag), *iflag);

  int numCols;
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(mat, &numCols, iflag), *iflag);

  // Put ones on the diagonal and zeroes on the off diagonals
  for (int col = 0; col != numCols; ++col)
  {
    for (int row = 0; row != numRows; ++row)
    { // Branch predictor should do well enough here, [[unlikely]] attribute in c++20?
      raw_values[lda * col + row] = col == row ? st::one() : st::zero();
    }
  }

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) vec,
                                 _MT_* vnrm, int* iflag)
{ 
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvec, vec, *iflag);

  int nvec = mvec->getNumVectors();
  Teuchos::ArrayView<_MT_> norms{vnrm, nvec};

  PHIST_TRY_CATCH(mvec->norm2(norms), *iflag);

  *iflag = PHIST_SUCCESS;
}
// TODO: All the functions that call Tpetra::MultiVector::scale()
//       should use Kokkos views instead of Teuchos ArrayViews
extern "C" void SUBR(mvec_normalize)(TYPE(mvec_ptr) vec,
                                     _MT_* vnrm, int* iflag)
{ 
  #include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvec, vec, *iflag);

  int nvec = mvec->getNumVectors();
  Teuchos::ArrayView<_MT_> norms{vnrm, nvec};
  Teuchos::Array<_ST_> scaling{nvec};

  PHIST_TRY_CATCH(mvec->norm2(norms), *iflag);

  for (int idx = 0; idx != nvec; ++idx)
    scaling[idx] = norms[idx] == mt::zero() ? mt::one() : mt::one() / norms[idx];
  
  PHIST_TRY_CATCH(mvec->scale(scaling), *iflag);

  *iflag = PHIST_SUCCESS; 
}

extern "C" void SUBR(mvec_scale)(TYPE(mvec_ptr) vec, 
                                 _ST_ scalar, int* iflag)
{  
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvec, vec, *iflag);
  PHIST_TRY_CATCH(mvec->scale(scalar), *iflag);
}

extern "C" void SUBR(mvec_vscale)(TYPE(mvec_ptr) vec, 
                                  const _ST_* scalar, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, mvec, vec, *iflag);

  int nvec = mvec->getNumVectors();
  Teuchos::ArrayView<_ST_> scalars{(_ST_*)scalar, nvec};

  PHIST_TRY_CATCH(mvec->scale(scalars), *iflag);

  *iflag = PHIST_SUCCESS; 
}

// vecOut = alpha * vecIn + beta * vecOut
extern "C" void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) vecIn,
                                    _ST_ beta, TYPE(mvec_ptr) vecOut, 
                                    int* iflag)
{
  #include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, vec, vecIn, *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, resultVec, vecOut, *iflag);

  // Check if we only want to copy vec into resultVec
  if (alpha == st::one() and beta == st::zero())
    *resultVec = *vec; 
  else
    PHIST_TRY_CATCH(resultVec->update(alpha, *vec, beta), *iflag);

  *iflag = PHIST_SUCCESS;
}

// vecOut[i] = alpha[i] * vecIn[i] + beta * vecOut[i]
extern "C" void SUBR(mvec_vadd_mvec)(const _ST_ alpha[], TYPE(const_mvec_ptr) vecIn,
                                     _ST_ beta, TYPE(mvec_ptr) vecOut, 
                                     int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t, vec, vecIn, *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t, resultVec, vecOut, *iflag);

  // Kokkos::parallel_for?
  for (int idx = 0; idx != vec->getNumVectors(); ++idx)
  {
    PHIST_TRY_CATCH(resultVec->getVectorNonConst(idx)->update(alpha[idx], 
                                                           *vec->getVector(idx), 
                                                           beta),
                    *iflag);
  }
  *iflag = PHIST_SUCCESS;
}


extern "C" void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) A,
                                      _ST_ beta,  TYPE(sdMat_ptr) B, 
                                      int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMatT_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) A,
                                       _ST_ beta,  TYPE(sdMat_ptr) B, 
                                       int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_times_mvec_communicate)(TYPE(const_sparseMat_ptr) A, 
                                                       TYPE(const_mvec_ptr) x, 
                                                       int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
    TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMatT_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, 
    TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_times_mvec_vadd_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A,
        const _ST_ shifts[], TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}


extern "C" void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) v, 
    TYPE(const_mvec_ptr) w, 
    _ST_* s, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
    TYPE(const_sdMat_ptr) C, 
    _ST_ beta, TYPE(mvec_ptr) W, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}


extern "C" void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
    TYPE(const_sdMat_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
    TYPE(const_sdMat_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sdMat_times_sdMatT)(_ST_ alpha, TYPE(const_sdMat_ptr) V, 
                                          TYPE(const_sdMat_ptr) W, 
                              _ST_ beta,        TYPE(sdMat_ptr) C,
                              int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) V, 
    TYPE(const_mvec_ptr) W, 
    _ST_ beta, TYPE(sdMat_ptr) C, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_QR)(TYPE(mvec_ptr) V, TYPE(sdMat_ptr) R, int* iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}

#ifdef PHIST_KERNEL_LIB_BUILTIN
extern "C" void SUBR(mvec_gather_mvecs)(TYPE(mvec_ptr) V, TYPE(const_mvec_ptr) W[], int nblocks, int *iflag)
{
  PHIST_MARK_AS_EXPERIMENTAL(__FUNCTION__)
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_scatter_mvecs)(TYPE(const_mvec_ptr) V, TYPE(mvec_ptr) W[], int nblocks, int *iflag)
{
  PHIST_MARK_AS_EXPERIMENTAL(__FUNCTION__)
  *iflag=PHIST_NOT_IMPLEMENTED;
}
#endif

//! mixed real/complex operation: split mvec into real and imag part.
//! if either reV or imV are NULL, it is not touched.
#ifdef IS_COMPLEX
# ifdef IS_DOUBLE
extern "C" void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, phist_Dmvec* reV, phist_Dmvec* imV, int *iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}
extern "C" void SUBR(mvec_combine)(TYPE(mvec_ptr) V, phist_Dconst_mvec_ptr reV, phist_Dconst_mvec_ptr imV, int *iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}
# else
extern "C" void SUBR(mvec_split)(TYPE(const_mvec_ptr) V, phist_Smvec* reV, phist_Smvec* imV, int *iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}
extern "C" void SUBR(mvec_combine)(TYPE(mvec_ptr) V, phist_Sconst_mvec_ptr reV, phist_Sconst_mvec_ptr imV, int *iflag)
{
  *iflag=PHIST_NOT_IMPLEMENTED;
}
# endif
#endif

#include "../common/kernels_no_inplace_VC.cpp"
#include "../common/kernels_no_VC_add_WD.cpp"
//#include "../common/kernels_no_carp.cpp"
#include "../common/kernels_no_gpu.cpp"
#include "../common/kernels_no_io.cpp"
#include "../common/kernels_no_fused.cpp"
#include "../common/default_context_def.hpp"