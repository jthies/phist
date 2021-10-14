/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \file phist_kernels_fused_decl.h 
//! \brief fused kernels that combine operations to reduce memory traffic

//! \defgroup fused fused kernels that combine operations to reduce memory traffic
//! \ingroup kernels
//!
//! advanced fused kernels for which default implementations         
//! are available but which may benefit from special attention int he kernel library.   
//!                                                                                     
//! Fused kernels compute two or more quantities in one kernel to allow e.g. for        
//! temporal cache locality. The output may include intermediate steps in a sequence    
//! of operations or independent quantities whose computation requires the same data.   
//! This is one of the forms of pipelining used in PHIST to achieve high node-level     
//! performance in memory-bounded situations.
//!
//! In order to avoid longish subroutine names, we have adopted a slightly
//! different nameing scheme for these operations:
//!
//! kernel                   | becomes short     |   mathematical formula
//! ------------------------ | ----------------- | -----------------------
//! sparseMat_times_mvec     |    => spmv        |   W=alpha*A*V+beta*W
//! mvec_times_sdMat         |    => mvsd        |   W=V*S
//! mvec_times_sdMat_inplace |    => mvsdi       |   V=V*S
//! mvecT_times_mvec         |    => mvTmv       |   S=V'W
//! mvec_dot_mvec            |    => mvdot       |   s[j]=V(:,j)'W(:,j)
//! mvec_add_mvec            |    => mvadd       |   W=alpha*V+beta*W
//!
//!@{

#ifdef __cplusplus
extern "C" {
#endif

//! \brief perform concatenated sparse matrix-vector products with two matrices
//!
//! y = alpha*(shifts1[j]*A1 + shifts2[j]*A2) + beta*Y
void SUBR(fused_spmv_pair)(_ST_ alpha, _ST_ const shifts1[], TYPE(const_sparseMat_ptr) A1,
                                       _ST_ const shifts2[], TYPE(const_sparseMat_ptr) A2,
                                       TYPE(const_mvec_ptr) X,
                           _ST_ beta, TYPE(mvec_ptr) Y, int* iflag);

//! \brief W=alpha*A*V + beta*W_in, WdotW[i] = W(:,i)'W(:,i), VdotW[i]=V(:,i)'W(:,i) with W_in=W on input.
//!
//! understands *iflag=PHIST_NO_GLOBAL_REDUCTION, in which case the dot products
//! are only carried out per MPI process.
//! The arrays VdotW and/or WdotW may be NULL, in which case the corresponding 
//! operation is not performed.
void SUBR(fused_spmv_mvdot)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr)  V,
                            _ST_ beta,                                     TYPE(mvec_ptr)  W,
                            _ST_* WdotW, _ST_* VdotW,
                            int* iflag);

//! \brief W=alpha*A*V + beta*W_in, WtW = W'W, VtW = V'W, with W=W_in on input.
//!
//! understands *iflag=PHIST_NO_GLOBAL_REDUCTION, in which case the inner products
//! are only carried out per MPI process.
//! The sdMats VtW and/or WtW may be NULL, in which case the corresponding 
//! operation is not performed.
//! Furthermore, W=NULL is allowed. This implies beta==0 and will not return the computed vector alpha*A*V
void SUBR(fused_spmv_mvTmv)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr)  V,
                             _ST_ beta,                               TYPE(mvec_ptr)        W,
                             TYPE(sdMat_ptr) WtW, TYPE(sdMat_ptr) VtW,
                             int* iflag);


//! \brief like fused_spmv_mvdot, and compute U=gamma*W+delta*U *after updating W*
//!
//! Any of U, WdotW, VdotW may be NULL and will not be touched in that case.
void SUBR(fused_spmv_mvdot_mvadd)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr)  V,
                                  _ST_ beta,                               TYPE(mvec_ptr)  W,
                                  _ST_ gamma, _ST_ delta,                  TYPE(mvec_ptr)  U,
                            _ST_* lWdotW, _ST_* lVdotW,
                            int* iflag);

//! D=alpha*V'*(W*C) + beta*D, W=W*C inplace
void SUBR(fused_mvsdi_mvTmv)(_ST_ alpha, TYPE(const_mvec_ptr)  V,
                                                            TYPE(mvec_ptr)        W,
                                                            TYPE(const_sdMat_ptr) C,
                                                _ST_ beta,  TYPE(sdMat_ptr)       D,
                                                int* iflag);


//! W=alpha*V*C + beta*W, WtW=W'W
void SUBR(fused_mvsd_mvTmv)(_ST_ alpha, TYPE(const_mvec_ptr)  V,
                                                  TYPE(const_sdMat_ptr) C,
                                      _ST_ beta,  TYPE(mvec_ptr)        W,
                                                  TYPE(sdMat_ptr)       WtW,
                                                  int* iflag);


#ifdef __cplusplus
} // extern "C"
#endif

//@}
