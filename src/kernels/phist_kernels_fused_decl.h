//! \defgroup 'fused' kernels computing two or more quantities.
//! note: these kernels are very uncommon, kernel libraries that do not have support for some or any of them may include
//! the suboptimal default implementations in common/kernels_no_fused*.cpp
//! We decided to add these operations because they appear in particular algorithms, not because they are of general use.
//! Optimized implementations of these kernels may achieve higher performance by data reuse.
//!@{

#ifdef __cplusplus
extern "C" {
#endif

//! \name fused spmv kernels (sparseMat_times_mvec *and* something else involving x and/or y) \addtogroup crsmat
//@{

//! W=alpha*A*V + beta*W, Wnrm[i] = ||W[i]||_2
void SUBR(sparseMat_times_mvec_fused_norm2)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr)  V,
                                            _ST_ beta,                               TYPE(mvec_ptr)        W,
                                                                                     _MT_*                 Wnrm,
                                            int* iflag);

//! W=alpha*A*V + beta*W, WdotV[i] = W[i]'V[i]
void SUBR(sparseMat_times_mvec_fused_dot)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr)  V,
                                          _ST_ beta,                               TYPE(mvec_ptr)        W,
                                                                                   _ST_*                 WdotV,
                                          int* iflag);

//! W=alpha*A*V + beta*W, WdotV[i] = W[i]'V[i], Wnrm[i] = ||W[i]||_2
void SUBR(sparseMat_times_mvec_fused_dot_norm2)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr)  V,
                                                _ST_ beta,                               TYPE(mvec_ptr)        W,
                                                            _ST_*                 WdotV, _MT_*                 Wnrm,
                                                int* iflag);

//! W=alpha*A*V + beta*W, D = W'W
void SUBR(sparseMat_times_mvec_fused_mvecT_times_mvec_self)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr)  V,
                                                            _ST_ beta,                               TYPE(mvec_ptr)        W,
                                                                                                     TYPE(sdMat_ptr)       D,
                                                            int* iflag);

//! W=alpha*A*V + beta*W, C = W'V
void SUBR(sparseMat_times_mvec_fused_mvecT_times_mvec_other)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr) V,
                                                            _ST_ beta,                                TYPE(mvec_ptr)        W,
                                                                                                      TYPE(sdMat_ptr)       C,
                                                            int* iflag);

//! W=alpha*A*V + beta*W, C = W'V, D = W'W
void SUBR(sparseMat_times_mvec_fused_mvecT_times_mvec_both)(_ST_ alpha, TYPE(const_sparseMat_ptr) A, TYPE(const_mvec_ptr)  V,
                                                            _ST_ beta,                               TYPE(mvec_ptr)        W,
                                                                        TYPE(sdMat_ptr)           C, TYPE(sdMat_ptr)       D,
                                                            int* iflag);
//@}

//! \name augmented kernels with two multi-vectors. \addtogroup mvec
//@{

//! D=alpha*V'*(W*C) + beta*D, W=W*C inplace
void SUBR(mvecT_times_mvec_times_sdMat_inplace)(_ST_ alpha, TYPE(const_mvec_ptr)  V,
                                                            TYPE(mvec_ptr)        W,
                                                            TYPE(const_sdMat_ptr) C,
                                                _ST_ beta,  TYPE(sdMat_ptr)       D,
                                                int* iflag);


//! W=alpha*V*C + beta*W, D=W'W
void SUBR(mvec_times_sdMat_aug)(_ST_ alpha, TYPE(const_mvec_ptr)  V,
                                                  TYPE(const_sdMat_ptr) C,
                                      _ST_ beta,  TYPE(mvec_ptr)        W,
                                                  TYPE(sdMat_ptr)       D,
                                                  int* iflag);

//! augmented kernel with two multi-vectors and two sdMats.
//! W <- = V*C + W*D
void SUBR(mvec_times_sdMat_add_mvec_times_sdMat)(TYPE(const_mvec_ptr) V, 
                                                 TYPE(const_sdMat_ptr) C,
                                                 TYPE(mvec_ptr) W, 
                                                 TYPE(const_sdMat_ptr) D,
                                                 int* iflag);


//@}

#ifdef __cplusplus
} //extern "C"
#endif

//@}
