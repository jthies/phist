#ifdef __cplusplus
extern "C" {
#endif

//! kernels to implement CARP for the matrix A-sigma[j]I, or
//! (A-sigma[j]I, Q; Q', 0], with Q orthonormal. The parallelization
//! and implementation of the `bordering' with Q are left to the
//! kernel lib, it may use e.g. CARP, multi-coloring etc., project out
//! all Q columns at once or one-by-one.

//! create data structures needed for subsequent calls to carp_sweep:
//!
//! input: 
//!
//! numShifts: number of shifts sigma 
//! sigma_r, sigma_i: possibly complex shifts sigma=sigma_r+i*sigma_i
//!
//! output:
//!
//! *work should be NULL on input and not touched while
//! carp_sweep is being called. If no longer needed, it should be cleaned
//! up using carp_destroy.
//!
//! If the shifts sigma change, carp_destroy and carp_setup should be
//! used to rebuild the working objects.
void SUBR(carp_setup)(TYPE(const_sparseMat_ptr) A, int numShifts, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        void** work, int* iflag);

//! forward/backward sweep of Kaczmarz/CARP algorithm (SSOR sweep on the normal equations),
//! with matrix A-sigma[j]*I applied to the j'th column of mvec X with rhs B. The
//! input argument work and shifts sigma_r,sigma_i must be unchanged from the _setup routine. 
//! For each shift (column of X,B)
//! a separate relaxation parameter omega[j] must be provided. In real arithmetic, this function
//! supports a complex shift and X vector, but only a real rhs B. In complex arithmetic, X_i is
//! ignored. If X_i is NULL in real arithmetic, sigma_i is ignored. Otherwise, sigma_i must be provided,
//! even if its values are 0.
void SUBR(carp_sweep)(TYPE(const_sparseMat_ptr) A,
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) Rhs,
        TYPE(mvec_ptr) X_r, TYPE(mvec_ptr) X_i,
        void* const work,
        _MT_ const * omega, int* iflag);

//!                                                                     
//! perform KACZ forward/backward sweep on the `augmented' system       
//!                                                                     
//! | A-sigma[j]I  Q ||X|   |Rhs|                                       
//! | Q'           0 ||q| = |0  |                                       
//!                                                                     
//! where Q'Q=I should be orthonormal. For the moment we only allow a   
//! real-valued Q in real arithmetic (this may change in the future).   
//!                                                                     
void SUBR(carp_sweep_aug)(TYPE(const_sparseMat_ptr) A,
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) Q,
        TYPE(const_mvec_ptr) Rhs,
        TYPE(mvec_ptr) X_r, TYPE(mvec_ptr) X_i,
        TYPE(sdMat_ptr) q_r, TYPE(sdMat_ptr) q_i,
        void* const work,
        _MT_ const * omega, int* iflag);

//! clean up data structures created by carp_setup
void SUBR(carp_destroy)(TYPE(const_sparseMat_ptr) A,
        void* work, int *iflag);


#ifdef __cplusplus
} //extern "C"
#endif


