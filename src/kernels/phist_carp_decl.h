#ifdef __cplusplus
extern "C" {
#endif

//! kernels to implement CARP for the matrix sigma[j]I-A.

//! create data structures needed for subsequent calls to carp_sweep: 
//! 
//! input: 
//!
//! numShifts: number of shifts sigma 
//! sigma_r, sigma_i: possibly complex shifts sigma=sigma_r+i*sigma_i
//!
//! output:
//! 
//! *nrms_ai2i and *work should be NULL on input and not touched while
//! carp_sweep is being called. If no longer needed, they should be cleaned
//! up using carp_destroy.
//!
//! If the shifts sigma change, carp_destroy and carp_setup should be
//! used to rebuild the working objects.
void SUBR(carp_setup)(TYPE(const_crsMat_ptr) A, int numShifts, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        _MT_ **nrms_ai2i, void** work, int* ierr);

//! forward/backward sweep of Kaczmarz/CARP algorithm (SSOR sweep on the normal equations),
//! with matrix sigma[j]*I-A applied to the columns of mvec X[j] with a single rhs B. The
//! input arguments nrms_ai2i and work must be unchanged from the _setup routine. For each
//! shift sigma[j], a separate relaxation parameter omega[j] must be provided.
void SUBR(carp_sweep)(TYPE(const_crsMat_ptr) A, int numShifts, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) Rhs,
        TYPE(mvec_ptr) X_r[], TYPE(mvec_ptr) X_i[],
        _MT_ const* nrm_ai2i, void* const work,
        _MT_ const * omega, int* ierr);

//! clean up data structures created by carp_setup
void SUBR(carp_destroy)(TYPE(const_crsMat_ptr) A, int numShifts,
        _MT_* nrms_ai2i, void* work, int *ierr);


#ifdef __cplusplus
} //extern "C"
#endif


