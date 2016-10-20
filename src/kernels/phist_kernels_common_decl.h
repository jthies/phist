#include "phist_macros.h"

//! \addtogroup kernels
//@{

/*! @file phist_kernels_common_decl.h 
    @brief simple derived functions implemented for all kernel libs centrally

*/

#ifdef __cplusplus
extern "C" {
#endif


//! \name getting data from objects
//@{

//! get global sparseMat size (number of rows) \ingroup crsmat
void SUBR(sparseMat_global_nrows)(TYPE(sparseMat_ptr) A, phist_gidx* s, int* iflag);

//! get global sparseMat size (number of columns) \ingroup crsmat
void SUBR(sparseMat_global_ncols)(TYPE(sparseMat_ptr) A, phist_gidx* s, int* iflag);

//! retrieve local length of the vectors in V \ingroup mvec
void SUBR(mvec_my_length)(TYPE(const_mvec_ptr) V, phist_lidx* len, int* iflag);

//! retrieve global length of the vectors in V \ingroup mvec
void SUBR(mvec_global_length)(TYPE(const_mvec_ptr) V, phist_gidx* len, int* iflag);

//! retrieve the comm used for MPI communication in V \ingroup mvec
void SUBR(mvec_get_comm)(TYPE(const_mvec_ptr) V, phist_const_comm_ptr* comm, int* iflag);

//@}

//! create a new mvec with the same dimensions (number of rows and columns) and
//! distribution (map)  as an existing one. The values of the new object are not
//! initialized explicitly, so if you want to clone the vector contents as well,
//! you will have to call mvec_set_block afterwards (or similar). *V must be NULL
//! on input.
void SUBR(mvec_clone_shape)(TYPE(mvec_ptr)* V, TYPE(const_mvec_ptr) V_in, int* iflag);

//! y[i]=alpha*(A*x+shift*x) + beta*y
void SUBR(sparseMat_times_mvec_add_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A,
        _ST_ shift, TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag);

#ifdef __cplusplus
} // extern "C"
#endif
