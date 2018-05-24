/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_macros.h"

//! \addtogroup kernels
//!@{

/*! @file phist_kernels_common_decl.h 
    @brief simple derived functions implemented for all kernel libs centrally

*/

#ifdef __cplusplus
extern "C" {
#endif


//! \name getting data from sparseMat
//!@{

//! get global sparseMat size (number of rows) \ingroup crsmat
void SUBR(sparseMat_global_nrows)(TYPE(sparseMat_ptr) A, phist_gidx* s, int* iflag);

//! get global sparseMat size (number of columns) \ingroup crsmat
void SUBR(sparseMat_global_ncols)(TYPE(sparseMat_ptr) A, phist_gidx* s, int* iflag);

//!@}

//! \name getting data from mvec
//!@{

//! retrieve local length of the vectors in V \ingroup mvec
void SUBR(mvec_my_length)(TYPE(const_mvec_ptr) V, phist_lidx* len, int* iflag);

//! retrieve global length of the vectors in V \ingroup mvec
void SUBR(mvec_global_length)(TYPE(const_mvec_ptr) V, phist_gidx* len, int* iflag);

//! retrieve the comm used for MPI communication in V \ingroup mvec
void SUBR(mvec_get_comm)(TYPE(const_mvec_ptr) V, phist_const_comm_ptr* comm, int* iflag);

//!@}

//! \brief create a new mvec with the same dimensions (number of rows and columns) and
//! distribution (map)  as an existing one. \ingroup mvec
//!
//!
//! The values of the new object are not
//! initialized explicitly, so if you want to clone the vector contents as well,
//! you will have to call mvec_set_block afterwards (or similar).
//!
//! *V must be NULL on input. The function accepts the same input values for iflag as mvec_create 
//! and will simply pass them on to the mvec_create function for the new object, so
//! in order to clone an mvec with exactly the same properties it may  be necessary
//! to manually add flags here like MVEC_REPLICATE_DEVICE_MEM
void SUBR(mvec_clone_shape)(TYPE(mvec_ptr)* V, TYPE(const_mvec_ptr) V_in, int* iflag);

//! \brief "fill" an mvec from a user-provided array. \ingroup mvec
//!
//!
//! If the mvec V has n local rows and k columns (as returned by mvec_my_length and mvec_num_vectors, resp.),
//! then ... <br>
//! If input_row_major==1, then data[i*lda+j], i=0..n-1, j=0..k-1
//! should contain the element to be placed in row i and column j of V. <br>
//! If input_row_major==0, then data[j*lda+i], i=0..n-1, j=0..k-1
//! should contain the element to be placed in row i and column j of V.<br>
//!
//! There is a corresponding function mvec_set_data for the reverse operation.
//!
void SUBR(mvec_set_data)(TYPE(mvec_ptr) V, const _ST_* data, phist_lidx lda, int input_row_major, int* iflag);

#ifdef __cplusplus
} // extern "C"
#endif
