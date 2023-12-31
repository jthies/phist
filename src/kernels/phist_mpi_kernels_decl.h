/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \file phist_mpi_kernels_decl.h 
//! \brief communication operations implemented for all kernel libraries

//! \addtogroup mpi_kernels
//!@{

//! \brief MPI_Isend for sending the local part of an mvec to another process
void SUBR(mvec_Isend)(TYPE(const_mvec_ptr) V, int dest, int tag, 
        MPI_Comm comm, MPI_Request* req, int* iflag);

//! \brief MPI_Irecv for receiving the local part of an mvec from another process
void SUBR(mvec_Irecv)(TYPE(const_mvec_ptr) V, int dest, int tag, 
        MPI_Comm comm, MPI_Request* req, int* iflag);

//! \brief synchronize values  of a small dense matrix among all processes of a given communicator \ingroup sdmat
//!
//!
//! This function is intended for testing purposes only and should *not* be used in any 
//! performance-relevant situations. Kernel libs not based on MPI *must* provide their
//! own implementation (as done by GHOST).
void SUBR(sdMat_sync_values)(TYPE(sdMat_ptr) V, phist_const_comm_ptr comm, int* iflag);

//!@}