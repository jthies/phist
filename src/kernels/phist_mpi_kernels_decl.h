//! MPI_Isend for an mvec

//! \todo documentation
void SUBR(mvec_Isend)(TYPE(const_mvec_ptr) V, int dest, int tag, 
        MPI_Comm comm, MPI_Request* req, int* iflag);

//! MPI_Irecv for an mvec

//! \todo documentation
void SUBR(mvec_Irecv)(TYPE(const_mvec_ptr) V, int dest, int tag, 
        MPI_Comm comm, MPI_Request* req, int* iflag);

// ghost has its own implementation of this function:
#ifndef PHIST_KERNEL_LIB_GHOST
//! synchronize values  of a small dense matrixamong all processes of a given communicator \ingroup sdmat

//! This function is inteded for testing purposes only and should *not* be used in any 
//! performance-relevant situations. Kernel libs not based on MPI *must* provide their
//! own implementation (as done by GHOST).
void SUBR(sdMat_sync_values)(TYPE(sdMat_ptr) V, const_comm_ptr_t comm, int* iflag);

#endif
