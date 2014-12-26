//! MPI_Isend for an mvec

//! \todo documentation
void SUBR(mvec_Isend)(TYPE(const_mvec_ptr) V, int dest, int tag, 
        MPI_Comm comm, MPI_Request* req, int* iflag);

//! MPI_Irecv for an mvec

//! \todo documentation
void SUBR(mvec_Irecv)(TYPE(const_mvec_ptr) V, int dest, int tag, 
        MPI_Comm comm, MPI_Request* req, int* iflag);
