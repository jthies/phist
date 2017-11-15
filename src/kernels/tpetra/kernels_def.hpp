/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

#include "../common/kernels_no_impl.cpp"


/*
extern "C"
{

void SUBR(sparse_mat_read_mm)(TYPE(sparseMat_ptr)* matrixPtr,
                              phist_const_comm_ptr vcomm,
                              const char* fileName, int* iflag)
{
    PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
    if (fileName == std::nullptr)
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
        resultMatrix = reader.readSparseFile(std::string{fileName}, comm);        
    )

    auto resultMatrixPtr = resultMatrix.release();
    *matrixPtr = static_cast<TYPE(sparseMat_ptr)>(resultMatrixPtr);
}

} // extern "C"
*/



