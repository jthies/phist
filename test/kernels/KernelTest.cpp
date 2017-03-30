/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "KernelTest.h"

#ifndef PHIST_HAVE_MPI
#define MPI_COMM_NULL 0
#endif
phist_comm_ptr KernelTest::comm_ = NULL;
bool KernelTest::haveS_ = false, KernelTest::haveD_ = false, KernelTest::haveC_ = false, KernelTest::haveZ_ = false;
MPI_Comm KernelTest::mpi_comm_ = MPI_COMM_NULL;
unsigned int KernelTest::rseed_ = 0;
int KernelTest::iflag_ = 0, KernelTest::mpi_rank_ = 0, KernelTest::mpi_size_ = 0;
bool KernelTest::isCuda_=false;

int KernelTest::staticKernelTestSetupCounter_ = 0;
