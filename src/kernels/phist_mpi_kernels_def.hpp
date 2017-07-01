/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include <limits>

// implemnts Isend (flag=0) and Irecv (flag=1)
void SUBR(mvec_transfer)(TYPE(const_mvec_ptr) V, int dest, int tag, 
        MPI_Comm comm, MPI_Request* req, int flag, int* iflag);

extern "C" void SUBR(mvec_Isend)(TYPE(const_mvec_ptr) V, int dest, int tag, 
        MPI_Comm comm, MPI_Request* req, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_OUT(PHIST_DEBUG,"do Isend\n");
  PHIST_CHK_IERR(SUBR(mvec_transfer)(V,dest,tag,comm,req,0,iflag),*iflag);
}

extern "C" void SUBR(mvec_Irecv)(TYPE(const_mvec_ptr) V, int dest, int tag, 
        MPI_Comm comm, MPI_Request* req, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_OUT(PHIST_DEBUG,"do Irecv\n");
  PHIST_CHK_IERR(SUBR(mvec_transfer)(V,dest,tag,comm,req,1,iflag),*iflag);
}


void SUBR(mvec_transfer)(TYPE(const_mvec_ptr) V, int dest, int tag, 
        MPI_Comm comm, MPI_Request* req, int flag,int* iflag)
{
#ifndef PHIST_HAVE_MPI
  *iflag = PHIST_NOT_IMPLEMENTED;
#else
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;
  ST* val;
  phist_lidx lda, nloc;
  int nvec;

  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&nloc,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nvec,iflag),*iflag);

  PHIST_CHK_IERR(SUBR(mvec_extract_view)((TYPE(mvec_ptr))V,&val,&lda,iflag),*iflag);
  int imax=std::numeric_limits<int>::max();

#if 1
  if (nloc*nvec>imax)
  {
    *iflag=PHIST_INTEGER_OVERFLOW;
    return;
  }
#ifdef PHIST_MVECS_ROW_MAJOR
  if (lda!=nvec) {
    PHIST_OUT(PHIST_ERROR,"general case not implemented, need lda (%d)=nvec (%d) right now\n"
                          "(file %s, line %d)\n",(int)lda,nvec,__FILE__,__LINE__);
#else
  if (lda!=nloc) {
    PHIST_OUT(PHIST_ERROR,"general case not implemented, need lda (%d)=nloc (%d) right now\n"
                          "(file %s, line %d)\n",(int)lda,(int)nloc,__FILE__,__LINE__);
#endif
    *iflag=PHIST_INVALID_INPUT;
    return;
  }
  int num_send=(int)(nloc*nvec);
  PHIST_OUT(PHIST_DEBUG,"tag: %d, num elements: %d\n",tag,num_send);
  if (flag==0)
  {
    MPI_Isend(val,num_send,st::mpi_type(),dest,tag,comm,req);
  }
  else
  {
    MPI_Irecv(val,num_send,st::mpi_type(),dest,tag,comm,req);
  }
#else
  int mpi_stride, mpi_chunksize, mpi_count;
  MPI_Datatype mpi_vtype;
  
  //TODO - we need the const-cast here becaue the kernel lib
  //       defines the vector to be non-const. Not sure if
  //       we should change the interface of this function or
  //       of mvec_extract_view

  if (lda>imax || nloc>imax)
  {
    *iflag=PHIST_INTEGER_OVERFLOW;
    return;
  }
  mpi_stride=(int)lda;
#ifdef PHIST_MVECS_ROW_MAJOR
  mpi_chunksize=nvec;
  mpi_count=nloc;
#else
  mpi_chunksize=nloc;
  mpi_count=nvec;
#endif

  MPI_Type_vector(mpi_count,
                  mpi_chunksize,
                  mpi_stride,
                  st::mpi_type(),
                  &mpi_vtype);
  MPI_Type_commit(&mpi_vtype);
  
  if (flag==0)
  {
    MPI_Isend(val,1,mpi_vtype,dest,tag,comm,req);
  }
  else
  {
    MPI_Irecv(val,1,mpi_vtype,dest,tag,comm,req);
  }
  // according to the MPI standard we can free the data type before the
  // communication is complete:
  MPI_Type_free(&mpi_vtype);
#endif
#endif
}

#ifndef PHIST_KERNEL_LIB_GHOST
// synchronize values  of a small dense matrixamong all processes of a given communicator 
// This function is not performace-critical, it is only used for testing purposes (at least
// it should be!)
extern "C" void SUBR(sdMat_sync_values)(TYPE(sdMat_ptr) V, phist_const_comm_ptr comm, int* iflag)
{
#ifndef PHIST_HAVE_MPI
  *iflag = PHIST_SUCCESS;
#else
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=0;

  MPI_Comm mpi_comm=MPI_COMM_WORLD;
  // this function is not yet correctly implemented for Epetra
#ifndef PHIST_KERNEL_LIB_EPETRA
  PHIST_CHK_IERR(phist_comm_get_mpi_comm(comm,&mpi_comm,iflag),*iflag);
#endif
  PHIST_CHK_IERR(SUBR(sdMat_from_device)(V,iflag),*iflag);
  phist_lidx lda;
  int nrows, ncols;
  ST* val;
  PHIST_CHK_IERR(SUBR(sdMat_extract_view)(V,&val,&lda,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(V,&nrows,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(V,&ncols,iflag),*iflag);
  int nr=nrows;
  int nc=ncols;
#ifdef PHIST_SDMATS_ROW_MAJOR
  nc=nrows; nr=ncols;
#endif
  if (nr==lda)
  {
    MPI_Bcast(val,nr*nc,st::mpi_type(),0,mpi_comm);
  }
  else
  {
    for (int j=0;j<nc;j++)
    {
      MPI_Bcast(val+j*lda,nr,st::mpi_type(),0,mpi_comm);
    }
  }
  PHIST_CHK_IERR(SUBR(sdMat_to_device)(V,iflag),*iflag);
#endif
}

#endif
