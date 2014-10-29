#include <limits>

// implemnts Isend (flag=0) and Irecv (flag=1)
void SUBR(mvec_transfer)(TYPE(const_mvec_ptr) V, int dest, int tag, 
        MPI_Comm comm, MPI_Request* req, int flag, int* ierr);

extern "C" void SUBR(mvec_Isend)(TYPE(const_mvec_ptr) V, int dest, int tag, 
        MPI_Comm comm, MPI_Request* req, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  PHIST_OUT(PHIST_DEBUG,"do Isend\n");
  PHIST_CHK_IERR(SUBR(mvec_transfer)(V,dest,tag,comm,req,0,ierr),*ierr);
}

extern "C" void SUBR(mvec_Irecv)(TYPE(const_mvec_ptr) V, int dest, int tag, 
        MPI_Comm comm, MPI_Request* req, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  PHIST_OUT(PHIST_DEBUG,"do Irecv\n");
  PHIST_CHK_IERR(SUBR(mvec_transfer)(V,dest,tag,comm,req,1,ierr),*ierr);
}


void SUBR(mvec_transfer)(TYPE(const_mvec_ptr) V, int dest, int tag, 
        MPI_Comm comm, MPI_Request* req, int flag,int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  ST* val;
  lidx_t lda, nloc;
  int nvec;

  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&nloc,ierr),*ierr);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&nvec,ierr),*ierr);

  PHIST_CHK_IERR(SUBR(mvec_extract_view)((TYPE(mvec_ptr))V,&val,&lda,ierr),*ierr);
  int imax=std::numeric_limits<int>::max();

#if 1
  if (nloc*nvec>imax)
  {
    *ierr=PHIST_INTEGER_OVERFLOW;
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
    *ierr=PHIST_INVALID_INPUT;
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
    *ierr=PHIST_INTEGER_OVERFLOW;
    return;
  }
  mpi_stride=(int)lda;
#ifdef MVECS_ROW_MAJOR
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
}