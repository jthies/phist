#include <limits>

extern "C" void SUBR(mvec_Isend)(TYPE(const_mvec_ptr) V, int dest, int tag, 
        MPI_Comm comm, MPI_Request* req, int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  ST* val;
  lidx_t lda, nloc;
  int nvec;
  int mpi_stride, mpi_chunksize, mpi_count;
  MPI_Datatype mpi_vtype;
  
  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&nloc,ierr),*ierr);
  //TODO - we need the const-cast here becaue the kernel lib
  //       defines the vector to be non-const. Not sure if
  //       we should change the interface of this function or
  //       of mvec_extract_view
  PHIST_CHK_IERR(SUBR(mvec_extract_view)((TYPE(mvec_ptr))V,&val,&lda,ierr),*ierr);
  int imax=std::numeric_limits<int>::max();
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
  
  MPI_Isend(val,1,mpi_vtype,dest,tag,comm,req);
  // according to the MPI standard we can free the data type before the
  // communication is complete:
  MPI_Type_free(&mpi_vtype);
}

extern "C" void SUBR(mvec_Irecv)(TYPE(const_mvec_ptr) V, int dest, int tag, 
        MPI_Comm comm, MPI_Request* req, int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  ST* val;
  lidx_t lda, nloc;
  int nvec;
  int mpi_stride, mpi_chunksize, mpi_count;
  MPI_Datatype mpi_vtype;
  
  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&nloc,ierr),*ierr);
  //TODO - we need the const-cast here becaue the kernel lib
  //       defines the vector to be non-const. Not sure if
  //       we should change the interface of this function or
  //       of mvec_extract_view
  PHIST_CHK_IERR(SUBR(mvec_extract_view)((TYPE(mvec_ptr))V,&val,&lda,ierr),*ierr);
  int imax=std::numeric_limits<int>::max();
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
  
  MPI_Irecv(val,1,mpi_vtype,dest,tag,comm,req);
  // according to the MPI standard we can free the data type before the
  // communication is complete:
  MPI_Type_free(&mpi_vtype);
}
