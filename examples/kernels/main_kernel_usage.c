#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "phist_typedefs.h"
#include "phist_kernels.h"
#include "phist_macros.h"
#include <stdio.h>

#ifdef PHIST_KERNEL_LIB_GHOST
#include "ghost.h" 
#endif

int main(int argc, char** argv)
  {
  int rank, num_proc;
  int iflag;

  const_comm_ptr_t comm = NULL;
  DsparseMat_ptr_t A;
  const_map_ptr_t row_map,range_map,domain_map;
  Dmvec_ptr_t x,y;
  
  double *x_val, *y_val;
  lidx_t nloc_x, nloc_y;
  int nvec_x,nvec_y;
  lidx_t lda_x, lda_y;
#ifdef PHIST_KERNEL_LIB_GHOST   
  const char* filename = "test.crs";
#else
  const char* filename = "test.mm";
#endif  
  int i;
  
  comm_ptr_t comm_world = NULL;
  
  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&iflag),iflag);

  PHIST_ICHK_IERR(phist_comm_create(&comm_world,&iflag),iflag);
  PHIST_ICHK_IERR(phist_DsparseMat_read_mm(&A,filename,comm_world,&iflag),iflag);
  PHIST_ICHK_IERR(phist_DsparseMat_get_range_map(A, &range_map, &iflag),iflag);
  PHIST_ICHK_IERR(phist_DsparseMat_get_domain_map(A, &domain_map, &iflag),iflag);

  PHIST_ICHK_IERR(phist_map_get_comm(range_map, &comm, &iflag),iflag);
  PHIST_ICHK_IERR(phist_comm_get_rank(comm, &rank, &iflag),iflag);
  PHIST_ICHK_IERR(phist_comm_get_size(comm, &num_proc, &iflag),iflag);

  PHIST_ICHK_IERR(phist_Dmvec_create(&x,domain_map,1,&iflag),iflag);
  PHIST_ICHK_IERR(phist_Dmvec_create(&y,range_map,1,&iflag),iflag);

  PHIST_OUT(3,"x pointer @ %p",x);
  PHIST_OUT(3,"y pointer @ %p",y);

  PHIST_ICHK_IERR(phist_Dmvec_my_length(x,&nloc_x,&iflag),iflag);
  PHIST_ICHK_IERR(phist_Dmvec_my_length(y,&nloc_y,&iflag),iflag);

  PHIST_ICHK_IERR(phist_Dmvec_num_vectors(x,&nvec_x,&iflag),iflag);
  PHIST_ICHK_IERR(phist_Dmvec_num_vectors(y,&nvec_y,&iflag),iflag);
  
  fprintf(stdout,"rank %d: x has local length %d and %d vectors\n",rank,(int)nloc_x,(int)nvec_x);
  //fprintf(stdout,"rank %d: y has local length %d and %d vectors\n",rank,(int)nloc_y,(int)nvec_y);
  
  PHIST_ICHK_IERR(phist_Dmvec_extract_view(x,&x_val,&lda_x,&iflag),iflag);
  PHIST_ICHK_IERR(phist_Dmvec_extract_view(y,&y_val,&lda_y,&iflag),iflag);
  
  // set x=1
  PHIST_ICHK_IERR(phist_Dmvec_put_value(x,1.0,&iflag),iflag);
  PHIST_ICHK_IERR(phist_Dmvec_put_value(y,-99.0,&iflag),iflag);

  // print initial vectors
  PHIST_OUT(1,"x vector data @ %p (lda %d)",x_val,lda_x);
  PHIST_OUT(1,"y vector data @ %p (lda %d)",y_val,lda_y);
  
  for (i=0;i<nloc_y;i++)
    {
    PHIST_OUT(0,"%d\t%16.8g\t%16.8g",i+1,x_val[i],y_val[i]);
    }

  // compute y=A*x
  PHIST_ICHK_IERR(phist_DsparseMat_times_mvec(1.0,A,x,0.0,y,&iflag),iflag);

  // print result
  PHIST_OUT(1,"after MVM:");
  for (i=0;i<nloc_y;i++)
    {
    PHIST_OUT(0,"%d\t%16.8g\t%16.8g",i+1,x_val[i],y_val[i]);
    }
  
  // delete everything
  PHIST_ICHK_IERR(phist_Dmvec_delete(x,&iflag),iflag);
  PHIST_ICHK_IERR(phist_Dmvec_delete(y,&iflag),iflag);
  PHIST_ICHK_IERR(phist_DsparseMat_delete(A,&iflag),iflag);

  PHIST_ICHK_IERR(phist_kernels_finalize(&iflag),iflag);
  }
