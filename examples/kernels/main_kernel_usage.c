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
  int ierr;

  const_comm_ptr_t comm;
  DcrsMat_ptr_t A;
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
  
  comm_ptr_t comm_world;
  
  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&ierr),ierr);

  PHIST_ICHK_IERR(phist_comm_create(&comm_world,&ierr),ierr);
#ifdef PHIST_KERNEL_LIB_GHOST
  PHIST_ICHK_IERR(phist_DcrsMat_read_bin(&A,filename,comm,&ierr),ierr);
  ghost_sparsemat_t* A_ghost = (ghost_sparsemat_t*)A;
  char *str;
  ghost_sparsemat_string(&str,A_ghost);
  printf("%s\n",str);
  free(str); str = NULL;
#else
  PHIST_ICHK_IERR(phist_DcrsMat_read_mm(&A,filename,comm,&ierr),ierr);
#endif  
  PHIST_ICHK_IERR(phist_DcrsMat_get_range_map(A, &range_map, &ierr),ierr);
  PHIST_ICHK_IERR(phist_DcrsMat_get_domain_map(A, &domain_map, &ierr),ierr);

  PHIST_ICHK_IERR(phist_map_get_comm(range_map, &comm, &ierr),ierr);
  PHIST_ICHK_IERR(phist_comm_get_rank(comm, &rank, &ierr),ierr);
  PHIST_ICHK_IERR(phist_comm_get_size(comm, &num_proc, &ierr),ierr);

  PHIST_ICHK_IERR(phist_Dmvec_create(&x,domain_map,1,&ierr),ierr);
  PHIST_ICHK_IERR(phist_Dmvec_create(&y,range_map,1,&ierr),ierr);

  PHIST_OUT(3,"x pointer @ %p",x);
  PHIST_OUT(3,"y pointer @ %p",y);

  PHIST_ICHK_IERR(phist_Dmvec_my_length(x,&nloc_x,&ierr),ierr);
  PHIST_ICHK_IERR(phist_Dmvec_my_length(y,&nloc_y,&ierr),ierr);

  PHIST_ICHK_IERR(phist_Dmvec_num_vectors(x,&nvec_x,&ierr),ierr);
  PHIST_ICHK_IERR(phist_Dmvec_num_vectors(y,&nvec_y,&ierr),ierr);
  
  fprintf(stdout,"rank %d: x has local length %d and %d vectors\n",rank,(int)nloc_x,(int)nvec_x);
  //fprintf(stdout,"rank %d: y has local length %d and %d vectors\n",rank,(int)nloc_y,(int)nvec_y);
  
  PHIST_ICHK_IERR(phist_Dmvec_extract_view(x,&x_val,&lda_x,&ierr),ierr);
  PHIST_ICHK_IERR(phist_Dmvec_extract_view(y,&y_val,&lda_y,&ierr),ierr);
  
  // set x=1
  PHIST_ICHK_IERR(phist_Dmvec_put_value(x,1.0,&ierr),ierr);
  PHIST_ICHK_IERR(phist_Dmvec_put_value(y,-99.0,&ierr),ierr);

  // print initial vectors
  PHIST_OUT(1,"x vector data @ %p (lda %d)",x_val,lda_x);
  PHIST_OUT(1,"y vector data @ %p (lda %d)",y_val,lda_y);
  
  for (i=0;i<nloc_y;i++)
    {
    PHIST_OUT(0,"%d\t%16.8g\t%16.8g",i+1,x_val[i],y_val[i]);
    }

  // compute y=A*x
  PHIST_ICHK_IERR(phist_DcrsMat_times_mvec(1.0,A,x,0.0,y,&ierr),ierr);

  // print result
  PHIST_OUT(1,"after MVM:");
  for (i=0;i<nloc_y;i++)
    {
    PHIST_OUT(0,"%d\t%16.8g\t%16.8g",i+1,x_val[i],y_val[i]);
    }
  
  // delete everything
  PHIST_ICHK_IERR(phist_Dmvec_delete(x,&ierr),ierr);
  PHIST_ICHK_IERR(phist_Dmvec_delete(y,&ierr),ierr);
  PHIST_ICHK_IERR(phist_DcrsMat_delete(A,&ierr),ierr);

  PHIST_ICHK_IERR(phist_kernels_finalize(&ierr),ierr);
  }
