#include "kernels/phist_kernels.h"
#include "tools/phist_macros.h"
#include <stdio.h>

int main(int argc, char** argv)
  {
  int rank, num_proc;
  int ierr;

  const_comm_ptr_t comm;
  DcrsMat_ptr_t A;
  const_map_ptr_t row_map,range_map,domain_map;
  Dmvec_ptr_t x,y;
  
  double *x_val, *y_val;
  int nloc_x, nloc_y;
  int nvec_x,nvec_y;
  int lda_x, lda_y;
   
  const char* filename = "test.mm";
  
  int i;
  
  comm_ptr_t comm_world;
  
  _PHIST_ERROR_HANDLER_(phist_kernels_init(&argc,&argv,&ierr),ierr);

  _PHIST_ERROR_HANDLER_(phist_comm_create(&comm_world,&ierr),ierr);

  _PHIST_ERROR_HANDLER_(phist_DcrsMat_read_mm(&A,filename,&ierr),ierr);
  
  _PHIST_ERROR_HANDLER_(phist_DcrsMat_get_range_map(A, &range_map, &ierr),ierr);
  _PHIST_ERROR_HANDLER_(phist_DcrsMat_get_domain_map(A, &domain_map, &ierr),ierr);

  _PHIST_ERROR_HANDLER_(phist_map_get_comm(range_map, &comm, &ierr),ierr);
  _PHIST_ERROR_HANDLER_(phist_comm_get_rank(comm, &rank, &ierr),ierr);
  _PHIST_ERROR_HANDLER_(phist_comm_get_size(comm, &num_proc, &ierr),ierr);

  _PHIST_ERROR_HANDLER_(phist_Dmvec_create(&x,domain_map,1,&ierr),ierr);
  _PHIST_ERROR_HANDLER_(phist_Dmvec_create(&y,range_map,1,&ierr),ierr);

  _PHIST_ERROR_HANDLER_(phist_Dmvec_my_length(x,&nloc_x,&ierr),ierr);
  _PHIST_ERROR_HANDLER_(phist_Dmvec_my_length(y,&nloc_y,&ierr),ierr);

  _PHIST_ERROR_HANDLER_(phist_Dmvec_num_vectors(x,&nvec_x,&ierr),ierr);
  _PHIST_ERROR_HANDLER_(phist_Dmvec_num_vectors(y,&nvec_y,&ierr),ierr);
  
  fprintf(stdout,"rank %d: x has local length %d and %d vectors\n",rank,nloc_x,nvec_x);
  fprintf(stdout,"rank %d: y has local length %d and %d vectors\n",rank,nloc_y,nvec_y);
  


  _PHIST_ERROR_HANDLER_(phist_Dmvec_extract_view(x,&x_val,&lda_x,&ierr),ierr);
  _PHIST_ERROR_HANDLER_(phist_Dmvec_extract_view(y,&y_val,&lda_y,&ierr),ierr);
  
  // set x=1
  _PHIST_ERROR_HANDLER_(phist_Dmvec_put_value(x,42.0,&ierr),ierr);
  _PHIST_ERROR_HANDLER_(phist_Dmvec_put_value(y,-99.0,&ierr),ierr);

  // compute y=A*x
  _PHIST_ERROR_HANDLER_(phist_DcrsMat_times_mvec(1.0,A,x,0.0,y,&ierr),ierr);

  // print result
  for (i=0;i<nloc_y;i++)
    {
    fprintf(stdout,"%d\t%16.8g\t%16.8g\n",i+1,x_val[i],y_val[i]);
    }
  
  // delete everything
  _PHIST_ERROR_HANDLER_(phist_Dmvec_delete(x,&ierr),ierr);
  _PHIST_ERROR_HANDLER_(phist_Dmvec_delete(y,&ierr),ierr);
  _PHIST_ERROR_HANDLER_(phist_DcrsMat_delete(A,&ierr),ierr);

  _PHIST_ERROR_HANDLER_(phist_kernels_finalize(&ierr),ierr);
  }
