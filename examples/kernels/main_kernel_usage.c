#include "kernels/essex_kernels.h"
#include "tools/essex_macros.h"
#include "../../test/helpers/essex_test_helpers.h"
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
  
  const char* filename = "test.mm";
  
  int i;
  
  essex_bad_object* troet;
  comm_ptr_t comm_world;
  
  _ESSEX_ERROR_HANDLER_(essex_kernels_init(&argc,&argv,&ierr),ierr);

  _ESSEX_ERROR_HANDLER_(essex_comm_create(&comm_world,&ierr),ierr);

  _ESSEX_ERROR_HANDLER_(essex_DcrsMat_read_mm(&A,filename,&ierr),ierr);
  
  _ESSEX_ERROR_HANDLER_(essex_DcrsMat_get_range_map(A, &range_map, &ierr),ierr);
  _ESSEX_ERROR_HANDLER_(essex_DcrsMat_get_domain_map(A, &domain_map, &ierr),ierr);

  _ESSEX_ERROR_HANDLER_(essex_map_get_comm(range_map, &comm, &ierr),ierr);
  _ESSEX_ERROR_HANDLER_(essex_comm_get_rank(comm, &rank, &ierr),ierr);
  _ESSEX_ERROR_HANDLER_(essex_comm_get_size(comm, &num_proc, &ierr),ierr);

  _ESSEX_ERROR_HANDLER_(essex_Dmvec_create(domain_map,1,&x,&ierr),ierr);
  _ESSEX_ERROR_HANDLER_(essex_Dmvec_create(range_map,1,&y,&ierr),ierr);

  _ESSEX_ERROR_HANDLER_(essex_Dmvec_my_length(x,&nloc_x,&ierr),ierr);
  _ESSEX_ERROR_HANDLER_(essex_Dmvec_my_length(y,&nloc_y,&ierr),ierr);
  
  fprintf(stdout,"rank %d: x has local length %d\n",rank,nloc_x);
  fprintf(stdout,"rank %d: y has local length %d\n",rank,nloc_y);
  


  _ESSEX_TEST_HANDLER_(essex_Dmvec_extract_view(troet,&x_val,0,&ierr),ierr,_ESSEX_BAD_CAST_);
  _ESSEX_ERROR_HANDLER_(essex_Dmvec_extract_view(x,&x_val,0,&ierr),ierr);
  _ESSEX_ERROR_HANDLER_(essex_Dmvec_extract_view(y,&y_val,0,&ierr),ierr);
  
  // set x=1
  _ESSEX_ERROR_HANDLER_(essex_Dmvec_put_value(x,1.0,&ierr),ierr);

  // compute y=A*x
  _ESSEX_ERROR_HANDLER_(essex_DcrsMat_X_mvec(1.0,A,x,0.0,y,&ierr),ierr);

  // print result
  for (i=0;i<nloc_y;i++)
    {
    fprintf(stdout,"%d\t%16.8g\n",i+1,y_val[i]);
    }
  
  // delete everything
  _ESSEX_ERROR_HANDLER_(essex_Dmvec_delete(x,&ierr),ierr);
  _ESSEX_ERROR_HANDLER_(essex_Dmvec_delete(y,&ierr),ierr);
  _ESSEX_ERROR_HANDLER_(essex_DcrsMat_delete(A,&ierr),ierr);
  

  _ESSEX_ERROR_HANDLER_(essex_kernels_finalize(&ierr),ierr);
  }
