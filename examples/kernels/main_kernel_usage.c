#include "kernels.h"
#include "macros.h"

int main()
  {
  int ierr;
  DcrsMat_ptr_t A;
  const_map_ptr_t row_map,range_map,domain_map;
  Dmvec_ptr_t x,y;
  
  double *xval, *yval;
  int nloc_x, nloc_y;
  
  const char* filename = "test.mm";
  
  ESSEX_ERROR_HANDLER(Dread_crsMat_mm(&A,filename,&ierr),ierr);
  
  ESSEX_ERROR_HANDLER(Dget_range_map(A, &range_map, &ierr),ierr);
  ESSEX_ERROR_HANDLER(Dget_domain_map(A, &domain_map, &ierr),ierr);

  ESSEX_ERROR_HANDLER(Dcreate_mvec)(domain_map,1,&x,&xval,&ierr),ierr);
  ESSEX_ERROR_HANDLER(Dcreate_mvec)(range_map,1,&y,&yval,&ierr),ierr);
  
  ESSEX_ERROR_HANDLER(get_my_length)(x,&nloc_x,&ierr),ierr);
  ESSEX_ERROR_HANDLER(get_my_length)(y,&nloc_y,&ierr),ierr);

  // set x=1
  ESSEX_ERROR_HANDLER(Dmvec_set_value(x,1.0,&ierr),ierr);

  // compute y=A*x
  ESSEX_ERROR_HANDLER(DcrsMat_X_mvec(1.0,A,x,0.0,y));

  // print result
  for (int i=0;i<nloc;i++)
    {
    fprintf(stdout,"%d\t%16.8g\n",i+1,yval[i]);
    ]
  
  // delete everything
  ESSEX_ERROR_HANDLER(delete_mvec(x,&ierr),ierr);
  ESSEX_ERROR_HANDLER(delete_mvec(y,&ierr),ierr);
  ESSEX_ERROR_HANDLER(delete_crsMat(A,&ierr),ierr);
  

  ESSEX_ERROR_HANDLER  
  
  }
