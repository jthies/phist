#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "phist_typedefs.h"
#include "phist_kernels.h"
#include "phist_operator.h"
#include "phist_macros.h"
#include "phist_lanczos.h"
#include <stdio.h>

#include "phist_gen_z.h"

int main(int argc, char** argv)
  {
  int rank, num_proc;
  int i, ierr;
  bool verbose;

  comm_ptr_t comm;
  ZcrsMat_ptr_t A;
  Zop_ptr_t A_op; // this is a wrapper for the CRS matrix which we pass to the actual solver
  
  const_map_ptr_t map; // map (element distribution) of vectors according to 
                       // the distribution of matrix rows
  Zmvec_ptr_t X; // multivector for getting the eigenvectors
  
  _MT_* evals; // Lanczos assumes a Hermitian matrix, so the eigenvectors are real
  _MT_* resid;

  _MT_ tol;     
  char* filename;
  
  int num_eigs,num_iters,max_iters;
  
  _PHIST_ERROR_HANDLER_(phist_kernels_init(&argc,&argv,&ierr),ierr);

  _PHIST_ERROR_HANDLER_(phist_comm_create(&comm,&ierr),ierr);

  _PHIST_ERROR_HANDLER_(phist_comm_get_rank(comm, &rank, &ierr),ierr);
  _PHIST_ERROR_HANDLER_(phist_comm_get_size(comm, &num_proc, &ierr),ierr);

  verbose= (rank==0);


  if (argc<2)
    {
    if (verbose) fprintf(stderr,"Usage: ./main_zLanczos <matrix market filename> [<num eigs>] [<tol>] [<max iters>]\n");
    return -1;
    }

  filename = argv[1];
  
  if (argc<3)
    {
    num_eigs=5;
    }
  else
    {
    num_eigs=atoi(argv[2]);
    }

  if (argc<4)
    {
    tol=1.0e-6;
    }
  else
    {
    tol=(_MT_)atof(argv[3]);
    }

  if (argc<5)
    {
    max_iters=250;
    }
  else
    {
    max_iters=atoi(argv[4]);
    }

  if (verbose) fprintf(stdout,"looking for %d Eigenpairs of largest magnitude\n",num_eigs);
  
  _PHIST_ERROR_HANDLER_(phist_ZcrsMat_read_mm(&A,filename,&ierr),ierr);
  
  _PHIST_ERROR_HANDLER_(phist_ZcrsMat_get_domain_map(A, &map, &ierr),ierr);

  _PHIST_ERROR_HANDLER_(phist_Zmvec_create(&X,map,num_eigs,&ierr),ierr);

  num_iters=max_iters;
  // create operator wrapper for computing Y=A*X using a CRS matrix
  A_op = (_TYPE_(op_ptr))malloc(sizeof(_TYPE_(op)));
  _PHIST_ERROR_HANDLER_(phist_Zop_wrap_crsMat(A_op,A,&ierr),ierr);
  
  // allocate memory for eigenvalues and residuals
  evals = (_ST_*)malloc(num_eigs*sizeof(_ST_));
  resid = (_MT_*)malloc(num_eigs*sizeof(_MT_));

  _PHIST_ERROR_HANDLER_(phist_Zlanczos(A_op,X,evals, 
        resid, tol,&num_iters,&num_eigs,&ierr),ierr);

  if (verbose)
    {
    fprintf(stdout,"Found %d eigenpairs after %d iterations\n",num_eigs,num_iters);
    }

  if (verbose && num_eigs>0)
    {
    fprintf(stdout,"Eigenvalue \t Ritz Residual \n");
    for (i=0;i<num_eigs;i++)
      {
      fprintf(stdout,"%16.8e\t%4.2e\n",evals[i],resid[i]);
      }
    }

  free(evals);
  free(resid);
  
  _PHIST_ERROR_HANDLER_(phist_ZcrsMat_delete(A,&ierr),ierr);
  _PHIST_ERROR_HANDLER_(phist_Zmvec_delete(X,&ierr),ierr);
  free(A_op);
  _PHIST_ERROR_HANDLER_(phist_kernels_finalize(&ierr),ierr);
  }
