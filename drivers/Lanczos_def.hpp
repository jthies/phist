
using namespace phist;
/*
typedef ScalarTraits< _ST_ >::scaLar_t ST;
typedef ScalarTraits< _ST_ >::magn_t MT;
*/

typedef _ST_ ST;
typedef _MT_ MT;
typedef _TYPE_(mvec_ptr) mvec_ptr_t;
typedef _TYPE_(const_mvec_ptr) const_mvec_ptr_t;

typedef _TYPE_(sdMat_ptr) sdMat_ptr_t;
typedef _TYPE_(const_sdMat_ptr) const_sdMat_ptr_t;

typedef _TYPE_(crsMat_ptr) crsMat_ptr_t;
typedef _TYPE_(const_crsMat_ptr) const_crsMat_ptr_t;

typedef _TYPE_(op_ptr) op_ptr_t;
typedef _TYPE_(const_op_ptr) const_op_ptr_t;

int main(int argc, char** argv)
  {
  int rank, num_proc;
  int i, ierr;
  bool verbose;

  comm_ptr_t comm;
  crsMat_ptr_t A;
  op_ptr_t A_op; // this is a wrapper for the CRS matrix which we pass to the actual solver
  
  const_map_ptr_t map; // map (element distribution) of vectors according to 
                       // the distribution of matrix rows
  mvec_ptr_t X; // multivector for getting the eigenvectors
  
  MT* evals; // Lanczos assumes a Hermitian matrix, so the eigenvectors are real
  MT* resid;

  MT tol;     
  char* filename;
  
  int num_eigs,num_iters,max_iters;
  
  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&ierr),ierr);

  PHIST_ICHK_IERR(phist_comm_create(&comm,&ierr),ierr);

  PHIST_ICHK_IERR(phist_comm_get_rank(comm, &rank,&ierr),ierr);
  PHIST_ICHK_IERR(phist_comm_get_size(comm, &num_proc,&ierr),ierr);

  verbose= (rank==0);


  if (argc<2)
    {
    if (verbose) std::cout << "Usage: ./main_zLanczos <matrix market filename> [<num eigs>] [<tol>] [<max iters>]"<<std::endl;;
    return 1;
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
    tol=(MT)atof(argv[3]);
    }

  if (argc<5)
    {
    max_iters=250;
    }
  else
    {
    max_iters=atoi(argv[4]);
    }

  if (verbose) std::cout << "looking for "<<num_eigs<<" Eigenpairs of largest magnitude"<<std::endl;
  
  PHIST_ICHK_IERR(_SUBR_(crsMat_read_mm)(&A,filename,&ierr),ierr);
  
  PHIST_ICHK_IERR(_SUBR_(crsMat_get_domain_map)(A, &map,&ierr),ierr);

  PHIST_ICHK_IERR(_SUBR_(mvec_create)(&X,map,num_eigs,&ierr),ierr);

  num_iters=max_iters;
  // create operator wrapper for computing Y=A*X using a CRS matrix
  A_op = (_TYPE_(op_ptr))malloc(sizeof(_TYPE_(op)));
  PHIST_ICHK_IERR(_SUBR_(op_wrap_crsMat)(A_op,A,&ierr),ierr);
  
  // allocate memory for eigenvalues and residuals
  evals = (MT*)malloc(num_eigs*sizeof(MT));
  resid = (MT*)malloc(num_eigs*sizeof(MT));

  _SUBR_(lanczos)(A_op,X,evals, 
        resid, tol,&num_iters,&num_eigs,&ierr);
  if (ierr!=0)
    {
    if (verbose) fprintf(stdout,"code %d returned from Lanczos",ierr);
    if (ierr<0) return ierr;
    }
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
  
  PHIST_ICHK_IERR(_SUBR_(crsMat_delete)(A,&ierr),ierr);
  PHIST_ICHK_IERR(_SUBR_(mvec_delete)(X,&ierr),ierr);
  free(A_op);
  PHIST_ICHK_IERR(phist_kernels_finalize(&ierr),ierr);
  return ierr;
  }
