#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "phist_macros.h"
#include "phist_enums.h"
#include "phist_kernels.h"
#include "phist_operator.h"
#include "phist_jdqr.h"
#include "phist_jadaOpts.h"

// ghost/spinChain stuff
#include "ghost.h"
#include "ghost/util.h"
#include "matfuncs.h"

GHOST_REGISTER_DT_D(my_datatype)

#include "phist_gen_d.h"
#include "phist_driver_utils.h"

int main(int argc, char** argv)
  {
     if ( argc < 2 )
      {
          printf("usage: spinChain <L> <howMany> <tol> <maxIter> <minBas> <maxBas> "
                 "                 <arno> <initShift>\n"
                 " where: L:         number of spins \n"
                 "        howMany:   number of (right-most) eigenvalues to seek\n"
                 "        tol:       convergence tolerance\n"
                 "        maxIter:   max. num JaDa iterations\n"
                 "        minBas:    number of vectors from which to restart\n"
                 "        maxBas:    max num vectors to create before restarting\n"
                 "        arno:      start with <minBas> Arnoldi steps\n"
                 "        initShift: if arno=0, start with fixed shift <iniShift>\n");
          exit(0); }

  int rank, num_proc;
  int i, ierr;
  int verbose;

  comm_ptr_t comm;
  op_ptr_t A_op; // this is a wrapper for the CRS matrix which we pass to the actual solver
  op_ptr_t B_op=NULL; // no mass matrix up to now
  
  const_map_ptr_t map; // map (element distribution) of vectors according to 
                       // the distribution of matrix rows
  mvec_ptr_t X; // multivector for getting the eigenvectors
  
  ST* evals; // for real non-symmetric matrices we can get complex pairs,
             // so we need twice the amount of memory to store the eigenvalues 
  MT* resid;
  
  int* is_cmplx=NULL; // only required for the real case for indicating complex EV

  phist_jadaOpts_t opts;

  
  int num_eigs,num_iters;
  
  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&ierr),ierr);

  PHIST_ICHK_IERR(phist_comm_create(&comm,&ierr),ierr);

  PHIST_ICHK_IERR(phist_comm_get_rank(comm, &rank,&ierr),ierr);
  PHIST_ICHK_IERR(phist_comm_get_size(comm, &num_proc,&ierr),ierr);

  verbose= (rank==0);

  int nSpins = atoi(argv[1]);
 
  phist_jadaOpts_setDefaults(&opts);
  
  opts.which = SR;
  if (argc>=3)
  {
    opts.numEigs=atoi(argv[2]);
  }

  if (argc>=4)
  {
    opts.convTol=(MT)atof(argv[3]);
  }

  if (argc>=5)
  {
    opts.maxIters=atoi(argv[4]);
  }

  if (argc>=6)
  {
    opts.minBas=atoi(argv[5]);
  }

  if (argc<7)
  {
    opts.maxBas=opts.minBas+20;
  }
  else
  {
    opts.maxBas=atoi(argv[6]);
  }

  if (argc>=8)
  {
    opts.arno=atoi(argv[7]);
  }

  if (argc>=9)
  {
    opts.initialShift=atof(argv[8]);
  }

  num_eigs = opts.numEigs;

  // setup matrix
#ifdef PHIST_KERNEL_LIB_GHOST
  ghost_context_t * ctx = NULL;
  ghost_mat_t * mat = NULL;
#else
  TYPE(crsMat_ptr) mat = NULL;
#endif

  ghost_midx_t DIM;
  ghost_midx_t conf_spinZ[3] = {nSpins,nSpins/2,0};
  SpinChainSZ( -2, &DIM, conf_spinZ, NULL);


  matfuncs_info_t info;
  SpinChainSZ( -1, NULL, NULL, &info);
  if ( my_datatype != info.datatype)
  {
    printf("error: datatyte does not match\n");
    exit(0);
  }
#ifdef PHIST_KERNEL_LIB_FORTRAN
  PHIST_ICHK_IERR(SUBR(crsMat_create_fromRowFunc)(&mat,
        info.nrows, info.ncols, info.row_nnz,
        &SpinChainSZ, &ierr), ierr);
#endif

#ifdef PHIST_KERNEL_LIB_GHOST
  ghost_error_t err=ghost_createContext(&ctx, info.nrows ,
      info.ncols,GHOST_CONTEXT_DEFAULT,NULL,MPI_COMM_WORLD,1.);
  if (err!=GHOST_SUCCESS)
  {
    PHIST_OUT(PHIST_ERROR,"error returned from createContext (file %s, line %d)",__FILE__,__LINE__);
  }

  ghost_mtraits_t mtraits = GHOST_MTRAITS_INIT(.format = GHOST_SPM_FORMAT_CRS, .flags = GHOST_SPM_DEFAULT, .datatype = my_datatype);
  mat = ghost_createMatrix(ctx,&mtraits,1);
  mat->fromRowFunc( mat, info.row_nnz , 0, &SpinChainSZ, 0);
  ghost_printMatrixInfo(mat);
#endif


  PHIST_ICHK_IERR(SUBR(crsMat_get_domain_map)(mat, &map,&ierr),ierr);

  PHIST_ICHK_IERR(SUBR(mvec_create)(&X,map,num_eigs+1,&ierr),ierr);
//  PHIST_ICHK_IERR(SUBR(mvec_random)(X,&ierr),ierr);
  // start with constant vector to make runs reproducible
  PHIST_ICHK_IERR(SUBR(mvec_put_value)(X,ONE,&ierr),ierr);

  PHIST_ICHK_IERR(SUBR(mvec_view_block)(X,&opts.v0,0,0,&ierr),ierr);
  
  // create operator wrapper for computing Y=A*X using a CRS matrix
  A_op = (op_ptr_t)malloc(sizeof(TYPE(op)));
  PHIST_ICHK_IERR(SUBR(op_wrap_crsMat)(A_op,mat,&ierr),ierr);

  
  
  // allocate memory for eigenvalues and residuals. We allocate
  // one extra entry because in the real case we may get that the
  // last EV to converge is a complex pair (requirement of JDQR)
  evals = (ST*)malloc((num_eigs+1)*sizeof(ST));
  resid = (MT*)malloc((num_eigs+1)*sizeof(MT));
  is_cmplx = (int*)malloc((num_eigs+1)*sizeof(int));

  // first column in X is currently used as starting vector of Arnoldi in jdqr. The first 
  // jmin vectors are constructed by an Arnoldi process for stability reasons.
  int nloc,lda; 
  ST* valX0;
  MT nrmX0[num_eigs+1];
  PHIST_ICHK_IERR(SUBR(mvec_my_length)(X,&nloc,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_extract_view)(X,&valX0,&lda,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_normalize)(X,nrmX0,&ierr),ierr);

  SUBR(jdqr)(A_op,B_op,X,evals,resid,is_cmplx, 
        opts,
        &num_eigs,&num_iters,
        &ierr);

  if (ierr!=0)
    {
    if (verbose) fprintf(stdout,"code %d returned from jdqr\n",ierr);
    if (ierr<0) return ierr;
    }
  if (verbose)
    {
    fprintf(stdout,"Found %d eigenpairs after %d iterations\n",num_eigs,num_iters);
    }

  MT expRes[MIN(num_eigs,1)];

  if (num_eigs>0)
    {
    // compute residuals explicitly
    TYPE(mvec_ptr) R=NULL,Xv=NULL;
    PHIST_ICHK_IERR(SUBR(mvec_create)(&R,map,num_eigs,&ierr),ierr);
    PHIST_ICHK_IERR(SUBR(mvec_view_block)(X,&Xv,0,num_eigs-1,&ierr),ierr);
    PHIST_ICHK_IERR(SUBR(crsMat_times_mvec)(ONE,mat,Xv,ZERO,R,&ierr),ierr);
#ifdef IS_COMPLEX
    PHIST_ICHK_IERR(SUBR(mvec_vadd_mvec)(evals,Xv,-ONE,R,&ierr),ierr);
#else
    // we have complex pairs as [v_r, v_i] and [lambda_r, lambda_i] right now.
    // To get the residual correct, create the block diagonal matrix D with   
    // D_j=[lambda_r, lambda_i; -lambda_i, lambda_r] for complex pairs and    
    // then compute A*X-X*D as the residual
    TYPE(sdMat_ptr) D=NULL;
    PHIST_ICHK_IERR(SUBR(sdMat_create)(&D,num_eigs,num_eigs,comm,&ierr),ierr);
    PHIST_ICHK_IERR(SUBR(sdMat_put_value)(D,ZERO,&ierr),ierr);
    ST *D_raw=NULL;
    lidx_t ldD;
    PHIST_ICHK_IERR(SUBR(sdMat_extract_view)(D,&D_raw,&ldD,&ierr),ierr);
    i=0;
    while (i<num_eigs)
      {
      D_raw[i*ldD+i]= evals[i];
      if (is_cmplx[i])
        {
        D_raw[(i+1)*ldD+(i+1)]=evals[i];
        D_raw[i*ldD+(i+1)]=-evals[i+1];
        D_raw[(i+1)*ldD+i]= evals[i+1];
        i++;
        }
      i++;
      }
    PHIST_ICHK_IERR(SUBR(mvec_times_sdMat)(ONE,Xv,D,-ONE,R,&ierr),ierr);
    PHIST_ICHK_IERR(SUBR(sdMat_delete)(D,&ierr),ierr);
#endif
    PHIST_ICHK_IERR(SUBR(mvec_norm2)(R,expRes,&ierr),ierr);
    PHIST_ICHK_IERR(SUBR(mvec_delete)(Xv,&ierr),ierr);
    PHIST_ICHK_IERR(SUBR(mvec_delete)(R,&ierr),ierr);
    }

  if (verbose && num_eigs>0)
    {
    fprintf(stdout,"  Eigenvalue\t\t\t\t\t\tRitz Residual\tExpl. Residual\n");
    fprintf(stdout,"======================================================================================\n");
    int i=0;
    while (i<num_eigs)
      {
#ifdef IS_COMPLEX
      fprintf(stdout,"%24.16e%+24.16ei\t%3.1e\t\t%3.1e\n",
        REAL(evals[i]),IMAG(evals[i]),resid[i],expRes[i]);
#else
      if (is_cmplx[i])
        {
        fprintf(stdout,"%24.16e%+24.16ei\t%3.1e\t\t%3.1e\n",evals[i],evals[i+1],resid[i],expRes[i]);
        fprintf(stdout,"%24.16e%+24.16ei\t%3.1e\t\t%3.1e\n",evals[i],-evals[i+1],resid[i+1],expRes[i+1]);
        i++;
        }
      else
        {
        fprintf(stdout,"%24.16e\t\t\t\t%3.1e\t\t%3.1e\n",evals[i],resid[i],expRes[i]);
        }
#endif
      i++;
      }
    }

  free(evals);
  free(resid);
  free(is_cmplx);
#ifdef PHIST_KERNEL_LIB_GHOST
  mat->destroy(mat);
  ghost_freeContext(ctx);
#else
  PHIST_ICHK_IERR(SUBR(crsMat_delete)(mat,&ierr),ierr);
#endif
  PHIST_ICHK_IERR(SUBR(mvec_delete)(X,&ierr),ierr);
  free(A_op);
  PHIST_ICHK_IERR(phist_kernels_finalize(&ierr),ierr);
  return ierr;
  }
