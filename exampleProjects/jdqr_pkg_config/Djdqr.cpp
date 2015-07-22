#include "phist_config.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "phist_macros.h"
#include "phist_enums.h"
#include "phist_kernels.h"
#include "phist_operator.h"
#include "phist_jdqr.h"
#include "phist_jadaOpts.h"
#include "phist_gen_d.h"
#include "phist_driver_utils.h"
#include "phist_ScalarTraits.hpp"
#include "phist_std_typedefs.hpp"

// JDQR reads an input matrix and tries to compute a given number of            
// exterior eigenpairs at either end of the spectrum. Input parameters          
// are (in this order on the command line)                                      
//                                                                              
// <matrix file name> :  a .mm or .bin file (matrix market or ghost binary CRS) 
// <num eigs>         : number of desired eigenpairs [10]                       
// <which>            : can be "LM", "SM", "LR", "SR" or "TARGET" (Largest/     
//                      Smallest Magnitude/Real part) ["LM"]                    
//                      If "TARGET" is chosen, use "shift" to set the target.   
//                      In this case, JDQR assumes that interior eigenvalues are
//                      sought and you should use CARP-CG as "innerSolver".     
// <tol>              : convergence tolerance [1.0e-10]                         
// <max iters>        : maximum number of iterations allowed [250]              
// <minBas>           : number of vectors to keep in the basis upon restart [10]
// <maxBas>           : number of vectors to generate before restarting         
//                      [minBas+20]                                             
// <arno>             : if 0, starts with an initial shift. If 1, starts with   
//                      <minBas> iterations of Arnoldi [1 unless which=TARGET]  
// <initial shift>    : given initial shift (particularly interesting if you    
//                      set <which> to "TARGET", in this case you should set    
//                      arno=0 as well).
int main(int argc, char** argv)
{
  int iflag=0;
  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&iflag),iflag);

  PHIST_MAIN_TASK_BEGIN
  
  int i;
  int rank, num_proc;
  int verbose;

  comm_ptr_t comm;
  sparseMat_ptr_t A;
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

  char* filename;
  
  int num_eigs,num_iters;
  
  PHIST_ICHK_IERR(phist_comm_create(&comm,&iflag),iflag);

  PHIST_ICHK_IERR(phist_comm_get_rank(comm, &rank,&iflag),iflag);
  PHIST_ICHK_IERR(phist_comm_get_size(comm, &num_proc,&iflag),iflag);

  verbose= (rank==0);

  if (argc<2)
  {
    if (verbose) fprintf(stdout,"Usage: %s <matrix> <num eigs> <which> <tol> <max it>\n"
                                "               <min bas> <max bas> <arno> <init shift> <solver type>\n\n"
                                "       \"matrix\" is the matrix file to be read/created (compulsory)\n\n"
                                "       all other args are optional but are interpreted in the order given, so only the last\n"
                                "       few can be left out.\n\n"
                                "TODO - documentation of options. For now, look in the driver file %s for more details.\n",
                                argv[0],__FILE__);
        SUBR(create_matrix)(NULL,NULL,"usage",&iflag);
    return 1;
  }

  filename = argv[1];
 
  phist_jadaOpts_setDefaults(&opts);
  
  if (argc>=3)
  {
    opts.numEigs=atoi(argv[2]);
  }

  if (argc>=4)
  {
    if (!strcmp(argv[3],"LM")) opts.which=LM;
    else if (!strcmp(argv[3],"SM")) opts.which=SM;
    else if (!strcmp(argv[3],"LR")) opts.which=LR;
    else if (!strcmp(argv[3],"SR")) opts.which=SR;
    else if (!strcmp(argv[3],"TARGET")) opts.which=TARGET;
    else
    {
      if (verbose) 
      {
        fprintf(stderr,"ERROR: the which parameter may be \"LM\",\"SM\",\"LR\" or \"SR\".\n");
        fprintf(stderr,"       run the program without args to get a usage message.\n");
      }
      return -1;
    }
  }
    
  if (argc>=5)
  {
    opts.convTol=(MT)atof(argv[4]);
  }

  if (argc>=6)
  {
    opts.maxIters=atoi(argv[5]);
  }

  if (argc>=7)
  {
    opts.minBas=atoi(argv[6]);
  }

  if (argc<8)
  {
    opts.maxBas=opts.minBas+20;
  }
  else
  {
    opts.maxBas=atoi(argv[7]);
  }

  if (argc>=9)
  {
    opts.arno=atoi(argv[8]);
  }

  if (argc>=10)
  {
    opts.initialShift=atof(argv[9]);
  }
  
  if (argc>=11)
  {
    if (!strcmp(argv[10],"GMRES")) opts.innerSolvType=GMRES;
    else if (!strcmp(argv[10],"CARP_CG")) opts.innerSolvType=CARP_CG;
    else
    {
      if (verbose) 
      {
        fprintf(stderr,"ERROR: the linSolv parameter may be \"GMRES\" or \"CARP_CG\" right now\n");
        fprintf(stderr,"       run the program without args to get a usage message.\n");
      }
      return -1;
    }
  }
  else
  {
    if (opts.which==TARGET) 
    {
      opts.innerSolvType=CARP_CG;
    }
    else
    {
      opts.innerSolvType=GMRES;
    }
  }

  num_eigs = opts.numEigs;

  iflag = PHIST_SPARSEMAT_REPARTITION;
  PHIST_ICHK_IERR(SUBR(create_matrix)(&A,comm,filename,&iflag),iflag);
  
  PHIST_ICHK_IERR(SUBR(sparseMat_get_domain_map)(A, &map,&iflag),iflag);

  PHIST_ICHK_IERR(SUBR(mvec_create)(&X,map,num_eigs+1,&iflag),iflag);
//  PHIST_ICHK_IERR(SUBR(mvec_random)(X,&iflag),iflag);
  // start with constant vector to make runs reproducible
  PHIST_ICHK_IERR(SUBR(mvec_put_value)(X,st::one(),&iflag),iflag);

  PHIST_ICHK_IERR(SUBR(mvec_view_block)(X,&opts.v0,0,0,&iflag),iflag);
  
  // create operator wrapper for computing Y=A*X using a CRS matrix
  A_op = (op_ptr_t)malloc(sizeof(TYPE(op)));
  PHIST_ICHK_IERR(SUBR(op_wrap_sparseMat)(A_op,A,&iflag),iflag);

  
  
  // allocate memory for eigenvalues and residuals. We allocate
  // one extra entry because in the real case we may get that the
  // last EV to converge is a complex pair (requirement of JDQR)
  evals = (ST*)malloc((num_eigs+1)*sizeof(ST));
  resid = (MT*)malloc((num_eigs+1)*sizeof(MT));
  is_cmplx = (int*)malloc((num_eigs+1)*sizeof(int));

  // first column in X is currently used as starting vector of Arnoldi in jdqr. The first 
  // jmin vectors are constructed by an Arnoldi process for stability reasons.
  lidx_t nloc,lda; 
  ST* valX0;
  MT nrmX0[num_eigs+1];
  PHIST_ICHK_IERR(SUBR(mvec_my_length)(X,&nloc,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_extract_view)(X,&valX0,&lda,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_normalize)(X,nrmX0,&iflag),iflag);

  SUBR(jdqr)(A_op,B_op,X,evals,resid,is_cmplx, 
        opts,
        &num_eigs,&num_iters,
        &iflag);

  if (iflag!=0)
    {
    if (verbose) fprintf(stdout,"code %d returned from jdqr\n",iflag);
    if (iflag<0) return iflag;
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
    PHIST_ICHK_IERR(SUBR(mvec_create)(&R,map,num_eigs,&iflag),iflag);
    PHIST_ICHK_IERR(SUBR(mvec_view_block)(X,&Xv,0,num_eigs-1,&iflag),iflag);
    PHIST_ICHK_IERR(SUBR(sparseMat_times_mvec)(st::one(),A,Xv,st::zero(),R,&iflag),iflag);
#ifdef IS_COMPLEX
    PHIST_ICHK_IERR(SUBR(mvec_vadd_mvec)(evals,Xv,-st::one(),R,&iflag),iflag);
#else
    // we have complex pairs as [v_r, v_i] and [lambda_r, lambda_i] right now.
    // To get the residual correct, create the block diagonal matrix D with   
    // D_j=[lambda_r, lambda_i; -lambda_i, lambda_r] for complex pairs and    
    // then compute A*X-X*D as the residual
    TYPE(sdMat_ptr) D=NULL;
    PHIST_ICHK_IERR(SUBR(sdMat_create)(&D,num_eigs,num_eigs,comm,&iflag),iflag);
    PHIST_ICHK_IERR(SUBR(sdMat_put_value)(D,st::zero(),&iflag),iflag);
    ST *D_raw=NULL;
    lidx_t ldD;
    PHIST_ICHK_IERR(SUBR(sdMat_extract_view)(D,&D_raw,&ldD,&iflag),iflag);
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
    PHIST_ICHK_IERR(SUBR(mvec_times_sdMat)(st::one(),Xv,D,-st::one(),R,&iflag),iflag);
    PHIST_ICHK_IERR(SUBR(sdMat_delete)(D,&iflag),iflag);
#endif
    PHIST_ICHK_IERR(SUBR(mvec_norm2)(R,expRes,&iflag),iflag);
    PHIST_ICHK_IERR(SUBR(mvec_delete)(Xv,&iflag),iflag);
    PHIST_ICHK_IERR(SUBR(mvec_delete)(R,&iflag),iflag);
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
        st::real(evals[i]),st::imag(evals[i]),resid[i],expRes[i]);
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
  
  PHIST_ICHK_IERR(SUBR(sparseMat_delete)(A,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_delete)(X,&iflag),iflag);
  free(A_op);

  PHIST_MAIN_TASK_END

  PHIST_ICHK_IERR(phist_kernels_finalize(&iflag),iflag);
  return iflag;
}