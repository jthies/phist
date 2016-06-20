#include "phist_config.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#ifdef PHIST_HAVE_OPENMP
#include <omp.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <stdbool.h>
#include <math.h>

#include "phist_macros.h"
#include "phist_enums.h"
#include "phist_kernels.h"
#include "phist_operator.h"
#include "phist_carp_cg.h"

// ghost/physics
#ifdef PHIST_KERNEL_LIB_GHOST
#include "ghost.h"
GHOST_REGISTER_DT_D(my_datatype)
#endif


#include "phist_gen_d.h"
#include "phist_driver_utils.h"

typedef enum {
FROM_FILE=0,
GRAPHENE=1,
ANDERSON=2
} problem_t;

#define TEST_SYSTEM

int test_rhs(ghost_gidx i, ghost_lidx j, void * val, void * aux)
{
  double* dval=(double*)val;
  *dval = (double)(i+1)/(double)(j+1);
  return 0;
}
/*
int test_xex_i(ghost_gidx i, ghost_lidx j, void * val, void * aux)
{
  double* dval=(double*)val;
  *dval = (double)(j+1)/(double)(i+1);
  return 0;
}
*/

// This driver routine creates the Graphene Hamiltonian,
// reads in a number of complex shifts from a textfile, 
// and attempts to solve some linear systems (sI-A)x=b  
// by calling carp_cg via the feastCorrectionSolver in- 
// terface. The idea is to test the imolementation for  
// later use in FEAST.                                  
int main(int argc, char** argv)
  {
     int iflag;
     if ( argc < 2 )
     {
       PHIST_SOUT(PHIST_ERROR,"Usage: carp_cg    <problem> <shiftFile> <num vecs>\n"
                                 "               <block size> <tol> <maxIter> <omega>\n"
                                 "  where: <problem> describes the input matrix (see below)\n"
                                 "         <shiftFile>: see \"exampleRuns/graphene_shifts.txt\" for an example.\n"
                                 "                      If omitted or 'none', no shift is used.\n"
                                 "         etc. all other args get default values if omitted\n");
       // print usage message for creating/reading a matrix
       SUBR(create_matrix)(NULL, NULL, "usage",&iflag);
       exit(1);
     }

  int rank, num_proc;
  int i;

  phist_comm_ptr comm;  
  phist_const_map_ptr map; // map (element distribution) of vectors according to 
                       // the distribution of matrix rows
  mvec_ptr B, *X_r, *X_i; // multivectors that will hold the RHS and the 
                            //complex solutions
  TYPE(mvec_ptr) X_r_ex0, X_i_ex0;// we construct B such that (sigma[0]*I-A)X_ex0=B,
                                  // so that we can check the error for at least one
                                  // of the systems.
  TYPE(mvec_ptr) err_r, err_i;  // for computing errors/residuals
  
  MT *resid, *err0;
  int nshifts,nrhs,blockSize,maxIter;
  MT tol;
  MT *sigma_r, *sigma_i;
  MT omega; /* relaxation parameter */
  
  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&iflag),iflag);

  PHIST_MAIN_TASK_BEGIN

  PHIST_ICHK_IERR(phist_comm_create(&comm,&iflag),iflag);

  PHIST_ICHK_IERR(phist_comm_get_rank(comm, &rank,&iflag),iflag);
  PHIST_ICHK_IERR(phist_comm_get_size(comm, &num_proc,&iflag),iflag);

  //int verbose= (rank==0);
  int iarg=1;
  const char* problem = argv[iarg++];
  
  char* shiftFileName=NULL;
  if (argc>iarg)
  {
    shiftFileName=argv[iarg++];
    if (!strncasecmp(shiftFileName,"none",4)) shiftFileName=NULL;
  }
  if (argc>iarg)
  {
    nrhs=atoi(argv[iarg++]);
  }
  else
  {
    nrhs=1;
  }
  
  // walk through array of vectors in blocks of <blockSize>,
  blockSize=nrhs;

  if (argc>iarg)
  {
    blockSize=atoi(argv[iarg++]);
  }

  blockSize=MIN(nrhs,blockSize);
 
  PHIST_SOUT(PHIST_VERBOSE,"Solve with %d rhs/shift and block size %d\n",
        nrhs,blockSize);
  
  if (argc>iarg)
  {
    tol=atof(argv[iarg++]);
  }
  else
  {
    tol=(MT)1.0e-12;
  }
  if (argc>iarg)
  {
    maxIter=atoi(argv[iarg++]);
  }
  else
  {
    maxIter=1000;
  }
  
  if (argc>iarg)
  {
    omega=atof(argv[iarg++]);
  }
  else
  {
    omega=1.0;
  }

    int num_complex=0;
    
  if (shiftFileName==NULL)
  {
    nshifts=1;
    sigma_r=(MT*)malloc(1*sizeof(MT));
    sigma_i=(MT*)malloc(1*sizeof(MT));
    *sigma_r=0.0;
    *sigma_i=0.0;
  }
  else
  {
    FILE *shiftFile=fopen(shiftFileName,"r");
  
    if (!shiftFile)
    {
      PHIST_SOUT(PHIST_ERROR,"could not open %s for reading the shifts\n",shiftFileName);
      exit(-1);
    }

    fscanf(shiftFile,"%d",&nshifts);
  
    if (nshifts<=0)
    {
      PHIST_SOUT(PHIST_ERROR,"first line of shift file should be nshifts>0 [int]\n"
                           "found: %d in file %s\n",nshifts,shiftFileName);
      exit(-1);
    }
  
    PHIST_SOUT(PHIST_VERBOSE,"using %d shifts:\n",nshifts);
    PHIST_SOUT(PHIST_VERBOSE,"================\n");
  
    sigma_r = (MT*)malloc(nshifts*sizeof(MT));
    sigma_i = (MT*)malloc(nshifts*sizeof(MT));
  
    for (i=0;i<nshifts;i++)
    {
      fscanf(shiftFile,"%lf %lf",&sigma_r[i],&sigma_i[i]);
      PHIST_SOUT(PHIST_VERBOSE,"sigma[%d]= %g + %gi\n",i,sigma_r[i],sigma_i[i]);
      if (sigma_i[i]!=(MT)0.0)
      {
        num_complex++;
      }
    }
  
    fclose(shiftFile);
  }
  
  if (!(num_complex==0||num_complex==nshifts))
  {
    PHIST_SOUT(PHIST_ERROR,
    "the way we construct the real-valued rhs\n"
    "to the complex linear system (sI-A)X=B below\n"
    "requires all shifts to have a non-zero imaginary\n"
    "part. If all shifts are real, we simply use a real\n"
    "valued exact solution\n"
    "(file %s, line %d)\n",__FILE__,__LINE__);
    return -1;
  }
  
  // setup matrix
  TYPE(sparseMat_ptr) mat = NULL;
  
  // this is in the tools/driver_utils.h header, a useful tool for
  // generating our favorite test matrices or reading them from a file:
       iflag=PHIST_SPARSEMAT_REPARTITION|PHIST_SPARSEMAT_OPT_CARP;
  PHIST_ICHK_IERR(SUBR(create_matrix)(&mat,comm,problem,&iflag),iflag);
  
  PHIST_ICHK_IERR(SUBR(sparseMat_get_domain_map)(mat, &map,&iflag),iflag);

  PHIST_ICHK_IERR(SUBR(mvec_create)(&B,map,nrhs,&iflag),iflag);
  
  // vectors to hold exact solution for system 0
  PHIST_ICHK_IERR(SUBR(mvec_create)(&X_r_ex0,map,nrhs,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_create)(&X_i_ex0,map,nrhs,&iflag),iflag);

  PHIST_ICHK_IERR(SUBR(mvec_create)(&err_r,map,nrhs,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_create)(&err_i,map,nrhs,&iflag),iflag);
  
  X_r = (TYPE(mvec_ptr)*)malloc(nshifts*sizeof(TYPE(mvec_ptr)));
  X_i = (TYPE(mvec_ptr)*)malloc(nshifts*sizeof(TYPE(mvec_ptr)));
  
  for (int i=0;i<nshifts;i++)
  {
    PHIST_ICHK_IERR(SUBR(mvec_create)(&X_r[i],map,nrhs,&iflag),iflag);
    PHIST_ICHK_IERR(SUBR(mvec_create)(&X_i[i],map,nrhs,&iflag),iflag);
    // start with zero vector to make runs reproducible
    PHIST_ICHK_IERR(SUBR(mvec_put_value)(X_r[i],ZERO,&iflag),iflag);
    PHIST_ICHK_IERR(SUBR(mvec_put_value)(X_i[i],ZERO,&iflag),iflag);
  }

/////////////////////////////////////////////////////////////
// generate linear systems                                 //
/////////////////////////////////////////////////////////////

#ifdef TEST_SYSTEM
  PHIST_ICHK_IERR(SUBR(mvec_put_func)(B,&test_rhs,NULL,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_put_value)(X_r_ex0,0.0,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_put_value)(X_i_ex0,0.0,&iflag),iflag);
#else
if (num_complex==0)
{
  PHIST_ICHK_IERR(SUBR(create_sol_and_rhs)(problem,mat,X_r_ex0,B,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_put_value)(X_i_ex0,ZERO,&iflag),iflag);
  // compute rhs B to match this exact solution for sigma[0]:
  PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(-sigma_r[0],X_r_ex0,ONE,B,&iflag),iflag);
}
else
{
// we assume that all shifts have a non-zero imaginary part.
// The complex linear systems j=0..nshifts-1 are given by
// (A-sigma[j]*I)(X_r[j]+i*X_i[j])=B+0i (i=sqrt(-1)).
//
// The equivalent real formulation is
//
// | A-sigma_r[j]I       sigma_i[j]I  | |X_r[j]|   |B|
// | -sigma_i[j]I        A-sigma_r[j]I | |X_i[j]| = |0|
//
// We first construct an exact solution for system 0:
// X_i_ex0 random, X_r_ex0 to satisfy the second block row,
// B from first row given X_r_ex0 and X_i_ex0.
  PHIST_ICHK_IERR(SUBR(mvec_random)(X_i_ex0,&iflag),iflag);
  
  // x_r = 1/sig_i(A-sig_r I)x_i
  PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(-sigma_r[0]/sigma_i[0],X_i_ex0,0.0,X_r_ex0,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(sparseMat_times_mvec)(1.0/sigma_i[0],mat,X_i_ex0,1.0,X_r_ex0,&iflag),iflag);
  
  // compute rhs B to match this exact solution for sigma[0]:
  PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(-sigma_r[0],X_r_ex0,0.0,B,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(sparseMat_times_mvec)(1.0,mat,X_r_ex0,1.0,B,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(sigma_i[0],X_i_ex0,1.0,B,&iflag),iflag);

}
#endif

  // for each shift one by one create the state, iterate and delete the state again
  TYPE(carp_cgState_ptr) carp_cgState;

  for (int i=0; i<nshifts; i++)
  {
    MT tmp_sigma_r[blockSize], tmp_sigma_i[blockSize];
    for (int j=0; j<blockSize; j++) 
    {
      tmp_sigma_r[j]=sigma_r[i];
      tmp_sigma_i[j]=sigma_i[i];
    }
    // at this point the carp_cg state will check if there are complex shifts and switch on the "RC" variant
    // in real arithmetic
    PHIST_ICHK_IERR(SUBR(carp_cgState_create)
            (&carp_cgState,mat,NULL,blockSize,tmp_sigma_r,tmp_sigma_i, &iflag),iflag);
            
    // relaxation factor
    for (int j=0;j<blockSize;j++) carp_cgState->omega_[j]=omega;

    PHIST_ICHK_IERR(SUBR(carp_cgState_reset)(carp_cgState,B,NULL,&iflag),iflag);

    SUBR(carp_cgState_iterate)(carp_cgState,X_r[i],X_i[i],tol,maxIter,false,&iflag);

    if (iflag>0)
    {
    PHIST_SOUT(PHIST_WARNING,"FAILURE: solver returned warning code %d\n",iflag);
    }
    else if (iflag<0)
    {
      PHIST_SOUT(PHIST_WARNING,"FAILURE: solver returned error code %d\n",iflag);
    }
    PHIST_ICHK_IERR(SUBR(carp_cgState_delete)(carp_cgState,&iflag),iflag);
  }
///////////////////////////////////////////////////////////////////
// compute residuals and error                                   //
///////////////////////////////////////////////////////////////////
if (iflag>=0)
{
#ifndef TEST_SYSTEM
PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(ONE,X_r_ex0,ZERO,err_r,&iflag),iflag);
PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(-ONE,X_r[0],ONE,err_r,&iflag),iflag);
PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(ONE,X_i_ex0,ZERO,err_i,&iflag),iflag);
PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(-ONE,X_i[0],ONE,err_i,&iflag),iflag);
double nrm_err0_1[nrhs], nrm_err0_2[nrhs];
PHIST_ICHK_IERR(SUBR(mvec_dot_mvec)(err_i,err_i,nrm_err0_1,&iflag),iflag);
PHIST_ICHK_IERR(SUBR(mvec_dot_mvec)(err_r,err_r,nrm_err0_2,&iflag),iflag);

PHIST_SOUT(PHIST_VERBOSE,"for first shift, first rhs:\n");
PHIST_SOUT(PHIST_VERBOSE,"err in re(x): %e\n",SQRT(nrm_err0_2[0]));
PHIST_SOUT(PHIST_VERBOSE,"       im(x): %e\n",SQRT(nrm_err0_1[0]));
#endif
}

// it may be interesting to look at the residual of the exact solution,
// because e.g. MATPDE3D provides an rhs based on an analytical solution,
// so Ax-b will say something about the discretization error.
if (num_complex==0)
{
  // our B is actually (A-sigma_r*I)X_r_ex0. Compute the residuals for X_r[0] and X_r_ex0.
  // Note that the `exact residual' is only 0 if B was 
  // constructed in such a way, if the benchmark problem provides an analytical solution, the
  // discretization error appears here (e.g. BENCH3D problems).
  TYPE(mvec_ptr) res0=err_r, res_ex0=err_i;
  PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(ONE,B,ZERO,res0,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(ONE,B,ZERO,res_ex0,&iflag),iflag);

  PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(sigma_r[0],X_r_ex0,ONE,res_ex0,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(sigma_r[0],X_r[0],ONE,res0,&iflag),iflag);

  PHIST_ICHK_IERR(SUBR(sparseMat_times_mvec)(-ONE,mat,X_r_ex0,ONE,res_ex0,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(sparseMat_times_mvec)(-ONE,mat,X_r[0],ONE,res0,&iflag),iflag);

  double nrm_res[nrhs],nrm_res_ex[nrhs],nrm_rhs[nrhs];
  PHIST_ICHK_IERR(SUBR(mvec_dot_mvec)(B,B,nrm_rhs,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_dot_mvec)(res0,res0,nrm_res,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_dot_mvec)(res_ex0,res_ex0,nrm_res_ex,&iflag),iflag);
  for (i=0;i<nrhs;i++) 
  {
    nrm_rhs[i]=SQRT(nrm_rhs[i]);
    nrm_res_ex[i]=SQRT(nrm_res_ex[i]);
    nrm_res[i]=SQRT(nrm_res[i]);
  }

  PHIST_SOUT(PHIST_INFO,"residual norms ||r||_2/||b||_2\n");
  for (int i=0; i<nrhs; i++)
  {
    PHIST_SOUT(PHIST_INFO,"%d\t%16.8e (%s)\n",i,nrm_res[i]/nrm_rhs[i],nrm_res[i]<=tol*nrm_rhs[i]?"SUCCESS":"FAILURE");
  }

#if PHIST_OUTLEV>=PHIST_DEBUG
  ST *x_val=NULL, *xex_val=NULL,*r_val=NULL,*rex_val=NULL;
  phist_lidx ldx, ldxex, ldr, ldrex,nloc;
  PHIST_ICHK_IERR(SUBR(mvec_extract_view)(X_r_ex0,&xex_val,&ldxex,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_extract_view)(X_r[0],&x_val,&ldx,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_extract_view)(err_r,&rex_val,&ldrex,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_extract_view)(err_i,&r_val,&ldr,&iflag),iflag);

  PHIST_ICHK_IERR(SUBR(mvec_my_length)(X_r_ex0,&nloc,&iflag),iflag);

  PHIST_SOUT(PHIST_DEBUG,"row\tx\tx_ex\tr\tr_ex\n");
  for (int i=0; i<8; i++)
  {
    if (i>=nloc) break;
    PHIST_SOUT(PHIST_DEBUG,"%d %16.8e %16.8e %16.8e %16.8e\n",
          i, xex_val[i*ldxex],x_val[i*ldx],r_val[i*ldr],rex_val[i*ldrex]);
  }
#endif
}

///////////////////////////////////////////////////////////////////
// clean up afterwards                                           //
///////////////////////////////////////////////////////////////////

  free(sigma_r);
  free(sigma_i);

  PHIST_ICHK_IERR(SUBR(sparseMat_delete)(mat,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_delete)(B,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_delete)(X_r_ex0,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_delete)(X_i_ex0,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_delete)(err_r,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_delete)(err_i,&iflag),iflag);
  for (int i=0;i<nshifts;i++)
  {
    PHIST_ICHK_IERR(SUBR(mvec_delete)(X_r[i],&iflag),iflag);
    PHIST_ICHK_IERR(SUBR(mvec_delete)(X_i[i],&iflag),iflag);
  }

  PHIST_MAIN_TASK_END

  PHIST_ICHK_IERR(phist_kernels_finalize(&iflag),iflag);
  return iflag;
  }
