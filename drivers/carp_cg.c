#include "phist_config.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "phist_macros.h"
#include "phist_enums.h"
#include "phist_kernels.h"
#include "phist_operator.h"
#include "phist_feastCorrectionSolver.h"

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

// This driver routine creates the Graphene Hamiltonian,
// reads in a number of complex shifts from a textfile, 
// and attempts to solve some linear systems (sI-A)x=b  
// by calling carp_cg via the feastCorrectionSolver in- 
// terface. The idea is to test the imolementation for  
// later use in FEAST.                                  
int main(int argc, char** argv)
  {
     if ( argc < 2 )
     {
          PHIST_SOUT(PHIST_ERROR,"Usage: carp_cg <problem> <shiftFile> <num vecs>\n"
                                 "               <block size> <tol> <maxIter>\n"
                                 "  where: <problem> can either be a filename or\n"
                                 "           \"graphene<L>\" or \"anderson<L>\"\n"
                                 "         <shiftFile>: see \"exampleRuns/graphene_shifts.txt\" for an example.\n"
                                 "                      If omitted, no shift is used.\n"
                                 "         etc. all other args get default values if omitted\n");
          exit(0);
     }

  int rank, num_proc;
  int i, ierr;
  int verbose;

  comm_ptr_t comm;  
  const_map_ptr_t map; // map (element distribution) of vectors according to 
                       // the distribution of matrix rows
  mvec_ptr_t B, *X_r, *X_i; // multivectors that will hold the RHS and the 
                            //complex solutions
  TYPE(mvec_ptr) X_r_ex0, X_i_ex0;// we construct B such that (sigma[0]*I-A)X_ex0=B,
                                  // so that we can check the error for at least one
                                  // of the systems.
  MT *resid, *err0;
  int nshifts,nrhs,blockSize,maxIter;
  MT tol;
  MT *sigma_r, *sigma_i;
  
  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&ierr),ierr);

  PHIST_ICHK_IERR(phist_comm_create(&comm,&ierr),ierr);

  PHIST_ICHK_IERR(phist_comm_get_rank(comm, &rank,&ierr),ierr);
  PHIST_ICHK_IERR(phist_comm_get_size(comm, &num_proc,&ierr),ierr);

  verbose= (rank==0);
  int iarg=1;
  const char* problem = argv[iarg++];
  
  char* shiftFileName=NULL;
  if (argc>iarg)
  {
    shiftFileName=argv[iarg++];
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
      if (sigma_i[i]==(MT)0.0)
      {
        PHIST_SOUT(PHIST_ERROR,"all shifts must have a non-zero imaginary part\n");
        exit(-1);
      }
    }
  
    fclose(shiftFile);
  }
    
  // setup matrix
  TYPE(crsMat_ptr) mat = NULL;
  
  // this is in the tools/driver_utils.h header, a useful tool for
  // generating our favorite test matrices or reading them from a file:
  PHIST_ICHK_IERR(SUBR(create_matrix)(&mat,problem,&ierr),ierr);
  
  PHIST_ICHK_IERR(SUBR(crsMat_get_domain_map)(mat, &map,&ierr),ierr);

  PHIST_ICHK_IERR(SUBR(mvec_create)(&B,map,nrhs,&ierr),ierr);
  
  // vectors to hold exact solution for system 0
  PHIST_ICHK_IERR(SUBR(mvec_create)(&X_r_ex0,map,nrhs,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_create)(&X_i_ex0,map,nrhs,&ierr),ierr);
  
  X_r = (TYPE(mvec_ptr)*)malloc(nshifts*sizeof(TYPE(mvec_ptr)));
  X_i = (TYPE(mvec_ptr)*)malloc(nshifts*sizeof(TYPE(mvec_ptr)));
  
  for (int i=0;i<nshifts;i++)
  {
    PHIST_ICHK_IERR(SUBR(mvec_create)(&X_r[i],map,nrhs,&ierr),ierr);
    PHIST_ICHK_IERR(SUBR(mvec_create)(&X_i[i],map,nrhs,&ierr),ierr);
    // start with zero vector to make runs reproducible
    PHIST_ICHK_IERR(SUBR(mvec_put_value)(X_r[i],ZERO,&ierr),ierr);
    PHIST_ICHK_IERR(SUBR(mvec_put_value)(X_i[i],ZERO,&ierr),ierr);
  }

/////////////////////////////////////////////////////////////
// generate linear systems                                 //
/////////////////////////////////////////////////////////////

// we assume that all shifts have a non-zero imaginary part.
// The complex linear systems j=0..nshifts-1 are given by
// (sigma[j]*I-A)(X_r[j]+i*X_i[j])=B+0i (i=sqrt(-1)).
//
// The equivalent real formulation is
//
// | sigma_r[j]I-A      -sigma_i[j]I  | |X_r[j]|   |B|
// | sigma_i[j]I        sigma_r[j]I-A | |X_i[j]| = |0|
//
// We first construct an exact solution for system 0:
// X_i_ex0 random, X_r_ex0 to satisfy the second block row,
// B from first row given X_r_ex0 and X_i_ex0.
  PHIST_ICHK_IERR(SUBR(mvec_random)(X_i_ex0,&ierr),ierr);
  
  // x_r = 1/sig_i(A-sig_r I)x_i
  PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(-sigma_r[0]/sigma_i[0],X_i_ex0,0.0,X_r_ex0,&ierr),ierr);
  if (sigma_i[0]!=(MT)0.0)
  {
    PHIST_ICHK_IERR(SUBR(crsMat_times_mvec)(1.0/sigma_i[0],mat,X_i_ex0,1.0,X_r_ex0,&ierr),ierr);
  }
  else
  {
    TYPE(mvec_ptr) tmp = X_i_ex0;
    X_i_ex0=X_r_ex0;
    X_r_ex0=tmp;
    PHIST_ICHK_IERR(SUBR(mvec_put_value)(X_i_ex0,0.0,&ierr),ierr);
  }
  
  // compute rhs B to match this exact solution for sigma[0]:
  PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(sigma_r[0],X_r_ex0,0.0,B,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(crsMat_times_mvec)(-1.0,mat,X_r_ex0,1.0,B,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(-sigma_i[0],X_i_ex0,1.0,B,&ierr),ierr);

/*
TYPE(sdMat_ptr) Rtmp=NULL;
PHIST_ICHK_IERR(SUBR(sdMat_create)(&Rtmp,nrhs,nrhs,NULL,&ierr),ierr);
PHIST_ICHK_IERR(SUBR(mvec_random)(B,&ierr),ierr);
PHIST_ICHK_IERR(SUBR(mvec_QR)(B,Rtmp,&ierr),ierr);
PHIST_ICHK_IERR(SUBR(sdMat_delete)(Rtmp,&ierr),ierr);
*/

///////////////////////////////////////////////////////////////////
// setup and solve via feastCorrectionSolver                     //
///////////////////////////////////////////////////////////////////
  TYPE(feastCorrectionSolver_ptr) fCorrSolver;
  PHIST_ICHK_IERR(SUBR(feastCorrectionSolver_create)
        (&fCorrSolver, mat, CARP_CG, blockSize, nshifts,sigma_r, sigma_i, &ierr),ierr);

  int numBlocks=nrhs/blockSize;

  PHIST_ICHK_IERR(SUBR(feastCorrectionSolver_run)
        (fCorrSolver, B, tol, maxIter,X_r, X_i, &ierr),ierr);

///////////////////////////////////////////////////////////////////
// compute residuals and error                                   //
///////////////////////////////////////////////////////////////////

PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(-ONE,X_r[0],ONE,X_r_ex0,&ierr),ierr);
PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(-ONE,X_i[0],ONE,X_i_ex0,&ierr),ierr);
double nrm_err0_1[nrhs], nrm_err0_2[nrhs];
PHIST_ICHK_IERR(SUBR(mvec_dot_mvec)(X_i_ex0,X_i_ex0,nrm_err0_1,&ierr),ierr);
PHIST_ICHK_IERR(SUBR(mvec_dot_mvec)(X_r_ex0,X_r_ex0,nrm_err0_2,&ierr),ierr);

PHIST_SOUT(PHIST_VERBOSE,"for first shift, first rhs:\n");
PHIST_SOUT(PHIST_VERBOSE,"err in re(x): %e\n",SQRT(nrm_err0_2[0]));
PHIST_SOUT(PHIST_VERBOSE,"       im(x): %e\n",SQRT(nrm_err0_1[0]));

///////////////////////////////////////////////////////////////////
// clean up afterwards                                           //
///////////////////////////////////////////////////////////////////

  PHIST_ICHK_IERR(SUBR(feastCorrectionSolver_delete)(fCorrSolver,&ierr),ierr);

  free(sigma_r);
  free(sigma_i);

  PHIST_ICHK_IERR(SUBR(crsMat_delete)(mat,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_delete)(B,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_delete)(X_r_ex0,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_delete)(X_i_ex0,&ierr),ierr);
  for (int i=0;i<nshifts;i++)
  {
    PHIST_ICHK_IERR(SUBR(mvec_delete)(X_r[i],&ierr),ierr);
    PHIST_ICHK_IERR(SUBR(mvec_delete)(X_i[i],&ierr),ierr);
  }
  PHIST_ICHK_IERR(phist_kernels_finalize(&ierr),ierr);
  return ierr;
  }
