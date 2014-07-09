#include "phist_config.h"

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
#include "phist_feastCorrectionSolver.h"

// ghost/physics
#ifdef PHIST_KERNEL_LIB_GHOST
#include "ghost.h"
GHOST_REGISTER_DT_D(my_datatype)
#endif
#include "matfuncs.h"


#include "phist_gen_d.h"
#include "phist_driver_utils.h"

// This driver routine creates the Graphene Hamiltonian,
// reads in a number of complex shifts from a textfile, 
// and attempts to solve some linear systems (sI-A)x=b  
// by calling carp_cg via the feastCorrectionSolver in- 
// terface. The idea is to test the imolementation for  
// later use in FEAST.                                  
int main(int argc, char** argv)
  {
     if ( argc < 4 )
     {
          printf("Usage: graphene_carp <nx> <ny> <shiftFile> <num vecs>\n"
                 "                     <tol> <maxIter>\n"
                 "       where: TODO - document/complete args, for now\n"
                 "       please look in the source file %s\n",__FILE__);
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
  
  int nshifts,nvec, maxIter;
  MT tol;
  MT *sigma_r, *sigma_i;
  
  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&ierr),ierr);

  PHIST_ICHK_IERR(phist_comm_create(&comm,&ierr),ierr);

  PHIST_ICHK_IERR(phist_comm_get_rank(comm, &rank,&ierr),ierr);
  PHIST_ICHK_IERR(phist_comm_get_size(comm, &num_proc,&ierr),ierr);

  verbose= (rank==0);

  // problem size
  ghost_idx_t WL[2];
  WL[0] = atoi(argv[1]);
  WL[1] = atoi(argv[2]);
  char* shiftFileName=argv[3];
  
  if (argc>=5)
  {
    nvec=atoi(argv[4]);
  }
  else
  {
    nvec=1;
  }
  
  if (argc>=6)
  {
    tol=atof(argv[5]);
  }
  else
  {
    tol=(MT)1.0e-12;
  }
  if (argc>=7)
  {
    maxIter=atoi(argv[6]);
  }
  else
  {
    maxIter=1000;
  }
  
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
  
  // setup matrix
#ifdef PHIST_KERNEL_LIB_GHOST
  ghost_context_t * ctx = NULL;
  ghost_sparsemat_t * mat = NULL;
#else
  TYPE(crsMat_ptr) mat = NULL;
#endif

  ghost_idx_t DIM = WL[0]*WL[1];
  matfuncs_info_t info;
  crsGraphene( -2, WL, NULL, NULL);
  crsGraphene( -1, NULL, NULL, &info);

#ifdef PHIST_KERNEL_LIB_FORTRAN
  PHIST_ICHK_IERR(SUBR(crsMat_create_fromRowFunc)(&mat,
        info.nrows, info.ncols, info.row_nnz,
        (void(*)(ghost_idx_t,ghost_idx_t*,ghost_idx_t*,void*))&crsGraphene, &ierr), ierr);
#elif defined(PHIST_KERNEL_LIB_GHOST)
  if ( my_datatype != info.datatype)
  {
    printf("error: datatyte does not match\n");
    exit(0);
  }

  ghost_error_t err=ghost_context_create(&ctx, info.nrows ,
      info.ncols,GHOST_CONTEXT_DEFAULT,NULL,0,MPI_COMM_WORLD,1.);
  if (err!=GHOST_SUCCESS)
  {
    PHIST_OUT(PHIST_ERROR,"error returned from createContext (file %s, line %d)",__FILE__,__LINE__);
  }

  ghost_sparsemat_traits_t mtraits = GHOST_SPARSEMAT_TRAITS_INITIALIZER;
  mtraits.datatype = my_datatype;
  ghost_sparsemat_create(&mat,ctx,&mtraits,1);
  ghost_sparsemat_src_rowfunc_t src = GHOST_SPARSEMAT_SRC_ROWFUNC_INITIALIZER;
  src.func = &crsGraphene;

  src.maxrowlen = info.row_nnz;
  mat->fromRowFunc(mat, &src);
  char *str;
  ghost_sparsemat_string(&str,mat);
  printf("%s\n",str);
  free(str); str = NULL;
#else
#error 'this driver is only available for the fortran and ghost kernel libs'
#endif

  PHIST_ICHK_IERR(SUBR(crsMat_get_domain_map)(mat, &map,&ierr),ierr);

  PHIST_ICHK_IERR(SUBR(mvec_create)(&B,map,nvec,&ierr),ierr);
  
  // vectors to hold exact solution for system 0
  PHIST_ICHK_IERR(SUBR(mvec_create)(&X_r_ex0,map,nvec,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_create)(&X_i_ex0,map,nvec,&ierr),ierr);
  
  X_r = (TYPE(mvec_ptr)*)malloc(nshifts*sizeof(TYPE(mvec_ptr)));
  X_i = (TYPE(mvec_ptr)*)malloc(nshifts*sizeof(TYPE(mvec_ptr)));
  
  for (int i=0;i<nshifts;i++)
  {
    PHIST_ICHK_IERR(SUBR(mvec_create)(&X_r[i],map,nvec,&ierr),ierr);
    PHIST_ICHK_IERR(SUBR(mvec_create)(&X_i[i],map,nvec,&ierr),ierr);
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
  PHIST_ICHK_IERR(SUBR(crsMat_times_mvec)(1.0/sigma_i[0],mat,X_i_ex0,1.0,X_r_ex0,&ierr),ierr);
  
  // compute rhs B to match this exact solution for sigma[0]:
  PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(sigma_r[0],X_r_ex0,0.0,B,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(crsMat_times_mvec)(-1.0,mat,X_r_ex0,1.0,B,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(-sigma_i[0],X_i_ex0,1.0,B,&ierr),ierr);

///////////////////////////////////////////////////////////////////
// setup and solve via feastCorrectionSolver                     //
///////////////////////////////////////////////////////////////////
  TYPE(feastCorrectionSolver_ptr) fCorrSolver;
  PHIST_ICHK_IERR(SUBR(feastCorrectionSolver_create)
        (&fCorrSolver, mat, CARP_CG, nvec, nshifts,sigma_r, sigma_i, &ierr),ierr);

  PHIST_ICHK_IERR(SUBR(feastCorrectionSolver_run)
        (fCorrSolver, B, tol, maxIter,X_r, X_i, &ierr),ierr);

///////////////////////////////////////////////////////////////////
// compute residuals and error                                   //
///////////////////////////////////////////////////////////////////

PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(-ONE,X_r[0],ONE,X_r_ex0,&ierr),ierr);
PHIST_ICHK_IERR(SUBR(mvec_add_mvec)(-ONE,X_i[0],ONE,X_i_ex0,&ierr),ierr);
double nrm_err0_1, nrm_err0_2;
PHIST_ICHK_IERR(SUBR(mvec_dot_mvec)(X_i_ex0,X_i_ex0,&nrm_err0_1,&ierr),ierr);
PHIST_ICHK_IERR(SUBR(mvec_dot_mvec)(X_r_ex0,X_r_ex0,&nrm_err0_2,&ierr),ierr);

PHIST_SOUT(PHIST_VERBOSE,"err in re(x): %e\n",SQRT(nrm_err0_2));
PHIST_SOUT(PHIST_VERBOSE,"       im(x): %e\n",SQRT(nrm_err0_1));

///////////////////////////////////////////////////////////////////
// clean up afterwards                                           //
///////////////////////////////////////////////////////////////////

  PHIST_ICHK_IERR(SUBR(feastCorrectionSolver_delete)(fCorrSolver,&ierr),ierr);

  free(sigma_r);
  free(sigma_i);
#ifdef PHIST_KERNEL_LIB_GHOST
  mat->destroy(mat);
  ghost_context_destroy(ctx);
#else
  PHIST_ICHK_IERR(SUBR(crsMat_delete)(mat,&ierr),ierr);
#endif
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
