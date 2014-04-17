#include "phist_config.h"

#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <string>
#include <sstream>
#include <iostream>

#include "phist_macros.h"
#include "phist_enums.h"
#include "phist_kernels.h"
#include "phist_operator.h"
#include "phist_subspacejada.h"

// ghost/spinChain stuff
#include "ghost.h"
#include "ghost/util.h"
extern "C" {
#include "matfuncs.h"
}
#ifdef PHIST_KERNEL_LIB_GHOST
extern "C" {
void init_mtraits(ghost_sparsemat_traits_t* mtraits);
}
#endif

GHOST_REGISTER_DT_D(my_datatype)


// double precision type
#include "phist_gen_d.h"
#include "phist_ScalarTraits.hpp"
#include "phist_std_typedefs.hpp"



// small c++ string helper function
bool endsWith(const std::string& s, const std::string& suffix)
{
  return s.rfind(suffix) == (s.size()-suffix.size());
}


int main(int argc, char** argv)
{
  int ierr;
  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&ierr),ierr);



  //------------------------------- parse input parameters ------------------------- 

  // check arguments
  if( argc < 2 )
  {
    PHIST_SOUT(PHIST_ERROR,"Usage: %s <number of spins> [<invariant subspace dimension>] [<which>] [<tol>] [<max iters> <JD block size> <min basis> <max basis> <GMRES block size> <GMRES subspace dimension> <initial shift> <initial shift iterations> <GMRES use IMGS> <GMRES abort when first converged in block>]\n", argv[0]);
    return 1;
  }

  // number of spins to generate the matrix
  int nSpins = 22;
  if( argc > 1 )
  {
    std::istringstream iss(argv[1]);
    iss >> nSpins;
    if( !(nSpins > 0 && nSpins % 2 == 0) )
    {
      PHIST_SOUT(PHIST_ERROR,"Error: number of spins must be positive even number!\n");
      return 1;
    }
  }

  const std::string filename_A(argv[1]);

  // number of eigenvalues to compute
  int nEig = 10;
  if( argc > 2 )
  {
    std::istringstream iss(argv[2]);
    iss >> nEig;
  }

  // which eigenvalues to compute
  eigSort_t which = SR;
  if( argc > 3 )
  {
    if( argv[3] == std::string("LM") )
      which = LM;
    else if( argv[3] == std::string("SM") )
      which = SM;
    else if( argv[3] == std::string("LR") )
      which = LR;
    else if( argv[3] == std::string("SR") )
      which = SR;
    else
    {
      PHIST_SOUT(PHIST_ERROR,"error parsing argument 3 <which> ('%s'), it should be one of LM, SM, LR or SR\n", argv[4]);
      return 1;
    }
  }

  // desired accuracy: residuum tolerance
  _MT_ tol = mt::eps()*10000;
  if( argc > 4 )
  {
    std::istringstream iss(argv[4]);
    iss >> tol;
    if( tol <= mt::eps() )
    {
      PHIST_SOUT(PHIST_WARNING,"specified tolerance %e is too small! (eps: %e)\n", tol, mt::eps());
    }
  }

  // maximum number of iterations
  int nIter = 250;
  if( argc > 5 )
  {
    std::istringstream iss(argv[5]);
    iss >> nIter;
  }

  // block size
  int blockDim = 4;
  if( argc > 6 )
  {
    std::istringstream iss(argv[6]);
    iss >> blockDim;
  }

  // min basis size
  int minBase = std::min(20,nEig+blockDim-1);
  if( argc > 7 )
  {
    std::istringstream iss(argv[7]);
    iss >> minBase;
  }

  // max basis size
  int maxBase = std::max(80, minBase+blockDim);
  if( argc > 8 )
  {
    std::istringstream iss(argv[8]);
    iss >> maxBase;
  }

  // inner GMRES block size
  int innerBlockDim = std::min(2, blockDim);
  if( argc > 9 )
  {
    std::istringstream iss(argv[9]);
    iss >> innerBlockDim;
  }

  // inner GMRES subspace dimension
  int innerMaxBase = 25;
  if( argc > 10 )
  {
    std::istringstream iss(argv[10]);
    iss >> innerMaxBase;
  }

  // initial shift
  _ST_ initialShift = st::zero();
  if( argc > 11 )
  {
    std::istringstream iss(argv[11]);
    iss >> initialShift;
  }

  // number of initial iterations with fixed initial shift
  int initialShiftIters = 0;
  if( argc > 12 )
  {
    std::istringstream iss(argv[12]);
    iss >> initialShiftIters;
  }

  // in the inner GMRES: use a iterated modified gram schmidt (more accurate, but quite costly!)
  bool innerIMGS = true;
  if( argc > 13 )
  {
    std::istringstream iss(argv[13]);
    iss >> innerIMGS;
  }

  // in the inner GMRES: abort solving a block when first system converged (effectively disables pipelining, but only "blocking" may be faster for now)
  bool innerGMRESabortAfterFirstConverged = false;
  if( argc > 14 )
  {
    std::istringstream iss(argv[14]);
    iss >> innerGMRESabortAfterFirstConverged;
  }



  //------------------------------- setup matrices and vectors --------------------- 
#ifdef PHIST_KERNEL_LIB_GHOST
  ghost_context_t * ctx;
  ghost_sparsemat_t *mat = NULL;
#else
  TYPE(crsMat_ptr) mat = NULL;
#endif
  ghost_idx_t DIM;
  //TODO - get from command line
  ghost_idx_t conf_spinZ[3] = {nSpins,nSpins/2,0};
  SpinChainSZ( -2, &DIM, conf_spinZ, NULL);

  matfuncs_info_t info;
  //crsGraphene( -1, NULL, NULL, &info);
  SpinChainSZ( -1, NULL, NULL, &info);

  if ( my_datatype != info.datatype)
  {
     printf("error: datatyte does not match\n");
     exit(0);
  }

#ifdef PHIST_KERNEL_LIB_FORTRAN
  PHIST_ICHK_IERR(SUBR(crsMat_create_fromRowFunc)(&mat,
        info.nrows, info.ncols, info.row_nnz,
        (void(*)(ghost_idx_t,ghost_idx_t*,ghost_idx_t*,void*))&SpinChainSZ, &ierr), ierr);
#endif

#ifdef PHIST_KERNEL_LIB_GHOST
  ghost_error_t err=ghost_context_create(&ctx, info.nrows , 
      info.ncols,GHOST_CONTEXT_DEFAULT,NULL,GHOST_SPARSEMAT_SRC_NONE,MPI_COMM_WORLD,1.);
  if (err!=GHOST_SUCCESS)
  {
    PHIST_OUT(PHIST_ERROR,"error returned from createContext (file %s, line %d)",__FILE__,__LINE__);
  }        

  ghost_sparsemat_traits_t mtraits;
  init_mtraits(&mtraits);
  ghost_sparsemat_create(&mat, ctx,&mtraits,1);
  ghost_sparsemat_src_rowfunc_t src = GHOST_SPARSEMAT_SRC_ROWFUNC_INITIALIZER;
  src.func = &SpinChainSZ;
  src.maxrowlen = info.row_nnz;
  mat->fromRowFunc(mat, &src);
  char *str;
  ghost_sparsemat_string(&str,mat);
  printf("%s\n",str);
  free(str); str = NULL;
#endif

  // create an operator from A
  op_ptr_t opA = new TYPE(op);
  PHIST_ICHK_IERR(SUBR(op_wrap_crsMat)(opA,mat,&ierr),ierr);

  // we need the domain map of the matrix
  const_comm_ptr_t comm = NULL;
  PHIST_ICHK_IERR(phist_map_get_comm(opA->domain_map,&comm,&ierr),ierr);

  // setup necessary vectors and matrices for the schur form
  mvec_ptr_t Q = NULL;
  PHIST_ICHK_IERR(SUBR(mvec_create)(&Q,opA->domain_map,nEig+blockDim-1,&ierr),ierr);
  sdMat_ptr_t R = NULL;
  PHIST_ICHK_IERR(SUBR(sdMat_create)(&R,nEig+blockDim-1,nEig+blockDim-1,comm,&ierr),ierr);
  _MT_ *resNorm = new _MT_[nEig+blockDim-1];

  // setup start vector (currently to (1 0 1 0 .. ) )
  mvec_ptr_t v0 = NULL;
  PHIST_ICHK_IERR(SUBR(mvec_create)(&v0,opA->domain_map,1,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_put_value)(v0,st::one(),&ierr),ierr);

  // used to calculate explicit residuals
  mvec_ptr_t res;
  PHIST_ICHK_IERR(SUBR(mvec_create)(&res,opA->domain_map,nEig+blockDim-1,&ierr),ierr);


  //------------------------------- run block JaDa algorithm ----------------------- 
  PHIST_ICHK_IERR(SUBR(subspacejada)(opA, NULL, v0, which, tol, nEig, &nIter, blockDim, minBase, maxBase, innerBlockDim, innerMaxBase, initialShiftIters, initialShift, innerIMGS, innerGMRESabortAfterFirstConverged, Q, R, resNorm, &ierr), ierr);
  int nConvergedEig = 0;
  for(int i = 0; i < nEig; i++)
    if(resNorm[i] <= tol)
      nConvergedEig++;
  PHIST_SOUT(PHIST_INFO, "subspacejada terminated after %d iterations and calculated %d eigenvalues;\n printing R:\n", nIter, nConvergedEig);
#ifdef PHIST_HAVE_MPI
  int me;
  PHIST_ICHK_IERR(ierr = MPI_Comm_rank(MPI_COMM_WORLD,&me), ierr);
  if( me == 0 )
  {
    PHIST_ICHK_IERR(SUBR( sdMat_print ) (R, &ierr), ierr);
  }
#else
  PHIST_ICHK_IERR(SUBR( sdMat_print ) (R, &ierr), ierr);
#endif
  PHIST_SOUT(PHIST_INFO, " residuum norm:");
  for(int i = 0; i < nEig; i++)
    PHIST_SOUT(PHIST_INFO, "\t%8.4e", resNorm[i]);

  // calculate real residual
  PHIST_ICHK_IERR(SUBR(crsMat_times_mvec)(st::one(),mat,Q,st::zero(),res,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_times_sdMat)(-st::one(),Q,R,st::one(),res,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_norm2)(res,resNorm,&ierr),ierr);
  PHIST_SOUT(PHIST_INFO, "\n explicit residuum norm:");
  for(int i = 0; i < nEig; i++)
    PHIST_SOUT(PHIST_INFO, "\t%8.4e", resNorm[i]);
  PHIST_SOUT(PHIST_INFO, "\n");




  //------------------------------- clear matrices and vectors --------------------- 
  // delete vectors and sdMats
  PHIST_ICHK_IERR(SUBR(mvec_delete)(res,&ierr),ierr);
  delete[] resNorm;
  PHIST_ICHK_IERR(SUBR(mvec_delete)(v0,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(sdMat_delete)(R,&ierr),ierr);
  PHIST_ICHK_IERR(SUBR(mvec_delete)(Q,&ierr),ierr);
  // clean up operator
  delete opA;
  // delete matrix
#ifdef PHIST_KERNEL_LIB_GHOST
  mat->destroy(mat);
  ghost_context_destroy(ctx);
#else
  PHIST_ICHK_IERR(SUBR(crsMat_delete)(mat,&ierr),ierr);
#endif

  PHIST_ICHK_IERR(phist_kernels_finalize(&ierr),ierr);
  return ierr;
}
