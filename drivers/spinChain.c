#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ghost.h>
#include <ghost_util.h>

// struct for JaDa settings
#include "phist_jadaOpts.h"
// single-vector JDQR solver
#include "phist_jdqr.h"
// for phist iterative solvers we need to provide
// an operator representation of the matrix
#include "phist_operator.h"

#include "matfuncs.h"

GHOST_REGISTER_DT_D(my_datatype)
//   my_datatype    =  ghost_datatype_int
//   my_datatype_t  = cast auf datatype ( float, double, .. )




int main( int argc, char* argv[] )
{



//    if( argc < 4 ) { printf("use: %s <M> <Number_of_Rand_Vecs> <dos_data_file>\n", argv[0]); exit(0); }

//TODO - allow settings to be passed on command line or somehow else
    int numEigs = 5;
    int numIters=300;

    ghost_vtraits_t vtraits = GHOST_VTRAITS_INIT(.flags = GHOST_VEC_LHS|GHOST_VEC_RHS, 
                                                 .datatype = my_datatype, 
                                                 .nvecs=numEigs+1);
    ghost_vtraits_t *vtraits1 = ghost_cloneVtraits(&vtraits);
    vtraits1->nvecs=1;
    
    ghost_mtraits_t mtraits = GHOST_MTRAITS_INIT(.format = GHOST_SPM_FORMAT_CRS, .flags = GHOST_SPM_DEFAULT, .datatype = my_datatype);

    int m,i,ierr;

    ghost_midx_t DIM;
    //TODO - get from command line
    ghost_midx_t conf_spinZ[3] = {20,10,0};
    SpinChainSZ( -2, &DIM, conf_spinZ, NULL);

    matfuncs_info_t info;
    //crsGraphene( -1, NULL, NULL, &info);
    SpinChainSZ( -1, NULL, NULL, &info);

    if ( my_datatype != info.datatype) {
	     printf("error: datatyte does not match\n");  exit(0);
	}

    ghost_mat_t * mat;
    ghost_context_t * ctx;
    ghost_vec_t* eigVecs; // eigenvector approximations
    ghost_vec_t* startVec; // for creating a starting vector for JaDa
    double eigVals[numEigs+1]; // eigenvalue approximations
    double resid[numEigs+1]; // residuals computed by JaDa
    int is_cmplx[numEigs+1]; // JDQR is intended for non-symmetric matrices as well,
                             // this flag array can be ignored here because the
                             // SpinChain matrix is symmetric
    double expRes[numEigs+11]; // explicit residual norms

    ghost_init(argc,argv);

    ctx = ghost_createContext( info.nrows , info.ncols ,GHOST_CONTEXT_DEFAULT,NULL,MPI_COMM_WORLD,1.);

    mat = ghost_createMatrix(ctx,&mtraits,1);

    mat->fromRowFunc( mat, info.row_nnz , 0, &SpinChainSZ, 0);

    ghost_printMatrixInfo(mat);

    // create starting vector (if you don't provide one, a random
    // start vector will be generated, but here we give something
    // so runs are reproducible)
    startVec = ghost_createVector(ctx,vtraits1);
    double one=1.0;
    startVec->fromScalar(startVec,&one);

    // create the eigenvectors
    eigVecs = ghost_createVector(ctx,&vtraits);

    // setup options for JDQR
    phist_jadaOpts_t opts;
    phist_jadaOpts_setDefaults(&opts);
    
    // TODO - expose JaDa parameters to the caller via command-line or input file
    opts.numEigs=numEigs;
    opts.which=SR;
    opts.convTol=1.0e-10;
    opts.maxIters=numIters;
    opts.minBas=10;
    opts.maxBas=30;
    opts.v0=startVec;

    // wrap the matrix into a phist operator
    Dop_t A_op;
    // no mass matrix
    Dop_t* B_op=NULL;
    
    phist_Dop_wrap_crsMat(&A_op,mat,&ierr);

    // run JDQR
    phist_Djdqr(&A_op,B_op,eigVecs,eigVals,resid,is_cmplx,
              opts,&numEigs,&numIters,
              &ierr);

    if (ghost_getRank(MPI_COMM_WORLD)==0)
    {
      fprintf(stdout, "Found %d eigenpair(s) after %d iterations.",numEigs,numIters);
      if (numEigs>0)
      {
        fprintf(stdout,"TODO - expl. residual not yet computed.\n");
        fprintf(stdout,"  Eigenvalue\t\t\t\tRitz Residual\tExpl. Residual\n");
        fprintf(stdout,"=======================================================\n");
        for (i=0;i<numEigs;i++)
        {
          if (is_cmplx[i])
          {
            fprintf(stdout,"WARNING: jdqr reports a complex eigenvalue, something is fishy.\n"
                           "         the eigenvalue will be printed in two consecutive rows\n"
                           "(real/imag part) anyway\n");
          }
          fprintf(stdout,"%24.16e\t%3.1e\t%3.1e\n",eigVals[i],resid[i],expRes[i]);
        }
      }
    }

    mat->destroy(mat);
    eigVecs->destroy(eigVecs);
    ghost_freeContext(ctx);
    free(vtraits1);

    ghost_finish();

    return EXIT_SUCCESS;
}
