#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ghost.h>
#include <ghost_util.h>

// struct for JaDa settings
#include "phist_jadaOpts.h"
// single-vector JDQR solver
#include "phist_jdqr.h"
//TROET
#include "phist_bgmres.h"
// for phist iterative solvers we need to provide
// an operator representation of the matrix
#include "phist_operator.h"

#include "matfuncs.h"

GHOST_REGISTER_DT_D(my_datatype)
//   my_datatype    =  ghost_datatype_int
//   my_datatype_t  = cast auf datatype ( float, double, .. )




int main( int argc, char* argv[] )
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

    int L=atoi(argv[1]);
    if ((double)(L/2.0)*2!=(double)L)
    {
        printf("parameter L must be even for this driver.\n");
        printf("Type ./spinChain (without args) to get a usage message\n");
        exit(-1);
    }
    int numEigs=5;
    int numIters;
    if (argc>2)
    {
        numEigs=atoi(argv[2]);
    }
    if (numEigs<=0)
    {
        printf("parameter numEigs must be >0 for this driver.\n");
        printf("Type ./spinChain (without args) to get a usage message\n");
        exit(-1);
    }
    ghost_vtraits_t vtraits = GHOST_VTRAITS_INIT(.flags = GHOST_VEC_LHS|GHOST_VEC_RHS, 
                                                 .datatype = my_datatype, 
                                                 .nvecs=numEigs+1);
    ghost_vtraits_t *vtraits1 = ghost_cloneVtraits(&vtraits);
    vtraits1->nvecs=1;
    
    ghost_mtraits_t mtraits = GHOST_MTRAITS_INIT(.format = GHOST_SPM_FORMAT_CRS, .flags = GHOST_SPM_DEFAULT, .datatype = my_datatype);

    int m,i,ierr;

    ghost_midx_t DIM;
    //TODO - get from command line
    ghost_midx_t conf_spinZ[3] = {L,L/2,0};
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
    ghost_vec_t* resVecs; // explicitly compute residuals after JaDa run
    ghost_vec_t* startVec; // for creating a starting vector for JaDa
    double eigVals[numEigs+1]; // eigenvalue approximations
    double resid[numEigs+1]; // residuals computed by JaDa
    int is_cmplx[numEigs+1]; // JDQR is intended for non-symmetric matrices as well,
                             // this flag array can be ignored here because the
                             // SpinChain matrix is symmetric
    double expRes[numEigs+1]; // explicit residual norms

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
    // allocate memory and set starting vector to 1
    startVec->fromScalar(startVec,&one);

    // create the eigenvectors
    eigVecs = ghost_createVector(ctx,&vtraits);
    // allocate memory for eigenvectors
    double zero=0.0;
    eigVecs->fromScalar(eigVecs,&zero);

    // setup options for JDQR
    phist_jadaOpts_t opts;
    phist_jadaOpts_setDefaults(&opts);
    
    // TODO - expose JaDa parameters to the caller via command-line or input file
    opts.numEigs=numEigs;
    opts.which=SR; // we assume that we always look for smallest real part eigs in this application
    if (argc>=4)
    {
      opts.convTol=atof(argv[3]);
    }
    if (argc>=5)
    {
        opts.maxIters=atoi(argv[4]);
    }
    if (argc>=6)
    {
        opts.minBas=atoi(argv[5]);
        if (argc>=7)
            {
                opts.maxBas=atoi(argv[6]);
            }
        else
            {
                opts.maxBas=opts.minBas+10;
            }
    }
    if (argc>=8)
    {
        opts.arno=atoi(argv[7]);
    }
    if (argc>=9)
    {
        opts.initialShift=atof(argv[8]);
    }
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
    // compute R=A*X-X*Lambda

    resVecs = ghost_createVector(ctx,&vtraits);
    resVecs->fromVec(resVecs,eigVecs,0);
    resVecs->vscale(resVecs,eigVals);
    double minus_one=-1.0;
    resVecs->scale(resVecs,&minus_one);
    int spmvmOpts = GHOST_SPMVM_AXPY;
    ghost_spmvm(ctx,resVecs,mat,eigVecs,&spmvmOpts);
    resVecs->dotProduct(resVecs,resVecs,expRes);
    for (i=0;i<numEigs;i++)
    {
      expRes[i]=sqrt(expRes[i]);
    }


    if (ghost_getRank(MPI_COMM_WORLD)==0)
    {
      fprintf(stdout, "Found %d eigenpair(s) after %d iterations.\n",numEigs,numIters);
      if (numEigs>0)
      {
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
          fprintf(stdout,"%24.16e\t\t%3.1e\t\t%3.1e\n",eigVals[i],resid[i],expRes[i]);
        }
      }
    }

    mat->destroy(mat);
    eigVecs->destroy(eigVecs);
    eigVecs->destroy(resVecs);
    ghost_freeContext(ctx);
    free(vtraits1);

    ghost_finish();

    return EXIT_SUCCESS;
}
