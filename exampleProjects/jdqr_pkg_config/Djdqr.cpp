/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"

#include "phist_macros.h"
#include "phist_enums.h"
#include "phist_kernels.h"
#include "phist_operator.h"
#include "phist_subspacejada.h"
#include "phist_jadaOpts.h"
#include "phist_gen_d.h"
#include "phist_driver_utils_decl.h"
#include "phist_ScalarTraits.hpp"
#include "phist_std_typedefs.hpp"

// demo app that creates a test matrix (Quantum harmonic oscillator) and computes some eigenmodes using subspacejada.
int main(int argc, char** argv)
{
  int iflag=0;
  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&iflag),iflag);

  PHIST_MAIN_TASK_BEGIN
  
  phist_comm_ptr comm;
  phist_const_map_ptr map; // map (element distribution) of vectors according to 
  sparseMat_ptr A;
  linearOp_ptr A_op; // this is a wrapper for the CRS matrix which we pass to the actual solver
  linearOp_ptr B_op=NULL; // no mass matrix used in this example
  mvec_ptr Q=NULL;
  sdMat_ptr R=NULL;
  
  phist_jadaOpts opts;

  const char* matrix="SCAMAC-Harmonic";
  
  int num_eigs,nV,num_iters;
  double max_err=0.;
  
  PHIST_ICHK_IERR(phist_comm_create(&comm,&iflag),iflag);

  phist_jadaOpts_setDefaults(&opts);
  
  opts.numEigs=10;
  opts.blockSize=2;
  opts.which=phist_SR;
  opts.convTol=1.0e-8;
  opts.maxIters=250;
  opts.minBas=20;
  opts.maxBas=40;
  opts.innerSolvType=phist_MINRES;

  num_eigs  = opts.numEigs;
  num_iters = 0;

  iflag = PHIST_SPARSEMAT_PERM_GLOBAL;
  PHIST_ICHK_IERR(SUBR(create_matrix)(&A,comm,matrix,&iflag),iflag);
  
  PHIST_ICHK_IERR(SUBR(sparseMat_get_domain_map)(A, &map,&iflag),iflag);

  nV = num_eigs + opts.blockSize-1;
  PHIST_ICHK_IERR(SUBR(mvec_create)(&Q,map,nV,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(mvec_random)(Q,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(sdMat_create)(&R,nV,nV,comm,&iflag),iflag);
  
  // create operator wrapper for computing Y=A*X using a CRS matrix
  A_op = new TYPE(linearOp);
  PHIST_ICHK_IERR(SUBR(linearOp_wrap_sparseMat)(A_op,A,&iflag),iflag);
  
  CT evals[nV];
  MT resid[nV];

  PHIST_SOUT(PHIST_INFO,"Jacobi-Davidson for left-most eigenvalues of a Quantum-mechanical Harmonic Oscillator...\n");
    PHIST_ICHK_NEG_IERR(SUBR(subspacejada)(A_op, B_op, opts,
          Q, R, evals, resid, &num_eigs, &num_iters, &iflag), iflag);

  PHIST_SOUT(PHIST_INFO,"Matrix: '%s'\nComputed %d out of %d desired eigenvalues\nafter %d iterations.\n", matrix, num_eigs, opts.numEigs, num_iters);
  if (num_eigs>0) PHIST_SOUT(PHIST_INFO, "Eigenvalue\tResidual\tError\n");
  for (int i=0; i<num_eigs; i++)
  {
    double err=double(i)-std::real(evals[i]);
    max_err = std::max(max_err, err);
    PHIST_SOUT(PHIST_INFO,"%8.4e\t%8.4e\t%8.4e\n",std::real(evals[i]),resid[i],err);
  }

  if (num_eigs<opts.numEigs || max_err>10.*opts.convTol)
  {
    PHIST_SOUT(PHIST_WARNING,"WARNING: Unexpectedly large error in eigenvalue(s).\n");
    iflag=+1;
  }
  
  PHIST_ICHK_IERR(SUBR(mvec_delete)(Q,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(sdMat_delete)(R,&iflag),iflag);
  PHIST_ICHK_IERR(SUBR(sparseMat_delete)(A,&iflag),iflag);
  delete A_op;

  PHIST_MAIN_TASK_END

  PHIST_ICHK_IERR(phist_kernels_finalize(&iflag),iflag);

  return iflag;
}
