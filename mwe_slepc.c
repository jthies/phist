/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  minimal working example with SLEPc
 *  \details Most of the code is borrowed from the SLEPc examples.
 *           The matrix is generated row-wise with ScaMaC routines,
 *           and assembled with the PETSc mechanism. Finally, a SLEPc eigensolver is called.
 *  \ingroup mwe
 */


static char help[] = "Standard eigenproblem based on the ScaMaC matrices.\n\n"
                     "The only command line argument should be:\n"
                     "  matrix name[,parameters], e.g., Hubbard,n_sites=10,n_fermions=5\n\n";

#include <stdbool.h>
#include <slepceps.h>
#include <scamac.h>

int main(int argc,char **argv) {
  Mat            A;           /* problem matrix */
  EPS            eps;         /* eigenproblem solver context */
  EPSType        type;
  PetscReal      error,tol,re,im;
  PetscScalar    kr,ki;
  Vec            xr,xi;
  PetscInt       n=10,i,II,Istart,Iend,nev,maxit,its,nconv;
  PetscErrorCode ierr;


  SlepcInitialize(&argc,&argv,(char*)0,help);

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n Symmetric eigenproblem\n\n"); CHKERRQ(ierr);

  char *examplestr = argv[1];

  ScamacErrorCode err;
  ScamacGenerator * my_gen = NULL;
  char *scamac_str = NULL;

  /* this is the first call to ScaMaC in this program */

  err = scamac_parse_argstr(examplestr, &my_gen, &scamac_str);
  if (err) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Problem with example string:\n%s\n",scamac_str); CHKERRQ(ierr);
    SETERRQ(PETSC_COMM_WORLD,1,"Abort due to error reported by ScaMaC");
  }
  free(scamac_str);
  scamac_str=NULL;

#if !defined(PETSC_USE_COMPLEX)
  if (scamac_generator_query_valtype(my_gen) == SCAMAC_VAL_COMPLEX) {
    SETERRQ(PETSC_COMM_WORLD,1,"ScaMaC example requires complex numbers, but your PETSc version does not support them.");
  }
#endif

  // remember for SLEPC solver
  bool is_symmetric;
  if ((scamac_generator_query_symmetry(my_gen) == SCAMAC_SYMMETRIC) || (scamac_generator_query_symmetry(my_gen) == SCAMAC_HERMITIAN)) {
    is_symmetric = true;
  } else {
    is_symmetric = false;
  }

  err = scamac_generator_check(my_gen, &scamac_str);
  if (err) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Problem with example parameters:\n%s\n",scamac_str); CHKERRQ(ierr);
    SETERRQ(PETSC_COMM_WORLD,1,"Abort due to error reported by ScaMaC");
  }
  free(scamac_str);
  scamac_str=NULL;

  err = scamac_generator_finalize(my_gen);
  if (err) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Problem: scamac_generator_finalize failed [%s]\n",scamac_error_desc(err)); CHKERRQ(ierr);
    SETERRQ(PETSC_COMM_WORLD,1,"Abort due to error reported by ScaMaC");
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Example\n-------\n%s\n",scamac_generator_query_name(my_gen)); CHKERRQ(ierr);

  err=scamac_generator_parameter_desc(my_gen, "desc", &scamac_str);
  if (err) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Problem: scamac_generator_parameter_desc failed [%s]\n",scamac_error_desc(err)); CHKERRQ(ierr);
    SETERRQ(PETSC_COMM_WORLD,1,"Abort due to error reported by ScaMaC");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nParameters\n----------\n%s\n\n",scamac_str); CHKERRQ(ierr);
  free(scamac_str);
  scamac_str=NULL;



  ScamacIdx nrow = scamac_generator_query_nrow(my_gen);
  ScamacIdx maxnzr = scamac_generator_query_maxnzrow(my_gen);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"matrix dimension: %D\n",nrow); CHKERRQ(ierr);

  /* Assemble the matrix */

  ScamacWorkspace * my_ws = NULL;
  err=scamac_workspace_alloc(my_gen, &my_ws);
  if (err) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Problem: scamac_workspace_alloc failed [%s]\n",scamac_error_desc(err)); CHKERRQ(ierr);
    SETERRQ(PETSC_COMM_WORLD,1,"Abort due to error reported by ScaMaC");
  }

  ScamacIdx *cind = NULL;
  double * val = NULL;

  err=scamac_alloc_cind_val(my_gen,SCAMAC_DEFAULT,&cind,&val);
  if (err) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Problem: scamac_alloc_cind_val failed [%s]\n",scamac_error_desc(err)); CHKERRQ(ierr);
    SETERRQ(PETSC_COMM_WORLD,1,"Abort due to error reported by ScaMaC");
  }


  ierr = MatCreate(PETSC_COMM_WORLD,&A); CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,nrow,nrow); CHKERRQ(ierr);
  ierr = MatSetFromOptions(A); CHKERRQ(ierr);

  MatMPIAIJSetPreallocation(A,maxnzr,NULL,maxnzr,NULL);
  MatSeqAIJSetPreallocation(A,maxnzr,NULL);
  MatSeqSBAIJSetPreallocation(A,1,maxnzr,NULL);

  ierr = MatSetUp(A); CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend); CHKERRQ(ierr);

  for (II=Istart; II<Iend; II++) {
    ScamacIdx nzr;
    err = scamac_generate_row(my_gen,my_ws,II,SCAMAC_DEFAULT,&nzr,cind,val);
    if (err) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Problem: scamac_generate_row failed [%s]\n",scamac_error_desc(err)); CHKERRQ(ierr);
      SETERRQ(PETSC_COMM_WORLD,1,"Abort due to error reported by ScaMaC");
    }

    ScamacIdx i;
    if (scamac_generator_query_valtype(my_gen) == SCAMAC_VAL_REAL) {
      for (i=0; i<nzr; i++) {
        // or use MatSetValues to add entries in one go
        err = MatSetValue(A,II,cind[i],val[i],INSERT_VALUES); CHKERRQ(ierr);
      }
    } else {// SCAMAC_VAL_COMPLEX
      for (i=0; i<nzr; i++) {
        // or use MatSetValues to add entries in one go
        err = MatSetValue(A,II,cind[i],val[2*i]+PETSC_i*val[2*i+1],INSERT_VALUES); CHKERRQ(ierr);
      }
    }
  }

  /* free ScaMaC workspace */
  free(cind);
  free(val);

  err = scamac_workspace_free(my_ws);
  if (err) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Problem: scamac_workspace_free failed [%s]\n",scamac_error_desc(err)); CHKERRQ(ierr);
    SETERRQ(PETSC_COMM_WORLD,1,"Abort due to error reported by ScaMaC");
  }

  err = scamac_generator_destroy(my_gen);
  if (err) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Problem: scamac_generator_destroy failed [%s]\n",scamac_error_desc(err)); CHKERRQ(ierr);
    SETERRQ(PETSC_COMM_WORLD,1,"Abort due to error reported by ScaMaC");
  }

  /* this was the last call to ScaMaC in this program */


  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  ierr = MatCreateVecs(A,NULL,&xr); CHKERRQ(ierr);
  ierr = MatCreateVecs(A,NULL,&xi); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nBegin iterations\n\n",n); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps); CHKERRQ(ierr);

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  ierr = EPSSetOperators(eps,A,NULL); CHKERRQ(ierr);
  if (is_symmetric) {
    ierr = EPSSetProblemType(eps,EPS_HEP); CHKERRQ(ierr);
  } else {
    ierr = EPSSetProblemType(eps,EPS_NHEP); CHKERRQ(ierr);
  }
  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps); CHKERRQ(ierr);
  /*
     Optional: Get some information from the solver and display it
  */
  ierr = EPSGetIterationNumber(eps,&its); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its); CHKERRQ(ierr);
  ierr = EPSGetType(eps,&type); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type); CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev); CHKERRQ(ierr);
  ierr = EPSGetTolerances(eps,&tol,&maxit); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit); CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Get number of converged approximate eigenpairs
  */
  ierr = EPSGetConverged(eps,&nconv); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv); CHKERRQ(ierr);

  if (nconv>0) {
    /*
       Display eigenvalues and relative errors
    */
    ierr = PetscPrintf(PETSC_COMM_WORLD,
                       "           k          ||Ax-kx||/||kx||\n"
                       "   ----------------- ------------------\n"); CHKERRQ(ierr);

    for (i=0; i<nconv; i++) {
      /*
        Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
        ki (imaginary part)
      */
      ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi); CHKERRQ(ierr);
      /*
         Compute the relative error associated to each eigenpair
      */
      ierr = EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error); CHKERRQ(ierr);

#if defined(PETSC_USE_COMPLEX)
      re = PetscRealPart(kr);
      im = PetscImaginaryPart(kr);
#else
      re = kr;
      im = ki;
#endif
      if (im!=0.0) {
        ierr = PetscPrintf(PETSC_COMM_WORLD," %9f%+9fi %12g\n",(double)re,(double)im,(double)error); CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g\n",(double)re,(double)error); CHKERRQ(ierr);
      }
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n"); CHKERRQ(ierr);
  }

  /*
     Free work space
  */
  ierr = EPSDestroy(&eps); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&xr); CHKERRQ(ierr);
  ierr = VecDestroy(&xi); CHKERRQ(ierr);
  ierr = SlepcFinalize();

  return ierr;

}
