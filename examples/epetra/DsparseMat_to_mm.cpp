/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"
#include "phist_tools.h"
#include "phist_kernels.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef PHIST_KERNEL_LIB_EPETRA
#include "Epetra_CrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#elif defined(PHIST_KERNEL_LIB_GHOST)
#include "ghost/sparsemat.h"
#endif

// in order to use the driver utils header, we first have to
// specify which data type we will use (real double (D) here)
#include "phist_gen_d.h"
#include "phist_driver_utils.h"

/*! this is a simple program that can be used to dump one of our
internal test cases as a MatrixMarket file. It is only intended
for the Epetra kernel lib, which as a simple mechanism for doing that
(it could be implemented in Tpetra as easily)
*/
int main(int argc, char** argv)
{
  int iflag;
  phist_comm_ptr comm_world = NULL;
  phist_DsparseMat_ptr A;
  phist_Dmvec_ptr X=NULL, B=NULL;
  phist_const_map_ptr map=NULL;

  PHIST_ICHK_IERR(phist_kernels_init(&argc,&argv,&iflag),iflag);

  // this generates the default comm object, typically MPI_COMM_WORLD for MPI-based kernel libs
  PHIST_ICHK_IERR(phist_comm_create(&comm_world,&iflag),iflag);

  // print some usage info if the program is called without arguments
  if (argc<3)
  {
    PHIST_SOUT(PHIST_ERROR,"Usage: %s <matrix> <output_file> [<num vecs>]\n"
                           "       \"matrix\" is the matrix to be created,\n\n"
                           "   and \"output_file\" is the name of the output ASCII MatrixMarket file to be written\n"
                           "       if \"num vecs\" (an integer) is provided, block vectors x and b are generated with x the solution of Ax=b. The vectors are stored in sol.mm and rhs.mm, respectively.\n",
                                argv[0]);
        // we use this function below to read/generate the matrix, and calling it
        // with the "usage" string makes it print a useful usage message.
        phist_Dcreate_matrix(NULL,NULL,"usage",&iflag);
    return 1;
  }

  const char* matname =argv[1];
  const char* filename=argv[2];
  int nvecs=0;
  if (argc>3) nvecs=atoi(argv[3]);

  // generate test case (or read file)
  PHIST_ICHK_IERR(phist_Dcreate_matrix(&A,comm_world,matname,&iflag),iflag);
  
  if (nvecs>0)
  {
    PHIST_ICHK_IERR(phist_DsparseMat_get_range_map(A,&map,&iflag),iflag);
    PHIST_ICHK_IERR(phist_Dmvec_create(&X,map,nvecs,&iflag),iflag);
    PHIST_ICHK_IERR(phist_Dmvec_create(&B,map,nvecs,&iflag),iflag);
    PHIST_ICHK_IERR(phist_Dcreate_sol_and_rhs(matname,A,X,B,&iflag),iflag);
  }

#ifdef PHIST_KERNEL_LIB_EPETRA
  Epetra_CrsMatrix* A_epetra = (Epetra_CrsMatrix*)A;
  PHIST_ICHK_IERR(iflag=EpetraExt::RowMatrixToMatrixMarketFile( filename, 
        *A_epetra, matname,"generated by PHIST, https://bitbucket.org/essex/phist"),iflag);

  if (nvecs>0)
  {
    Epetra_MultiVector* X_epetra = (Epetra_MultiVector*)X;
    Epetra_MultiVector* B_epetra = (Epetra_MultiVector*)B;
  PHIST_ICHK_IERR(iflag=EpetraExt::MultiVectorToMatrixMarketFile( "sol.mm", 
        *X_epetra, matname,"solution vectors, generated by PHIST, https://bitbucket.org/essex/phist"),iflag);
  PHIST_ICHK_IERR(iflag=EpetraExt::MultiVectorToMatrixMarketFile( "rhs.mm", 
        *B_epetra, matname,"rhs vectors, generated by PHIST, https://bitbucket.org/essex/phist"),iflag);
  }
#elif defined(PHIST_KERNEL_LIB_GHOST)
  ghost_sparsemat* ghost_A = (ghost_sparsemat*)A;
  ghost_sparsemat_to_mm((char*)filename,ghost_A);
#else
  PHIST_SOUT(PHIST_ERROR,"This driver is only meant for the ghost and epetra kernel libs");
  return -99;
#endif
    
  PHIST_ICHK_IERR(phist_DsparseMat_delete(A,&iflag),iflag);
  if (X!=NULL)
  {
    PHIST_ICHK_IERR(phist_Dmvec_delete(X,&iflag),iflag);
  }
  if (B!=NULL)
  {
    PHIST_ICHK_IERR(phist_Dmvec_delete(B,&iflag),iflag);
  }
  PHIST_ICHK_IERR(phist_kernels_finalize(&iflag),iflag);
  }
