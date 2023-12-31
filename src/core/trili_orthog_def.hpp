/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! This file can be used to replace our own orthog routine
//! by one that uses an ortho manager from Belos. This is
//! mostly for testing our own code but can also be used to
//! quickly incorporate new developments in their code.

#ifdef HAVE_TRILINOS_ORTHO_MANAGER
#ifndef DOXYGEN

# include "phist_rcp_helpers.hpp"
# include "phist_operator.h"
# include "phist_BelosOperatorTraits.hpp"

# ifdef PHIST_KERNEL_LIB_TPETRA
#  include "BelosTpetraAdapter.hpp"
#  include "phist_tpetra_typedefs.hpp"
# elif defined(PHIST_KERNEL_LIB_EPETRA)
#  include "Epetra_MultiVector.h"
#  include "BelosEpetraAdapter.hpp"

// we include the cpp file here because epetra_helpers is not
// compiled into a library anywhere??
#  include "epetra_helpers.cpp"

#  include "Epetra_SerialComm.h"
#  include "Epetra_SerialDenseMatrix.h"
#  include "Epetra_CrsMatrix.h"
# endif

# include "BelosOrthoManager.hpp"
// we can try any of the methods implemented in Belos:
#ifdef BELOS_HAVE_TSQR
# include "BelosTsqrOrthoManager.hpp"
#endif
# include "BelosDGKSOrthoManager.hpp"
# include "BelosICGSOrthoManager.hpp"
# include "BelosIMGSOrthoManager.hpp"

#endif //DOXYGEN
#endif /* HAVE_TRILINOS_ORTHO_MANAGER */


void SUBR(trili_orthog)(TYPE(const_mvec_ptr) V,
                     TYPE(mvec_ptr) W,
                     TYPE(const_linearOp_ptr) B,
                     TYPE(sdMat_ptr) R1,
                     TYPE(sdMat_ptr) R2,
                     int numSweeps,
                     int* rankVW,
                     int* iflag)
  {
#ifndef HAVE_TRILINOS_ORTHO_MANAGER
  *iflag = PHIST_NOT_IMPLEMENTED;
#else
  PHIST_ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"


  int m,k;
  int rankW;
  *iflag=0;
  MT rankTol=mt::sqrt(mt::eps());
  
  // all vectors and matrices must be allocated.
  if (V==NULL || W==NULL || R1==NULL || R2==NULL) 
    {
    *iflag=-1;
    return;
    }

  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&m,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(W,&k,iflag),*iflag);
  
  if (k==0) // no vectors to be orthogonalized
    {
    return;
    }

  if (m==0)
    {
    PHIST_CHK_NEG_IERR(SUBR(mvec_QR)(W,R1,iflag),*iflag);
    return;
    }

#ifdef PHIST_TESTING

  // check that all array dimensions are correct
  PHIST_CHK_IERR(*iflag=(int)(V==NULL || W==NULL || R1==NULL || R2==NULL),*iflag);

  phist_lidx n,nrR1,ncR1,nrR2,ncR2,tmp;

  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&n,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_my_length)(W,&tmp,iflag),*iflag);

  PHIST_CHK_IERR(*iflag=((n==tmp)?0:-1),*iflag);

  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(R1,&nrR1,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(R1,&ncR1,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(R2,&nrR2,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(R2,&ncR2,iflag),*iflag);

  // Q (and W) must match R1, R1 must be square
  PHIST_CHK_IERR(*iflag=((k==nrR1)?0:-1),*iflag);
  PHIST_CHK_IERR((*iflag=(k==ncR1)?0:-1),*iflag);
  // V must match R2, and V*R2 must match W
  PHIST_CHK_IERR((*iflag=(m==nrR2)?0:-1),*iflag);
  PHIST_CHK_IERR((*iflag=(k==ncR2)?0:-1),*iflag);

//  PHIST_DEB("orthog: V is %dx%d,  W is %dx%d\n"
//                        "       R1 is %dx%d, R2 is %dx%d\n",n,m,n,k,nrR1,ncR1,nrR2,ncR2);
#endif

#ifdef PHIST_KERNEL_LIB_TPETRA
typedef phist::tpetra::Traits< ST >::mvec_t MV;
typedef phist::tpetra::Traits< ST >::mvec_t BelosMV;
typedef phist::tpetra::Traits< ST >::mvec_t MAT;
typedef phist::tpetra::Traits< ST >::mvec_t BelosMAT;
typedef phist::tpetra::Traits< ST >::Teuchos_sdMat_t TeuchosMAT;
#elif defined(PHIST_KERNEL_LIB_EPETRA)
typedef Epetra_MultiVector MV;
typedef MV MAT;
typedef Epetra_MultiVector BelosMV;
typedef BelosMV BelosMAT;
typedef Teuchos::SerialDenseMatrix<int,ST > TeuchosMAT;
#elif defined(PHIST_KERNEL_LIB_GHOST)
typedef ghost_vec_t MV;
typedef MV MAT;
typedef phist::GhostMV BelosMV;
typedef BelosMV BelosMAT;
typedef Teuchos::SerialDenseMatrix<int,ST > TeuchosMAT;
#else
#error "not implemented"
#endif

Teuchos::RCP<const BelosMV> V_2 = phist::rcp((const MV*)V,false);
Teuchos::RCP<BelosMV> W_2 = phist::rcp((MV*)W,false);
Teuchos::RCP<BelosMAT> R1_1 = phist::rcp((MAT*)R1,false);
Teuchos::RCP<BelosMAT> R2_1 = phist::rcp((MAT*)R2,false);
#ifdef PHIST_KERNEL_LIB_TPETRA
Teuchos::RCP<TeuchosMAT> R1_2 = 
        phist::tpetra::Traits< ST >::CreateTeuchosViewNonConst(R1_1,iflag);
Teuchos::RCP<TeuchosMAT> R2_2 = 
        phist::tpetra::Traits< ST >::CreateTeuchosViewNonConst(R2_1,iflag);
#elif defined(PHIST_KERNEL_LIB_EPETRA)
        Teuchos::RCP<TeuchosMAT> R1_2 = CreateTeuchosViewNonConst(R1_1,iflag);
        Teuchos::RCP<TeuchosMAT> R2_2 = CreateTeuchosViewNonConst(R2_1,iflag);
#else
#error "not implemented"
#endif

  // try the TSQR method
  Teuchos::RCP<Belos::OrthoManager< ST ,MV> > ortho;


  Teuchos::RCP<const Teuchos::ParameterList> default_params, fast_params;

//#define TRY_TSQR
#define TRY_IMGS

  Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp
        (new Teuchos::ParameterList());

#ifdef TRY_TSQR
typedef Belos::TsqrOrthoManager<ST,MV> orthoMan_t;
Teuchos::RCP<orthoMan_t> tsqr = Teuchos::rcp(new orthoMan_t("hist/orthog/tsqr"));
ortho = tsqr;

  // if the matrix is rank-deficient, fill the 'missing' columns with
  // random vectors and return the rank of the original matrix:
  params->set("randomizeNullSpace",true);
  // this is the tolerance for determining rank deficiency. The tolerance
  // is relative to the largest singular value of R in the QR decomp of V,
  // i.e. the 2-norm of R.
//  params->set("relativeRankTolerance",rankTol);

  default_params = tsqr->getValidParameters();
  fast_params = tsqr->getFastParameters();

#elif defined(TRY_IMGS)

typedef Belos::IMGSOrthoManager<ST,MV,TYPE(linearOp)> orthoMan_t;
Teuchos::RCP<orthoMan_t> imgs = Teuchos::rcp(new orthoMan_t("hist/orthog/imgs"));
if (B != NULL) {
  imgs = Teuchos::rcp(new orthoMan_t("hist/orthog/imgs", Teuchos::rcp(B, false)));
}
ortho = imgs;

  default_params = imgs->getValidParameters();
  fast_params = imgs->getFastParameters();

#endif

  // set default parameters but avoid discarding our entries from above:
  params->setParametersNotAlreadySet(*default_params);

  Teuchos::RCP<orthoMan_t> cast = Teuchos::rcp_dynamic_cast<orthoMan_t>(ortho);
  cast->setParameterList(params);

  Teuchos::ArrayView< Teuchos::RCP<const BelosMV > > V_array( &V_2, 1 );
  Teuchos::Array< Teuchos::RCP<TeuchosMAT > > R2_array( 1, R2_2 );
  *rankVW=ortho->projectAndNormalize( *W_2, R2_array, R1_2, V_array );

#ifdef PHIST_KERNEL_LIB_EPETRA
# if PHIST_OUTLEV>=PHIST_DEBUG
  if (B != NULL)
  {
    const Epetra_MultiVector *Vview = (const Epetra_MultiVector *)V;
    Epetra_MultiVector *Wview = (Epetra_MultiVector *)V;

    int m=Vview->NumVectors();
    int k=Wview->NumVectors();

    PHIST_CAST_PTR_FROM_VOID(const Epetra_CrsMatrix,B_ptr,B->A,*iflag);

    Epetra_MultiVector tmp(Wview->Map(), k);

    B_ptr->Multiply(false, *Wview, tmp);
    Epetra_SerialDenseMatrix C(m,k);
    Epetra_SerialComm comm;
    Epetra_LocalMap tinyMap(m,0,comm);
    Epetra_MultiVector Cview(View,tinyMap,C.A(),C.LDA(),k);
    Cview.Multiply('T','N',1.0,*Vview,tmp,0.0);

    PHIST_DEB("V'BV matrix\n");
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < k; j++) {
        PHIST_DEB("%f ", C[j][i]);
      }
      PHIST_DEB("\n");;
    }
  }
# endif
#endif
  *iflag = 0; //ncols-rankW;// return positive number if rank not full.
#endif /* HAVE_TRILINOS_ORTHO_MANAGER */
  }
