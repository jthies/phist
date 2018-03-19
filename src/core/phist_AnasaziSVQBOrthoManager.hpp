/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

#include "phist_config.h"

#ifdef PHIST_HAVE_ANASAZI

#ifndef _ST_
#error "you must to include a 'phist_gen_X' header before this filee!"
#endif

#include "phist_trilinos_type_config.h"

#ifdef PHIST_TRILINOS_TYPE_AVAIL

// note: we should *not* include any phist headers here because they may override the
// settings from the previously included "phist_gen_X.h" file. So instead we have to 
// include all headers that this file needs in the source files phist_belos*.cpp.

#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_RCP.hpp"

namespace Anasazi {

  template<>
  int SVQBOrthoManager<_ST_, ::phist::BelosMV< _ST_ >, ::phist::types< _ST_ >::linearOp>::findBasis(
                ::phist::BelosMV< _ST_ > &X, Teuchos::RCP<::phist::BelosMV< _ST_ > > MX, 
                Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,_ST_ > > > C, 
                Teuchos::RCP<Teuchos::SerialDenseMatrix<int,_ST_ > > B,
                Teuchos::Array<Teuchos::RCP<const ::phist::BelosMV< _ST_ > > > Q,
                bool normalize_in) const 
  {
    typedef phist::types< _ST_ >   pt;
    typedef phist::kernels< _ST_ > pk;
    typedef phist::core< _ST_ >    pc;
    typedef phist::ScalarTraits< _ST_ > st;
    typedef phist::ScalarTraits< _MT_ > mt;

    // the original Belos::ICGSOrthoManager implementation randomizes the null-space, so we'll do that, too.
    int iflag=0;
    int iflag_in=PHIST_ORTHOG_RANDOMIZE_NULLSPACE;
    
    if (!normalize_in)
    {
      PHIST_SOUT(PHIST_ERROR,"findBasis without normalizing X -> not implemented'\n");
      throw phist::Exception(PHIST_NOT_IMPLEMENTED);
    }

    Teuchos::RCP< phist::BelosMV< _ST_ > > X_or_MX=MX;
    Teuchos::RCP< const phist::BelosMV< _ST_ > > Qi=Teuchos::null;
    pt::const_mvec_ptr vQi=nullptr;
    
    if (MX==Teuchos::null) X_or_MX=Teuchos::rcp(&X,false);
    if (Q.size()>0) {Qi=Q[0]; vQi=Qi->get();}

    phist_const_comm_ptr comm;
    pk::mvec_get_comm(X.get(),&comm,&iflag);
    pt::sdMat_ptr Bphist=NULL, Cphist=NULL, XtMX=NULL;
    int ncolsX = MVT::GetNumberVecs(X);
    int ncolsQ = -1;

    pk::sdMat_create(&Bphist,ncolsX,ncolsX,comm,&iflag);
    pk::sdMat_create(&XtMX,ncolsX,ncolsX,comm,&iflag);
    SdMatOwner< _ST_ > _Bphist(Bphist),_XtMX(XtMX), _Cphist;

    if (Q.size()>1)
    {
      // when orthogonalizing X against Q_i, the previous efforts may be compromized and <Q_j,X>!=0 for j<i may occur,
      // so we would actually have to add another loop here. I don't know exactly which Belos solvers require this 
      // feature, if someone complains we can still add it, of course.
      PHIST_SOUT(PHIST_ERROR,"findBasis with several Q -> not implemented'\n");
      throw phist::Exception(PHIST_NOT_IMPLEMENTED);
    }

    int rank_of_X=ncolsX;
    
    {
      int ncolsQi = (Qi!=Teuchos::null)? MVT::GetNumberVecs(*Qi) : 0;
      if (ncolsQi!=ncolsQ)
      {
        ncolsQ=ncolsQi;
        if (Cphist) pk::sdMat_delete(Cphist,&iflag);
        if (ncolsQ>0) pk::sdMat_create(&Cphist,ncolsQ,ncolsX,comm,&iflag);
        // make sure after the loop the last Cphist is deleted
        _Cphist.set(Cphist);
      }
      // (re-)compute X'*M*X because X is updated for every i
      pk::mvecT_times_mvec(st::one(),X.get(),X_or_MX->get(),st::zero(),XtMX,&iflag);
      int iflag=iflag_in;
      int rankQiX; // rank of [Q[i],X] before orthog (we'll randomize the null-space, so unless iflag=-8 is returned,
                   // [Q[i] X] has full rank afterwards)
      _MT_ rankTol=mt::rankTol();
      _MT_ orthoEps=mt::eps();
      int numSweeps=3;
      try {
      iflag=PHIST_ORTHOG_TRIANGULAR_R1;
      phist::core< _ST_ >::orthog_impl
          (vQi,X.get(),_Op.get(),X_or_MX->get(),
          XtMX,Bphist,Cphist,
          numSweeps,&rankQiX,rankTol,orthoEps,
          &iflag);
      } catch (phist::Exception e)
      {
        // note: this function returns iflag=-8 if the nullspace of X could not be augmented by
        // random vectors. iflag=1 means that [Q X] did not have full rank, but Belos wants us 
        // to return the rank *of X after randomizing the null space*.
        if (iflag==-8) rank_of_X=rankQiX-ncolsQi;
      }

      if (Cphist!=nullptr && C[0]!=Teuchos::null) MVT::CopySdMatToTeuchos(Cphist, *C[0]);
      if (B!=Teuchos::null) MVT::CopySdMatToTeuchos(Bphist, *B);
      if (iflag==1) PHIST_SOUT(PHIST_VERBOSE,"orthog suggests rank([Q,X])=%d on input.\n",rankQiX);
    }
    return rank_of_X;
  }


} // namespace Anasazi

#undef PHIST_TRILINOS_TYPE_AVAIL
#endif // PHIST_TRILINOS_TYPE_AVAIL

#endif // PHIST_HAVE_ANASAZI
