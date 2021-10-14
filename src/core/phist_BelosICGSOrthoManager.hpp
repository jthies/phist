/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
//! \file phist_BelosICGSOrthoManager.hpp
//! \brief orthogonalization routine in Belos

#include "phist_config.h"

#ifdef PHIST_HAVE_BELOS

#ifndef _ST_
#error "you must to include a 'phist_gen_X' header before this file!"
#endif

#ifndef DOXYGEN
#include "phist_trilinos_type_config.h"
#endif //DOXYGEN

#ifdef PHIST_TRILINOS_TYPE_AVAIL

// note: we should *not* include any phist headers here because they may override the
// settings from the previously included "phist_gen_X.h" file. So instead we have to 
// include all headers that this file needs in the source files phist_belos*.cpp.

#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_RCP.hpp"

namespace Belos {

  //! overloads the crucial orthogonalization routine in Belos     
  //! for our wrapped mvec objects, replacing it by phist's orthog.
  //! The reason is that Belos' implementation makes heavy use of  
  //! single-column views, which are unnecessary and particularly  
  //! inefficient for kernel libraries using row-major storage.    
  template<>
  int ICGSOrthoManager< _ST_, ::phist::BelosMV< _ST_ >, ::phist::types< _ST_ >::linearOp >::
   projectAndNormalizeWithMxImpl (BelosMV< _ST_ > &X,
                                   ::Teuchos::RCP<::phist::BelosMV< _ST_ > > MX,
                                   ::Teuchos::Array<::Teuchos::RCP<::Teuchos::SerialDenseMatrix<int, _ST_ > > > C,
                                   ::Teuchos::RCP<::Teuchos::SerialDenseMatrix<int, _ST_ > > B,
                                   ::Teuchos::ArrayView<::Teuchos::RCP<const ::phist::BelosMV< _ST_ > > > Q) const
  {
    typedef phist::types< _ST_ >   pt;
    typedef phist::kernels< _ST_ > pk;
    typedef phist::core< _ST_ >    pc;
    typedef phist::ScalarTraits< _ST_ > st;

    // the original Belos::ICGSOrthoManager implementation randomizes the null-space, so we'll do that, too.
    int iflag=0;
    int iflag_in=PHIST_ORTHOG_RANDOMIZE_NULLSPACE;
    
    Teuchos::RCP< phist::BelosMV< _ST_ > > X_or_MX=MX;
    if (MX==Teuchos::null) X_or_MX=Teuchos::rcp(&X,false);
    
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
      PHIST_SOUT(PHIST_ERROR,"projectAndNormalize with several Q -> not implemented'\n");
      throw phist::Exception(PHIST_NOT_IMPLEMENTED);
    }

    if (Q.size()==0)
    {
      PHIST_SOUT(PHIST_ERROR,"projectAndNormalize without Q -> only 'normalize (not implemented)'\n");
      throw phist::Exception(PHIST_NOT_IMPLEMENTED);
    }

    int rank_of_X=ncolsX;
    
    // in contrast to orthog, this routine allows orthogonalizing against
    // several blocks. Hence the outer loop
    for (int i=0; i<Q.size(); i++)
    {
      
      int ncolsQi = MVT::GetNumberVecs(*Q[i]);
      if (ncolsQi!=ncolsQ)
      {
        ncolsQ=ncolsQi;
        if (Cphist) pk::sdMat_delete(Cphist,&iflag);
        pk::sdMat_create(&Cphist,ncolsQ,ncolsX,comm,&iflag);
        // make sure after the loop the last Cphist is deleted
        _Cphist.set(Cphist);
      }
      // (re-)compute X'*M*X because X is updated for every i
      pk::mvecT_times_mvec(st::one(),X.get(),X_or_MX->get(),st::zero(),XtMX,&iflag);
      int iflag=iflag_in;
      int rankQiX; // rank of [Q[i],X] before orthog (we'll randomize the null-space, so unless iflag=-8 is returned,
                   // [Q[i] X] has full rank afterwards)
      _MT_ rankTol=_MT_(sing_tol_);
      _MT_ orthoEps=_MT_(blk_tol_);
      int numSweeps=max_ortho_steps_;
      try {
      iflag|=PHIST_ORTHOG_TRIANGULAR_R1;
      phist::core< _ST_ >::orthog_impl
          (Q[i]->get(),X.get(),_Op.get(),X_or_MX->get(),
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

      MVT::CopySdMatToTeuchos(Cphist, *C[i]);
      MVT::CopySdMatToTeuchos(Bphist, *B);
      if (iflag==1) PHIST_SOUT(PHIST_VERBOSE,"orthog suggests rank([Q,X])=%d on input.\n",rankQiX);
    }
    return rank_of_X;
  }
} /* namespace Belos */
#undef PHIST_TRILINOS_TYPE_AVAIL
#endif /* PHIST_TRILINOS_TYPE_AVAIL */
#endif /* PHIST_HAVE_BELOS */

