/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef BELOS_PHIST_ADAPTER_HPP
#define BELOS_PHIST_ADAPTER_HPP

#include "phist_config.h"

#ifdef PHIST_HAVE_BELOS

#include "phist_types.hpp"
#include "phist_kernels.hpp"
#include "phist_ScalarTraits.hpp"

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Array.hpp>

#include <BelosConfigDefs.hpp>
#include <BelosTypes.hpp>
#include <BelosMultiVecTraits.hpp>

#include "phist_BelosMV.hpp"

namespace Belos {

using ::phist::BelosMV;


  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::MultiVecTraits for phist mvecs.
  //
  ////////////////////////////////////////////////////////////////////

  /*!  \brief Template specialization of Belos::MultiVecTraits class using the ghost_densemat class.

    This interface will ensure that any phist mvec will be accepted by the Belos
    templated solvers, regardless of what kernel library is used.
    
    Note: Belos requires 'scattered' views, i.e. views of a certain (possibly discontinuous and permuted) subset of 
    columns of an mvec, we throw an exception in this interface if such a request is made.
    */
  template<class ST>
  class MultiVecTraits<ST, BelosMV<ST>  >
  {
  public:

    typedef ::phist::ScalarTraits<ST> st;
    typedef ::phist::types<ST> tt;
    typedef ::phist::kernels<ST> kt;
    
    typedef typename tt::mvec_ptr mvec_ptr;
    typedef typename tt::const_mvec_ptr const_mvec_ptr;
    typedef typename tt::sdMat_ptr sdMat_ptr;
    typedef typename tt::const_sdMat_ptr const_sdMat_ptr;
    
    typedef typename st::magn_t magn_t;
    typedef typename Teuchos::SerialDenseMatrix<int,ST> Teuchos_sdMat;
    
    static inline Teuchos::RCP<BelosMV<ST> > mvec_rcp(mvec_ptr V,bool own_mem) {return ::phist::mvec_rcp<ST>(V,own_mem);}
    static inline Teuchos::RCP<const BelosMV<ST> > mvec_rcp(const_mvec_ptr V,bool own_mem) {return ::phist::mvec_rcp<ST>(V,own_mem);}

  static int GetNumberVecs( const BelosMV<ST>& mv )
  {
    int iflag=0;
    int nvecs;
    kt::mvec_num_vectors(mv.get(),&nvecs,&iflag);
    return nvecs;
  }

    static int GetVecLength( const BelosMV<ST>& mv )
    {
      return GetGlobalLength(mv);
    }

  static phist_gidx GetGlobalLength( const BelosMV<ST>& mv )
  {
    int iflag=0;
    phist_gidx nglob;
    kt::mvec_global_length(mv.get(),&nglob,&iflag);
    return nglob;
  }



  static Teuchos::RCP<BelosMV<ST> > Clone( const BelosMV<ST>& mv, const int numvecs )
  {
    int iflag=0;
    phist_const_map_ptr map=nullptr;
    kt::mvec_get_map(mv.get(),&map,&iflag);
    mvec_ptr new_mvec=nullptr;
    kt::mvec_create(&new_mvec,map,numvecs,&iflag);
    return mvec_rcp(new_mvec,true);
  }

static Teuchos::RCP<BelosMV<ST> > CloneCopy( const BelosMV<ST>& mv )
{
  int iflag=0;
  auto result = Clone(mv,GetNumberVecs(mv));
  kt::mvec_add_mvec(st::one(),mv.get(),st::zero(),result->get(),&iflag);
  return result;
}

    static Teuchos::RCP<BelosMV<ST> > CloneCopy( const BelosMV<ST>& mv, const std::vector<int>& index )
    {
      int iflag=0;
      int nvecs=(int)index.size();
      Teuchos::RCP<BelosMV<ST>> result = Clone(mv,nvecs);

      bool contig=true;
      for (int j=1; j<nvecs; j++) contig&=(index[j]==index[j-1]+1);

      if (contig)
      {
        int jmin = index[0];
        int jmax = index[nvecs-1];
        kt::mvec_get_block(mv.get(),result->get(),jmin,jmax,&iflag);
      }
      else
      {
        mvec_ptr result_j=nullptr;
        // do it one by one
        for (int j=0; j<nvecs; j++)
        {
          int jj=index[j];
          kt::mvec_view_block(result->get(),&result_j,j,j,&iflag);
          kt::mvec_get_block(mv.get(),result_j,jj,jj,&iflag); 
        }
        if (result_j!=nullptr) kt::mvec_delete(result_j,&iflag);
      }
      return result;
    }

    static Teuchos::RCP<BelosMV<ST> > 
    CloneCopy (const BelosMV<ST>& mv, 
         const Teuchos::Range1D& index)
    {
      int iflag=0;
      int jmin=index.lbound();
      int jmax=index.ubound();
      int nvec=jmax-jmin+1;
      auto result=Clone(mv,nvec);
      kt::mvec_get_block(mv.get(),result->get(),jmin,jmax,&iflag);
    return result;
    }


static Teuchos::RCP<BelosMV<ST> > CloneViewNonConst( BelosMV<ST>& mv, const std::vector<int>& index )
{
  TEUCHOS_TEST_FOR_EXCEPTION(index.size() == 0,std::invalid_argument,
          "Belos::MultiVecTraits<ST,BelosMV<ST>>::CloneView(mv,index): numvecs must be greater than zero.");
  int nvec=(int)index.size();
  bool contig=true;
  for (int j=2; j<nvec; j++) contig&=(index[j] == index[j-1]+1);
        
  TEUCHOS_TEST_FOR_EXCEPTION(!contig,std::invalid_argument,
          "Belos::MultiVecTraits<ST,BelosMV<ST>>::CloneView(mv,index): can only view contiguous blocks with PHIST.");
    int iflag=0;
    mvec_ptr result=NULL;
    kt::mvec_view_block(mv.get(),&result,index[0],index[nvec-1],&iflag);
      return mvec_rcp(result,true);
  }


    static Teuchos::RCP<BelosMV<ST> > 
    CloneViewNonConst (BelosMV<ST>& mv, 
           const Teuchos::Range1D& index)
    {
      int iflag=0;
      mvec_ptr result=nullptr;
      kt::mvec_view_block(mv.get(),&result,(int)index.lbound(),(int)index.ubound(),&iflag);
      return mvec_rcp(result,true);
    }


    static Teuchos::RCP<const BelosMV<ST> > CloneView(const BelosMV<ST>& mv, const std::vector<int>& index )
    {
      PHIST_ENTER_FCN(__FUNCTION__);    
      return Teuchos::rcp_dynamic_cast<const BelosMV<ST> >
        (CloneViewNonConst(const_cast<BelosMV<ST>&>(mv),index));
    }

    static Teuchos::RCP<const BelosMV<ST> > 
    CloneView (const BelosMV<ST>& mv, 
         const Teuchos::Range1D& index)
    {
      PHIST_ENTER_FCN(__FUNCTION__);    
      return Teuchos::rcp_dynamic_cast<const BelosMV<ST> >
        (CloneViewNonConst(const_cast<BelosMV<ST>&>(mv),index));
    }


    static bool HasConstantStride( const BelosMV<ST>& mv )
    { return true; }

    static void MvTimesMatAddMv( ST alpha, const BelosMV<ST>& A, 
                                 const Teuchos_sdMat& B, 
                                 ST beta, BelosMV<ST>& mv )
    {
      int iflag=0;
      // copy B to a phist sdMat
      sdMat_ptr Bphist=NULL;
      kt::sdMat_create(&Bphist,B.numRows(),B.numCols(),nullptr,&iflag);
      ST* B_raw=nullptr;
      ghost_lidx lda;
      kt::sdMat_extract_view(Bphist,&B_raw,&lda,&iflag);
      for (int j=0; j<B.numCols(); j++)
      {
        for (int i=0; i<B.numRows(); i++)
        {
          B_raw[j*lda+i]=B(i,j);
        }
      }

      kt::mvec_times_sdMat(alpha,A.get(),Bphist,beta,mv.get(),&iflag);
      kt::sdMat_delete(Bphist,&iflag);
    }

    // compute mv = alpha*A + beta*B. This function is abused in Belos by aliasing mv and A 
    // or B, e.g. A = 0*A+1*B instead of A=B. We therefore have to be a bit careful with
    // the memcopy we use if either alpha or beta are 0.
    static void MvAddMv( ST alpha, const BelosMV<ST>& A, ST beta, const BelosMV<ST>& B, BelosMV<ST>& mv )
    {
      int iflag=0;
      if (alpha==st::zero() && beta==st::zero())
      {
        kt::mvec_put_value(mv.get(),st::zero(),&iflag);
      }
      if (mv.get()==A.get())
      {
        kt::mvec_add_mvec(beta,B.get(),alpha,mv.get(),&iflag);
      }
      else if (mv.get()==B.get())
      {
        kt::mvec_add_mvec(alpha,A.get(),beta,mv.get(),&iflag);
      }
      else
      {
        ST gamma=st::zero();
        if (alpha!=st::zero())
        {
          kt::mvec_add_mvec(alpha,A.get(),gamma,mv.get(),&iflag);
          gamma=st::one();
        }
        if (beta!=st::zero())
        {
          kt::mvec_add_mvec(beta,B.get(),gamma,mv.get(),&iflag);
        }
      }
    }

    static void MvScale ( BelosMV<ST>& mv, ST alpha )
    {
      int iflag=0;
      kt::mvec_scale(mv.get(),alpha,&iflag);
    }

    static void MvScale ( BelosMV<ST>& mv, const std::vector<ST>& alphas )
    {
      int iflag=0;
      kt::mvec_vscale(mv.get(),&alphas[0],&iflag);
    }

    // C=alpha*A*B
    static void MvTransMv( ST alpha, const BelosMV<ST>& A, const BelosMV<ST>& B, Teuchos_sdMat& C)
    {
      int iflag=0;
      sdMat_ptr C_tmp=NULL;
      int nrC=GetNumberVecs(A),ncC=GetNumberVecs(B);
      phist_const_comm_ptr comm;
      kt::mvec_get_comm(A.get(),&comm,&iflag);
      kt::sdMat_create(&C_tmp,nrC,ncC,comm,&iflag);
      
      kt::mvecT_times_mvec(alpha,A.get(),B.get(),st::one(),C_tmp,&iflag);

      ST* C_raw=NULL;
      phist_lidx lda;
      kt::sdMat_extract_view(C_tmp,&C_raw,&lda,&iflag);
      
      for (int j=0; j<C.numCols(); j++)
        for (int i=0; i<C.numRows(); i++)
        {
          C(i,j)=C_raw[j*lda+i];
        }
      kt::sdMat_delete(C_tmp,&iflag);
    }

    static void MvDot( const BelosMV<ST>& A, const BelosMV<ST>& B, std::vector<ST> &dots)
    {
      int iflag=0;
      kt::mvec_dot_mvec(A.get(),B.get(),&dots[0],&iflag);          
    }

    static void MvNorm(const BelosMV<ST>& mv, std::vector<magn_t> &normvec, NormType type=TwoNorm)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(type != TwoNorm,std::invalid_argument,
          "Belos::MultiVecTraits<ST,BelosMV<ST>>::MvNorm(mv,normvec): MvNorm only accepts TwoNorm up to now.");
      int iflag=0;
      kt::mvec_norm2(mv.get(),&normvec[0],&iflag);          
    }

    static void SetBlock( const BelosMV<ST>& A, const std::vector<int>& index, BelosMV<ST>& mv )
    {
      int iflag=0;
      int nvecs=(int)index.size();

      bool contig=true;
      for (int j=1; j<nvecs; ++j) contig&=(index[j]==index[j-1]+1);

      if (contig)
      {
        int jmin = index[0];
        int jmax = index[nvecs-1];
        kt::mvec_set_block(mv.get(),A.get(),jmin,jmax,&iflag);
      }
      else
      {
        mvec_ptr result_j=nullptr, input_j=nullptr;
        // do it one by one
        for (int j=0; j<nvecs; j++)
        {
          int jj=index[j];
          kt::mvec_view_block(mv.get(),&result_j,jj,jj,&iflag);
          kt::mvec_view_block((mvec_ptr)A.get(),&input_j,j,j,&iflag);
          kt::mvec_set_block(result_j,input_j,0,0,&iflag); 
        }
        if (result_j!=nullptr) kt::mvec_delete(result_j,&iflag);
        if (input_j!=nullptr) kt::mvec_delete(input_j,&iflag);
      }
    }

    static void
    SetBlock (const BelosMV<ST>& A, 
        const Teuchos::Range1D& index, 
        BelosMV<ST>& mv)
    {
      int iflag=0;
      kt::mvec_set_block(mv.get(),A.get(),(int)index.lbound(),(int)index.ubound(),&iflag);
    }

    static void
    Assign (const BelosMV<ST>& A, 
      BelosMV<ST>& mv)
    {
      int nvA=GetNumberVecs(A);
      int nvmv=GetNumberVecs(mv);

  TEUCHOS_TEST_FOR_EXCEPTION(nvA!=nvmv,std::invalid_argument,
          "Belos::MultiVecTraits<ST,BelosMV<ST>>::Assign: mvecs must have same number of columns.");

      int iflag=0;
      kt::mvec_get_block(A.get(),mv.get(),0,nvA-1,&iflag);
    }


    static void MvRandom( BelosMV<ST>& mv )
    { 
      int iflag=0;
      kt::mvec_random(mv.get(),&iflag);
    }

    static void MvInit( BelosMV<ST>& mv, ST alpha = st::zero() )
    {
      int iflag=0;
      kt::mvec_put_value(mv.get(),alpha,&iflag);
    }

    static void MvPrint( const BelosMV<ST>& mv, std::ostream& os )
    {
      // TODO - the stream argument is ignored, ghost always prints to stdout
      int iflag=0;
      kt::mvec_print(mv.get(),&iflag);
    }


};

} // end of Belos namespace 

#endif // PHIST_HAVE_BELOS

#endif // end of file BELOS_GHOST_ADAPTER_HPP

