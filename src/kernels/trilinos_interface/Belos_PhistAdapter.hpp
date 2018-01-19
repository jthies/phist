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
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#ifdef PHIST_HAVE_BELOS

#include "phist_types.hpp"
#include "phist_kernels.hpp"

#include <ghost.h>
#include "phist_config.h"
#include "phist_typedefs.h"
#include "phist_kernels.h"
#include "phist_ScalarTraits.hpp"
#include "./typedefs.hpp"
#include "phist_tools.h"
#include "phist_macros.h"
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Array.hpp>

#include <BelosConfigDefs.hpp>
#include <BelosTypes.hpp>
#include <BelosMultiVecTraits.hpp>

#include "phist_GhostMV.hpp"
#include "phist_rcp_helpers.hpp"

#include <ghost/densemat.h>
#include <ghost/densemat_rm.h>
#include <ghost/densemat_cm.h>

// this file is mostly copied from the Belos Tpetra adapter implementation in Trilinos 11.2.4

#ifndef CHK_GERR
#define CHK_GERR(CALL,RETURNVALUE) \
  { \
    if (GHOST_SUCCESS!=CALL) \
    { \
      PHIST_SOUT(PHIST_ERROR,"ghost call %s failed (%s,%d)",#CALL,__FILE__,__LINE__);\
      return RETURNVALUE;\
    }\
  }
#endif



namespace Belos {

using ::phist::GhostMV;

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::MultiVecTraits for ghost_densemat.
  //
  ////////////////////////////////////////////////////////////////////

  /*!  \brief Template specialization of Belos::MultiVecTraits class using the ghost_densemat class.

    This interface will ensure that any ghost_densemat will be accepted by the Belos
    templated solvers.  
    
    NOTE: the implementation here is a bit outdated, with our cxx_bindings now it could be done much nicer using 
    templates directly instead of all the ifs. Also, there is a new mechanism in Belos to provide a multivector class
    via a factory, that may also be a good idea.
    */
  template<class Scalar>
  class MultiVecTraits<Scalar, GhostMV >
  {
  public:

    typedef ::phist::ScalarTraits<Scalar> st;
    typedef ::phist::types<Scalar> tt;
    typedef ::phist::kernels<Scalar> kt;
    
    typedef typename tt::mvec_ptr mvec_ptr;
    typedef typename tt::sdMat_ptr sdMat_ptr;
    
    typedef typename st::magn_t magn_t;
    typedef typename Traits<Scalar>::Teuchos_sdMat_t Teuchos_sdMat_t;

  static int GetNumberVecs( const GhostMV& mv )
  {
    int iflag=0;
    int nvecs;
    kt::mvec_num_vectors(mv.get(),&nvecs,&iflag);
    return nvecs;
  }

    static int GetVecLength( const GhostMV& mv )
    {
      return GetGlobalLength(mv);
    }

  static phist_gidx GetGlobalLength( const GhostMV& mv )
  {
    int iflag=0;
    phist_gidx nglob;
    kt::mvec_global_length(mv.get(),&nglob,&iflag);
    return nglob;
  }



  static Teuchos::RCP<GhostMV > Clone( const GhostMV& mv, const int numvecs )
  {
    int iflag=0;
    phist_const_map_ptr map=nullptr;
    kt::mvec_get_map(mv.get(),&map,&iflag);
    mvec_ptr new_mvec=nullptr;
    kt::mvec_create(&new_mvec,map,numvecs,&iflag);
    return phist::rcp((ghost_densemat*)new_mvec,true);
  }

static Teuchos::RCP<GhostMV > CloneCopy( const GhostMV& mv )
{
  int iflag=0;
  auto result = Clone(mv,GetNumberVecs(mv));
  kt::mvec_add_mvec(st::one(),mv.get(),st::zero(),result->get(),&iflag);
  return result;
}

    static Teuchos::RCP<GhostMV > CloneCopy( const GhostMV& mv, const std::vector<int>& index )
    {
      int iflag=0;
      int nvecs=(int)index.size();
      Teuchos::RCP<GhostMV> result = Clone(mv,nvecs);

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

    static Teuchos::RCP<GhostMV > 
    CloneCopy (const GhostMV& mv, 
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


static Teuchos::RCP<GhostMV > CloneViewNonConst( GhostMV& mv, const std::vector<int>& index )
{
  TEUCHOS_TEST_FOR_EXCEPTION(index.size() == 0,std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,GhostMV>::CloneView(mv,index): numvecs must be greater than zero.");
  int nvec=(int)index.size();
  bool contig=true;
  for (int j=2; j<nvec; j++) contig&=(index[j] == index[j-1]+1);
        
  TEUCHOS_TEST_FOR_EXCEPTION(!contig,std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,GhostMV>::CloneView(mv,index): can only view contiguous blocks with PHIST.");
    int iflag=0;
    mvec_ptr result=NULL;
    kt::mvec_view_block(mv.get(),&result,index[0],index[nvec-1],&iflag);
      return phist::rcp((ghost_densemat*)result,true);
  }


    static Teuchos::RCP<GhostMV > 
    CloneViewNonConst (GhostMV& mv, 
           const Teuchos::Range1D& index)
    {
      int iflag=0;
      mvec_ptr result=nullptr;
      kt::mvec_view_block(mv.get(),&result,(int)index.lbound(),(int)index.ubound(),&iflag);
      return phist::rcp((ghost_densemat*)result,true);
    }


    static Teuchos::RCP<const GhostMV > CloneView(const GhostMV& mv, const std::vector<int>& index )
    {
      PHIST_ENTER_FCN(__FUNCTION__);    
      return Teuchos::rcp_dynamic_cast<const GhostMV >
        (CloneViewNonConst(const_cast<GhostMV&>(mv),index));
    }

    static Teuchos::RCP<const GhostMV > 
    CloneView (const GhostMV& mv, 
         const Teuchos::Range1D& index)
    {
      PHIST_ENTER_FCN(__FUNCTION__);    
      return Teuchos::rcp_dynamic_cast<const GhostMV >
        (CloneViewNonConst(const_cast<GhostMV&>(mv),index));
    }


    static bool HasConstantStride( const GhostMV& mv )
    { return true; }

    static void MvTimesMatAddMv( Scalar alpha, const GhostMV& A, 
                                 const Teuchos_sdMat_t& B, 
                                 Scalar beta, GhostMV& mv )
    {
      ghost_densemat* Bphist=createGhostCopyOfTeuchosSDM(B);
      int iflag=0;
      kt::mvec_times_sdMat(alpha,A.get(),Bphist,beta,mv.get(),&iflag);
      kt::mvec_delete(Bphist,&iflag);
    }

    // compute mv = alpha*A + beta*B. This function is abused in Belos by aliasing mv and A 
    // or B, e.g. A = 0*A+1*B instead of A=B. We therefore have to be a bit careful with
    // the memcopy we use if either alpha or beta are 0.
    static void MvAddMv( Scalar alpha, const GhostMV& A, Scalar beta, const GhostMV& B, GhostMV& mv )
    {
      int iflag=0;
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
        kt::mvec_add_mvec(alpha,A.get(),st::zero(),mv.get(),&iflag);
        kt::mvec_add_mvec(beta,B.get(),st::one(),mv.get(),&iflag);
      }
    }

    static void MvScale ( GhostMV& mv, Scalar alpha )
    {
      int iflag=0;
      kt::mvec_scale(mv.get(),alpha,&iflag);
    }

    static void MvScale ( GhostMV& mv, const std::vector<Scalar>& alphas )
    {
      int iflag=0;
      kt::mvec_vscale(mv.get(),&alphas[0],&iflag);
    }

    // C=alpha*A*B
    static void MvTransMv( Scalar alpha, const GhostMV& A, const GhostMV& B, Teuchos_sdMat_t& C)
    {
      int iflag=0;
      sdMat_ptr C_tmp=NULL;
      int nrC=GetNumberVecs(A),ncC=GetNumberVecs(B);
      phist_const_comm_ptr comm;
      kt::mvec_get_comm(A.get(),&comm,&iflag);
      kt::sdMat_create(&C_tmp,nrC,ncC,comm,&iflag);
      
      kt::mvecT_times_mvec(alpha,A.get(),B.get(),st::one(),C_tmp,&iflag);

      Scalar* C_raw=NULL;
      phist_lidx lda;
      kt::sdMat_extract_view(C_tmp,&C_raw,&lda,&iflag);
      
      for (int j=0; j<C.numCols(); j++)
        for (int i=0; i<C.numRows(); i++)
        {
          C(i,j)=C_raw[j*lda+i];
        }
      kt::sdMat_delete(C_tmp,&iflag);
    }

    static void MvDot( const GhostMV& A, const GhostMV& B, std::vector<Scalar> &dots)
    {
      int iflag=0;
      kt::mvec_dot_mvec(A.get(),B.get(),&dots[0],&iflag);          
    }

    static void MvNorm(const GhostMV& mv, std::vector<magn_t> &normvec, NormType type=TwoNorm)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(type != TwoNorm,std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,GhostMV>::MvNorm(mv,normvec): MvNorm only accepts TwoNorm up to now.");
      int iflag=0;
      kt::mvec_norm2(mv.get(),&normvec[0],&iflag);          
    }

    static void SetBlock( const GhostMV& A, const std::vector<int>& index, GhostMV& mv )
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
    SetBlock (const GhostMV& A, 
        const Teuchos::Range1D& index, 
        GhostMV& mv)
    {
      int iflag=0;
      kt::mvec_set_block(mv.get(),A.get(),(int)index.lbound(),(int)index.ubound(),&iflag);
    }

    static void
    Assign (const GhostMV& A, 
      GhostMV& mv)
    {
      int nvA=GetNumberVecs(A);
      int nvmv=GetNumberVecs(mv);

  TEUCHOS_TEST_FOR_EXCEPTION(nvA!=nvmv,std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,GhostMV>::Assign: mvecs must have same number of columns.");

      int iflag=0;
      kt::mvec_get_block(A.get(),mv.get(),0,nvA-1,&iflag);
    }


    static void MvRandom( GhostMV& mv )
    { 
      int iflag=0;
      kt::mvec_random(mv.get(),&iflag);
    }

    static void MvInit( GhostMV& mv, Scalar alpha = st::zero() )
    {
      int iflag=0;
      kt::mvec_put_value(mv.get(),alpha,&iflag);
    }

    static void MvPrint( const GhostMV& mv, std::ostream& os )
    {
      // TODO - the stream argument is ignored, ghost always prints to stdout
      ghost_densemat* _mv = const_cast<ghost_densemat*>(mv.get());
      char* the_string=nullptr;
      ghost_densemat_string(&the_string,_mv);
      os << the_string << std::endl;
      delete [] the_string;
    }

  // private helper function
  static ghost_densemat* createGhostCopyOfTeuchosSDM
        (const Teuchos_sdMat_t& M)
  {
    PHIST_ENTER_FCN(__FUNCTION__);
      int nrows=M.numRows();
      int ncols=M.numCols();

      ghost_densemat* Mghost=nullptr;
      int iflag=0;      
#ifdef PHIST_HAVE_SP
      if (st::type_char()=='S')
      {
        phist_SsdMat_create((void**)&Mghost,nrows,ncols,nullptr,&iflag);
      }
      else if (st::type_char()=='C')
      {
        phist_CsdMat_create((void**)&Mghost,nrows,ncols,nullptr,&iflag);
      }
      else 
#endif
      if (st::type_char()=='D')
      {
        phist_DsdMat_create((void**)&Mghost,nrows,ncols,nullptr,&iflag);
      }
      else if (st::type_char()=='Z')
      {
        phist_ZsdMat_create((void**)&Mghost,nrows,ncols,nullptr,&iflag);
      }
    if (iflag!=PHIST_SUCCESS)
    {
      PHIST_SOUT(PHIST_ERROR,"phist_XsdMat_create_view returned non-zero error code %d\n"
                             "(file %s, line %d)\n",iflag,__FILE__,__LINE__);
    }
    // copy data from Teuchos matrix
    if (Mghost->traits.storage!=GHOST_DENSEMAT_COLMAJOR)
    {
       throw "sdMat not col-major!";
    }
    Scalar* M_raw=(Scalar*)Mghost->val;
    ghost_lidx lda=Mghost->stride;
    for (int j=0; j<M.numCols(); j++)
        for (int i=0; i<M.numRows(); i++)
        {
          M_raw[j*lda+i]=M(i,j);
        }
    return Mghost;
  }

};

} // end of Belos namespace 

#endif // PHIST_HAVE_BELOS

#endif // end of file BELOS_GHOST_ADAPTER_HPP

