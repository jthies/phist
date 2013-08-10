#ifndef BELOS_GHOST_ADAPTER_HPP
#define BELOS_GHOST_ADAPTER_HPP

#include "ghost.h"

#include <phist_ScalarTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Array.hpp>

#include <BelosConfigDefs.hpp>
#include <BelosTypes.hpp>
#include <BelosMultiVecTraits.hpp>
#include <BelosOperatorTraits.hpp>

#ifdef HAVE_BELOS_TSQR
#  include <Ghost_TsqrAdaptor.hpp>
#endif // HAVE_BELOS_TSQR

// this file is mostly copied from the Belos Tpetra adapter implementation in Trilinos 11.2.4

namespace Belos {

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::MultiVecTraits for ghost_vec_t.
  //
  ////////////////////////////////////////////////////////////////////

  /*!  \brief Template specialization of Belos::MultiVecTraits class using the ghost_vec_t class.

    This interface will ensure that any ghost_vec_t will be accepted by the Belos
    templated solvers.  */
  template<class Scalar>
  class MultiVecTraits<Scalar, ghost_vec_t >
  {
  public:

    static Teuchos::RCP<ghost_vec_t > Clone( const ghost_vec_t& mv, const int numvecs )
    { 
      // TODO - this actually copies the data as well
      return Teuchos::rcp( mv.extract(const_cast<ghost_vec_t*>(&mv),0,nvecs) );
    }

    static Teuchos::RCP<ghost_vec_t > CloneCopy( const ghost_vec_t& mv )
    {
      return Teuchos::rcp( mv.extract(const_cast<ghost_vec_t*>(&mv),0,nvecs) );
    }

    static Teuchos::RCP<ghost_vec_t > CloneCopy( const ghost_vec_t& mv, const std::vector<int>& index )
    { 
      TEUCHOS_TEST_FOR_EXCEPTION(index.size() == 0,std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,ghost_vec_t>::CloneCopy(mv,index): numvecs must be greater than zero.");

      TEUCHOS_TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::runtime_error,
          "Belos::MultiVecTraits<Scalar,ghost_vec_t>::CloneCopy(mv,index): indices must be >= zero.");

      TEUCHOS_TEST_FOR_EXCEPTION( (size_t)*std::max_element(index.begin(),index.end()) >= mv.getNumVectors(), std::runtime_error,
          "Belos::MultiVecTraits<Scalar,ghost_vec_t>::CloneCopy(mv,index): indices must be < mv.getNumVectors().");

      bool contig=true;
      for (typename std::vector<int>::size_type j=1; j<index.size(); ++j) {
        if (index[j] != index[j-1]+1) {
          // not contiguous; short circuit
          contig=false;
          break;
        }
      }
    
    TEUCHOS_TEST_FOR_EXCEPTION(contig==false, std::runtime_error,
            "Belos::MultiVecTraits<Scalar,ghost_vec_t>::CloneCopy(mv,index): only contiguous range implemented so far.");

      // contiguous
      return Teuchos::rcp( mv.extract(const_cast<ghost_vec_t*>(&mv),index[0],index.size()) );
    }

    static Teuchos::RCP<ghost_vec_t > 
    CloneCopy (const ghost_vec_t& mv, 
	       const Teuchos::Range1D& index)
    {
      const bool validRange = index.size() > 0 && 
	index.lbound() >= 0 && 
	index.ubound() < GetNumberVecs(mv);
      if (! validRange)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<Scalar, ghost_vec_t<...> >::"
	    "CloneCopy(mv,index=[" << index.lbound() << ", " << index.ubound() 
	     << "]): ";
	  TEUCHOS_TEST_FOR_EXCEPTION(index.size() == 0, std::invalid_argument,
			     os.str() << "Empty index range is not allowed.");
	  TEUCHOS_TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
			     os.str() << "Index range includes negative "
			     "index/ices, which is not allowed.");
	  // Range1D bounds are signed; size_t is unsigned.
	  TEUCHOS_TEST_FOR_EXCEPTION(index.ubound() >= GetNumberVecs(mv),
			     std::invalid_argument, 
			     os.str() << "Index range exceeds number of vectors " 
			     << GetNumberVecs() << " in the input multivector.");
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			     os.str() << "Should never get here!");
	}
      return Teuchos::rcp( mv.extract(const_cast<ghost_vec_t*>(&mv),index.lbound(),index.ubound()-index.lbound()+1) );
    }


    static Teuchos::RCP<ghost_vec_t > CloneViewNonConst( ghost_vec_t& mv, const std::vector<int>& index )
    {
      TEUCHOS_TEST_FOR_EXCEPTION(index.size() == 0,std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,ghost_vec_t>::CloneView(mv,index): numvecs must be greater than zero.");
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,ghost_vec_t>::CloneView(mv,index): indices must be >= zero.");
      TEUCHOS_TEST_FOR_EXCEPTION( (size_t)*std::max_element(index.begin(),index.end()) >= mv.getNumVectors(), std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,ghost_vec_t>::CloneView(mv,index): indices must be < mv.getNumVectors().");
#endif
      bool contig=true;
      for (typename std::vector<int>::size_type j=1; j<index.size(); ++j) {
        if (index[j] != index[j-1]+1) {
          contig=false;
          break;
        }
      }

    TEUCHOS_TEST_FOR_EXCEPTION(contig==false, std::runtime_error,
            "Belos::MultiVecTraits<Scalar,ghost_vec_t>::CloneCopy(mv,index): only contiguous range implemented so far.");

      // contiguous
      return Teuchos::rcp( mv.view(&mv,index[0],index.size()) );
      }


    static Teuchos::RCP<ghost_vec_t > 
    CloneViewNonConst (ghost_vec_t& mv, 
		       const Teuchos::Range1D& index)
    {
      // NOTE (mfh 11 Jan 2011) We really should check for possible
      // overflow of int here.  However, the number of columns in a
      // multivector typically fits in an int.
      const int numCols = static_cast<int> (mv.getNumVectors());
      const bool validRange = index.size() > 0 && 
	index.lbound() >= 0 && index.ubound() < numCols;
      if (! validRange)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<Scalar, ghost_vec_t<...> >::"
	    "CloneViewNonConst(mv,index=[" << index.lbound() << ", " 
	     << index.ubound() << "]): ";
	  TEUCHOS_TEST_FOR_EXCEPTION(index.size() == 0, std::invalid_argument,
			     os.str() << "Empty index range is not allowed.");
	  TEUCHOS_TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
			     os.str() << "Index range includes negative "
			     "index/ices, which is not allowed.");
	  TEUCHOS_TEST_FOR_EXCEPTION(index.ubound() >= numCols, std::invalid_argument, 
			     os.str() << "Index range exceeds number of "
			     "vectors " << numCols << " in the input "
			     "multivector.");
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			     os.str() << "Should never get here!");
	}
      return Teuchos::rcp( mv.view(&mv,index.lbound(),index.ubound()-index.lbound()+1) );
    }


    static Teuchos::RCP<const ghost_vec_t > CloneView(const ghost_vec_t& mv, const std::vector<int>& index )
    {
    return Teuchos::rcp_dynamic_cast<const ghost_vec_t >(CloneViewNonConst      
        (std::const_cast<ghost_vec_t&>(mv),index);
    }

    static Teuchos::RCP<const ghost_vec_t > 
    CloneView (const ghost_vec_t& mv, 
	       const Teuchos::Range1D& index)
    {
    return Teuchos::rcp_dynamic_cast<const ghost_vec_t >(CloneViewNonConst      
        (std::const_cast<ghost_vec_t&>(mv),index));
    }

    static int GetVecLength( const ghost_vec_t& mv )
    { 
    return mv.traits->TODO;
    }

    static int GetNumberVecs( const ghost_vec_t& mv )
    { 
    return mv.traits->nvecs;
    }

    static bool HasConstantStride( const ghost_vec_t& mv )
    { return true; }

    static void MvTimesMatAddMv( Scalar alpha, const ghost_vec_t& A, 
                                 const Teuchos::SerialDenseMatrix<int,Scalar>& B, 
                                 Scalar beta, ghost_vec_t& mv )
    {
      // create local map
      ghost_vtraits_t dmtraits = GHOST_VTRAITS_INIT
        (.flags = GHOST_VEC_LHS, .nrows=B.numRows(), .nvecs=B.numCols(), .datatype=phist::ScalarTraits<scalar_type>::ghost_dt);
      // create view of Teuchos matrix as ghost_vec_t
      // TODO: ghost doesn't have this construction mechanism yet
      // B.values(),B.stride(),B.numRows(),B.numCols()
      // multiply
      // TODO - ghost tries to allocate the result vector, which is undesired here.
      ghost_gemm("T",&A,&B,&mv,&alpha,&beta,GHOST_GEMM_NO_REDUCE);
              
    }

    static void MvAddMv( Scalar alpha, const ghost_vec_t& A, Scalar beta, const ghost_vec_t& B, ghost_vec_t& mv )
    { 
      //TODO - do this in one step if alpha and beta are not 0
      mv.zero(&mv);
      if (alpha==Teuchos::ScalarTraits<Scalar>::zero())
        {
        mv.axpy(&mv,A,(void*)&alpha);
        }
      if (beta==Teuchos::ScalarTraits<Scalar>::zero())
        {
        mv.axpy(&mv,B,(void*)&beta);
        }
    }

    static void MvScale ( ghost_vec_t& mv, Scalar alpha )
    { mv.scale(&mv,(void*)&alpha); }

    static void MvScale ( ghost_vec_t& mv, const std::vector<Scalar>& alphas )
    { 
    TODO
    }

    static void MvTransMv( Scalar alpha, const ghost_vec_t& A, const ghost_vec_t& B, Teuchos::SerialDenseMatrix<int,Scalar>& C)
    {
    TODO 
      // form alpha * A^H * B, then copy into SDM
      // we will create a multivector C_mv from a a local map
      // this map has a serial comm, the purpose being to short-circuit the MultiVector::reduce() call at the end of MultiVector::multiply()
      // otherwise, the reduced multivector data would be copied back to the GPU, only to turn around and have to get it back here.
      // this saves us a round trip for this data.
/*
      const int numRowsC = C.numRows(),
                numColsC = C.numCols(),
                strideC  = C.stride();
      Teuchos::SerialComm<int> scomm;
      // create local map with serial comm
      Tpetra::Map<LO,GO,Node> LocalMap(numRowsC, 0, Teuchos::rcpFromRef< const Teuchos::Comm<int> >(scomm), Tpetra::LocallyReplicated, A.getMap()->getNode());
      // create local multivector to hold the result
      const bool INIT_TO_ZERO = true;
      ghost_vec_t C_mv(Teuchos::rcpFromRef(LocalMap),numColsC, INIT_TO_ZERO);
      // multiply result into local multivector
      C_mv.multiply(Teuchos::CONJ_TRANS,Teuchos::NO_TRANS,alpha,A,B,Teuchos::ScalarTraits<Scalar>::zero());
      // get comm
      Teuchos::RCP< const Teuchos::Comm<int> > pcomm = A.getMap()->getComm();
      // create arrayview encapsulating the Teuchos::SerialDenseMatrix
      Teuchos::ArrayView<Scalar> C_view(C.values(),strideC*numColsC);
      if (pcomm->getSize() == 1) {
        // no accumulation to do; simply extract the multivector data into C
        // extract a copy of the result into the array view (and therefore, the SerialDenseMatrix)
        C_mv.get1dCopy(C_view,strideC);
      }  
      else {
        // get a const host view of the data in C_mv
        Teuchos::ArrayRCP<const Scalar> C_mv_view = C_mv.get1dView();
        if (strideC == numRowsC) {
          // sumall into C
          Teuchos::reduceAll<int,Scalar>(*pcomm,Teuchos::REDUCE_SUM,numColsC*numRowsC,C_mv_view.getRawPtr(),C_view.getRawPtr());
        }
        else {
          // sumall into temp, copy into C
          Teuchos::Array<Scalar> destBuff(numColsC*numRowsC);
          Teuchos::reduceAll<int,Scalar>(*pcomm,Teuchos::REDUCE_SUM,numColsC*numRowsC,C_mv_view.getRawPtr(),destBuff.getRawPtr());
          for (int j=0; j < numColsC; ++j) {
            for (int i=0; i < numRowsC; ++i) {
              C_view[strideC*j+i] = destBuff[numRowsC*j+i];
            }
          }
        }
      }
    */
    }

    static void MvDot( const ghost_vec_t& A, const ghost_vec_t& B, std::vector<Scalar> &dots)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(A.getNumVectors() != B.getNumVectors(),std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,ghost_vec_t>::MvDot(A,B,dots): A and B must have the same number of vectors.");

      TEUCHOS_TEST_FOR_EXCEPTION(dots.size() < (typename std::vector<int>::size_type)A.traits.nvecs,std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,ghost_vec_t>::MvDot(A,B,dots): dots must have room for all dot products.");

      Teuchos::ArrayView<Scalar> av(dots)
      A.dotProduct(const_cast<ghost_vec_t*>(&A),std::const_cast<ghost_vec_t*>(&B),dots.getRawPtr());
    }

    static void MvNorm(const ghost_vec_t& mv, std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &normvec, NormType type=TwoNorm)
    { 
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(normvec.size() < (typename std::vector<int>::size_type)mv.traits.nvecs,std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,ghost_vec_t>::MvNorm(mv,normvec): normvec must have room for all norms.");
#endif
      Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> av(normvec);
      TEUCHOS_TEST_FOR_EXCEPTION(type != TwoNorm,std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,ghost_vec_t>::MvNorm(mv,normvec): normvec only accepts TwoNorm up to now.");
      
      switch (type) {
        case OneNorm:
          break;
        case TwoNorm:
          // TODO - ghost doesn't support the dotProduct (and others) for multiple vecs.
          // TODO - we can't put in a magnitude type here
          mv.dotProduct(std::const_cast<ghost_vec_t*>(&mv),
                        std::const_cast<ghost_vec_t*>(&mv),
                        av.getRawPtr());
          for (int i=0;i<av.size();i++)
            {
            av[i]=std::sqrt(av[i]);
            }
          break;
        case InfNorm:
          break;
      }
    }

TODO - haven't looked at the rest, yet

    static void SetBlock( const ghost_vec_t& A, const std::vector<int>& index, ghost_vec_t& mv )
    {
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION((typename std::vector<int>::size_type)A.getNumVectors() < index.size(),std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,ghost_vec_t>::SetBlock(A,index,mv): index must be the same size as A.");
#endif
      Teuchos::RCP<ghost_vec_t > mvsub = CloneViewNonConst(mv,index);
      if ((typename std::vector<int>::size_type)A.getNumVectors() > index.size()) {
        Teuchos::RCP<const ghost_vec_t > Asub = A.subView(Teuchos::Range1D(0,index.size()-1));
        (*mvsub) = (*Asub);
      }
      else {
        (*mvsub) = A;
      }
      mvsub = Teuchos::null;
    }

    static void
    SetBlock (const ghost_vec_t& A, 
	      const Teuchos::Range1D& index, 
	      ghost_vec_t& mv)
    {

      // Range1D bounds are signed; size_t is unsigned.
      // Assignment of ghost_vec_t is a deep copy.

      // ghost_vec_t::getNumVectors() returns size_t.  It's
      // fair to assume that the number of vectors won't overflow int,
      // since the typical use case of multivectors involves few
      // columns, but it's friendly to check just in case.
      const size_t maxInt = static_cast<size_t> (Teuchos::OrdinalTraits<int>::max());
      const bool overflow = maxInt < A.getNumVectors() && maxInt < mv.getNumVectors();
      if (overflow)
	{
	  std::ostringstream os;
	  os <<	"Belos::MultiVecTraits<Scalar, ghost_vec_t<Scalar, ..."
	    "> >::SetBlock(A, index=[" << index.lbound() << ", " 
	     << index.ubound() << "], mv): ";
	  TEUCHOS_TEST_FOR_EXCEPTION(maxInt < A.getNumVectors(), std::range_error,
			     os.str() << "Number of columns in the input multi"
			     "vector 'A' (a size_t) overflows int.");
	  TEUCHOS_TEST_FOR_EXCEPTION(maxInt < mv.getNumVectors(), std::range_error,
			     os.str() << "Number of columns in the output multi"
			     "vector 'mv' (a size_t) overflows int.");
	}
      // We've already validated the static casts above.
      const int numColsA = static_cast<int> (A.getNumVectors());
      const int numColsMv = static_cast<int> (mv.getNumVectors());
      // 'index' indexes into mv; it's the index set of the target.
      const bool validIndex = index.lbound() >= 0 && index.ubound() < numColsMv;
      // We can't take more columns out of A than A has.
      const bool validSource = index.size() <= numColsA;

      if (! validIndex || ! validSource)
	{
	  std::ostringstream os;
	  os <<	"Belos::MultiVecTraits<Scalar, ghost_vec_t<Scalar, ..."
	    "> >::SetBlock(A, index=[" << index.lbound() << ", " 
	     << index.ubound() << "], mv): ";
	  TEUCHOS_TEST_FOR_EXCEPTION(index.lbound() < 0, std::invalid_argument,
			     os.str() << "Range lower bound must be nonnegative.");
	  TEUCHOS_TEST_FOR_EXCEPTION(index.ubound() >= numColsMv, std::invalid_argument,
			     os.str() << "Range upper bound must be less than "
			     "the number of columns " << numColsA << " in the "
			     "'mv' output argument.");
	  TEUCHOS_TEST_FOR_EXCEPTION(index.size() > numColsA, std::invalid_argument,
			     os.str() << "Range must have no more elements than"
			     " the number of columns " << numColsA << " in the "
			     "'A' input argument.");
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
	}
      typedef Teuchos::RCP<ghost_vec_t > MV_ptr;
      typedef Teuchos::RCP<const ghost_vec_t > const_MV_ptr;

      // View of the relevant column(s) of the target multivector mv.
      // We avoid view creation overhead by only creating a view if
      // the index range is different than [0, (# columns in mv) - 1].
      MV_ptr mv_view;
      if (index.lbound() == 0 && index.ubound()+1 == numColsMv)
	mv_view = Teuchos::rcpFromRef (mv); // Non-const, non-owning RCP
      else
	mv_view = CloneViewNonConst (mv, index);

      // View of the relevant column(s) of the source multivector A.
      // If A has fewer columns than mv_view, then create a view of
      // the first index.size() columns of A.
      const_MV_ptr A_view;
      if (index.size() == numColsA)
	A_view = Teuchos::rcpFromRef (A); // Const, non-owning RCP
      else
	A_view = CloneView (A, Teuchos::Range1D(0, index.size()-1));

      // Assignment of ghost_vec_t objects via operator=()
      // assumes that both arguments have compatible Maps.  If
      // HAVE_TPETRA_DEBUG is defined at compile time, operator=()
      // will throw an std::runtime_error if the Maps are
      // incompatible.
      *mv_view = *A_view; 
    }

    static void
    Assign (const ghost_vec_t& A, 
	    ghost_vec_t& mv)
    {

      // Range1D bounds are signed; size_t is unsigned.
      // Assignment of ghost_vec_t is a deep copy.

      // ghost_vec_t::getNumVectors() returns size_t.  It's
      // fair to assume that the number of vectors won't overflow int,
      // since the typical use case of multivectors involves few
      // columns, but it's friendly to check just in case.
      const size_t maxInt = static_cast<size_t> (Teuchos::OrdinalTraits<int>::max());
      const bool overflow = maxInt < A.getNumVectors() && maxInt < mv.getNumVectors();
      if (overflow)
	{
	  std::ostringstream os;
	  os <<	"Belos::MultiVecTraits<Scalar, ghost_vec_t<Scalar, ..."
	    "> >::Assign(A, mv): ";
	  TEUCHOS_TEST_FOR_EXCEPTION(maxInt < A.getNumVectors(), std::range_error,
			     os.str() << "Number of columns in the input multi"
			     "vector 'A' (a size_t) overflows int.");
	  TEUCHOS_TEST_FOR_EXCEPTION(maxInt < mv.getNumVectors(), std::range_error,
			     os.str() << "Number of columns in the output multi"
			     "vector 'mv' (a size_t) overflows int.");
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
	}
      // We've already validated the static casts above.
      const int numColsA = static_cast<int> (A.getNumVectors());
      const int numColsMv = static_cast<int> (mv.getNumVectors());
      if (numColsA > numColsMv)
	{
	  std::ostringstream os;
	  os <<	"Belos::MultiVecTraits<Scalar, ghost_vec_t<Scalar, ..."
	    "> >::Assign(A, mv): ";
	  TEUCHOS_TEST_FOR_EXCEPTION(numColsA > numColsMv, std::invalid_argument,
			     os.str() << "Input multivector 'A' has " 
			     << numColsA << " columns, but output multivector "
			     "'mv' has only " << numColsMv << " columns.");
	  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
	}
      // Assignment of ghost_vec_t objects via operator=()
      // assumes that both arguments have compatible Maps.  If
      // HAVE_TPETRA_DEBUG is defined at compile time, operator=()
      // will throw an std::runtime_error if the Maps are
      // incompatible.
      if (numColsA == numColsMv)
	mv = A;
      else
	{
	  Teuchos::RCP<ghost_vec_t > mv_view = 
	    CloneViewNonConst (mv, Teuchos::Range1D(0, numColsA-1));
	  *mv_view = A;
	}
    }


    static void MvRandom( ghost_vec_t& mv )
    { 
      mv.randomize(); 
    }

    static void MvInit( ghost_vec_t& mv, Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() )
    { mv.putScalar(alpha); }

    static void MvPrint( const ghost_vec_t& mv, std::ostream& os )
    { 
      Teuchos::FancyOStream fos(Teuchos::rcp(&os,false));
      mv.describe(fos,Teuchos::VERB_EXTREME);
    }

#ifdef HAVE_BELOS_TSQR
    /// \typedef tsqr_adaptor_type
    /// \brief TsqrAdaptor specialization for ghost_vec_t
    ///
    typedef GhostTsqrAdaptor<Scalar> tsqr_adaptor_type;
#endif // HAVE_BELOS_TSQR
  }; 

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::OperatorTraits for Tpetra::Operator.
  //
  ////////////////////////////////////////////////////////////////////

  /// \brief Partial specialization of OperatorTraits for Tpetra::Operator.
  /// Note: this is only going to work if MV is the class used by the kernel
  /// library, for instance, you can't compile phist with PHIST_KERNEL_LIB=
  /// tpetra and then use ghost vectors in Belos.
  template <class Scalar, class MV> 
  class OperatorTraits <Scalar, MV, phist::ScalarTraits<Scalar>::op_t >
  {
  public:

    typedef phist::ScalarTraits<Scalar>::op_t phist_op_t;
    typedef phist::ScalarTraits<Scalar>::mvec_t phist_mvec_t;

    static void 
    Apply (const phist_op_t& Op, 
	   const MV& X,
	   MV& Y,
	   ETrans trans=NOTRANS)
    {
    TEUCHOS_TEST_FOR_EXCEPTION(trans!=NOTRANS,std::invalid_argument,
          "Belos::OperatorTraits<Scalar,MV,phist_op_t>:: Apply: only implemented for trans=NOTRANS up to now.");
    int ierr;
    Scalar alpha = Teuchos::ScalarTraits<Scalar>::one();
    Scalar beta = Teuchos::ScalarTraits<Scalar>::zero();
    op.apply(alpha,(const phist_mvec_t*)&X, beta, (phist_mvec_t*)&Y,&ierr);
    }

    static bool
    HasApplyTranspose (const phist_op_t& Op)
    {
      return false;
    }
  };

} // end of Belos namespace 

#endif 
// end of file BELOS_GHOST_ADAPTER_HPP
