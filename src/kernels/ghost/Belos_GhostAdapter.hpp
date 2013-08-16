#ifndef BELOS_GHOST_ADAPTER_HPP
#define BELOS_GHOST_ADAPTER_HPP

#include "ghost.h"

#include <phist_ScalarTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Array.hpp>

#include <BelosConfigDefs.hpp>
#include <BelosTypes.hpp>
#include <BelosMultiVecTraits.hpp>

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

      TEUCHOS_TEST_FOR_EXCEPTION( (size_t)*std::max_element(index.begin(),index.end()) >= mv.traits.nvecs, std::runtime_error,
          "Belos::MultiVecTraits<Scalar,ghost_vec_t>::CloneCopy(mv,index): indices must be < mv.traits.nvecs.");

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
      TEUCHOS_TEST_FOR_EXCEPTION( (size_t)*std::max_element(index.begin(),index.end()) >= mv.traits.nvecs, std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,ghost_vec_t>::CloneView(mv,index): indices must be < mv.traits.nvecs.");
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
      const int numCols = static_cast<int> (mv.traits.nvecs);
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
    return Teuchos::rcp_dynamic_cast<const ghost_vec_t >
        (CloneViewNonConst(std::const_cast<ghost_vec_t&>(mv),index));
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
// TOOD - this is not available in ghost,
// but as it is only used for confirming 
// that vectors are compatible for certain
// operations in Belos, I think we can
// cheat here and simply return the local vec length instead.
    return mv.traits->nrows;
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
      // create view of Teuchos matrix as ghost_vec_t
      ghost_vec_t* Bghost=createGhostViewOfTeuchosSDM(B);
      // multiply
      ghost_gemm("T",&A,Bghost,&mv,&alpha,&beta,GHOST_GEMM_NO_REDUCE);
      
      ghost_freeVector(Bghost);
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

    // C=alpha*A*B
    static void MvTransMv( Scalar alpha, const ghost_vec_t& A, const ghost_vec_t& B, Teuchos::SerialDenseMatrix<int,Scalar>& C)
    {
    ghost_vec_t* Cghost=createGhostViewOfTeuchosSDM(C);
    ghost_gemm("T",std::const_cast<ghost_vec_t*>(&A),
                   std::const_cast<ghost_vec_t*>(&B),
                   Cghost,
                   &alpha, &beta,
                   GHOST_GEMM_ALLREDUCE);
    ghost_freeVector(Cghost);
    }

    static void MvDot( const ghost_vec_t& A, const ghost_vec_t& B, std::vector<Scalar> &dots)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(A.traits.nvecs != B.traits.nvecs,std::invalid_argument,
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
          "Belos::MultiVecTraits<Scalar,ghost_vec_t>::MvNorm(mv,normvec): MvNorm only accepts TwoNorm up to now.");
      
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

    static void SetBlock( const ghost_vec_t& A, const std::vector<int>& index, ghost_vec_t& mv )
    {
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION((typename std::vector<int>::size_type)A.traits.nvecs < index.size(),std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,ghost_vec_t>::SetBlock(A,index,mv): index must be the same size as A.");
#endif
      Teuchos::RCP<ghost_vec_t > mvsub = CloneViewNonConst(mv,index);
      ghost_mvec_t* Asub = std::const_cast<ghost_mvec_t*>(&A);
      if ((typename std::vector<int>::size_type)A.traits.nvecs > index.size()) {
        Asub = CloneViewNonConst(*Asub,Teuchos::Range1D(0,index.size()-1)).get();
      }
    //TODO - copying vector data in ghost???
    //(*mvsub) = (*Asub);
    }

    static void
    SetBlock (const ghost_vec_t& A, 
	      const Teuchos::Range1D& index, 
	      ghost_vec_t& mv)
    {
    
      // We've already validated the static casts above.
      const int numColsA = A.traits.nvecs;
      const int numColsMv = mv.traits.nvecs;
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

//      TODO
//      *mv_view = *A_view; 
    }

    static void
    Assign (const ghost_vec_t& A, 
	    ghost_vec_t& mv)
    {

      // Range1D bounds are signed; size_t is unsigned.
      // Assignment of ghost_vec_t is a deep copy.

      const int numColsA = A.traits.nvecs;
      const int numColsMv = mv.traits.nvecs;
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
      //TODO
      if (numColsA == numColsMv)
        //TODO
        mv = A;
      else
        {
          Teuchos::RCP<ghost_vec_t > mv_view = 
            CloneViewNonConst (mv, Teuchos::Range1D(0, numColsA-1));
          //TODO
          *mv_view = A;
        }
    }


//TODO
    static void MvRandom( ghost_vec_t& mv )
    { 
      mv.randomize(); 
    }
//TODO
    static void MvInit( ghost_vec_t& mv, Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero() )
    { mv.putScalar(alpha); }
//TODO
    static void MvPrint( const ghost_vec_t& mv, std::ostream& os )
    { 
      Teuchos::FancyOStream fos(Teuchos::rcp(&os,false));
      mv.describe(fos,Teuchos::VERB_EXTREME);
    }

  ghost_vec_t* createGhostViewOfTeuchosSDM
        (const Teuchos::SerialDenseMatrix<lidx_t,Scalar>& M)
  {    
      ghost_vtraits_t dmtraits = GHOST_VTRAITS_INIT(
                .flags = GHOST_VEC_DEFAULTS, 
                .nrows=M.numRows(),
                .nvecs=M.numCols(), 
                .nrowspadded=M.stride(),
                .datatype=phist::ScalarTraits<Scalar>::ghost_dt);
                                                                                                                                                  
      ghost_vec_t* Mghost=ghost_createVector(&dmtraits);
      Mghost->viewPlain(Mghost, (void*)M.Values(),
                Mghost->traits->nrows, Mghost->traits->nvec,
                0,0,ghost->traits->nrowspadded);
  return Mghost;
  } 


#ifdef HAVE_BELOS_TSQR
    /// \typedef tsqr_adaptor_type
    /// \brief TsqrAdaptor specialization for ghost_vec_t
    ///
    typedef GhostTsqrAdaptor<Scalar> tsqr_adaptor_type;
#endif // HAVE_BELOS_TSQR
  }; 

} // end of Belos namespace 

#endif 
// end of file BELOS_GHOST_ADAPTER_HPP
