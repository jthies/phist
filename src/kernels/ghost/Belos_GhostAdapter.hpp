#ifndef BELOS_GHOST_ADAPTER_HPP
#define BELOS_GHOST_ADAPTER_HPP

#include "ghost.h"
#include "phist_typedefs.h"
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

#ifdef HAVE_BELOS_TSQR
#  include <Ghost_TsqrAdaptor.hpp>
#endif // HAVE_BELOS_TSQR


// this file is mostly copied from the Belos Tpetra adapter implementation in Trilinos 11.2.4

namespace Belos {

using ::phist::GhostMV;

  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::MultiVecTraits for ghost_densemat_t.
  //
  ////////////////////////////////////////////////////////////////////

  /*!  \brief Template specialization of Belos::MultiVecTraits class using the ghost_densemat_t class.

    This interface will ensure that any ghost_densemat_t will be accepted by the Belos
    templated solvers.  */
  template<class Scalar>
  class MultiVecTraits<Scalar, GhostMV >
  {  
  public:

    typedef ::phist::ScalarTraits<Scalar> st;
    typedef typename st::magn_t magn_t;
    typedef typename Traits<Scalar>::Teuchos_sdMat_t Teuchos_sdMat_t;

    static Teuchos::RCP<GhostMV > Clone( const GhostMV& mv, const int numvecs )
    {
    ENTER_FCN(__FUNCTION__);    
      ghost_densemat_t* _mv = const_cast<GhostMV&>(mv).get();
      ghost_densemat_traits_t vtraits = _mv->traits;
      // copy the data even if the input vector is itself a view
      // (bitwise NAND operation to unset the view flag if set)
      vtraits.flags = (ghost_densemat_flags_t)((int)vtraits.flags & ~(int)GHOST_DENSEMAT_VIEW);
      vtraits.ncols=numvecs;
      ghost_densemat_t* mv_clone;
      ghost_densemat_create(&mv_clone,_mv->context,vtraits);
      // this allocates the memory for the vector
      Scalar z=st::zero();
      mv_clone->fromScalar(mv_clone,(void*)&z);
      return phist::rcp(mv_clone,true);
    }

    static Teuchos::RCP<GhostMV > CloneCopy( const GhostMV& mv )
    {
    ENTER_FCN(__FUNCTION__);    
      ghost_densemat_t* _mv = const_cast<GhostMV&>(mv).get();
      ghost_densemat_traits_t vtraits = _mv->traits;
      // copy the data even if the input vector is itself a view
      // (bitwise NAND operation to unset the view flag if set)
      vtraits.flags = (ghost_densemat_flags_t)((int)vtraits.flags & ~(int)GHOST_DENSEMAT_VIEW);
      ghost_densemat_t* mv_clone;
      ghost_densemat_create(&mv_clone,_mv->context,vtraits);
      mv_clone->fromVec(mv_clone,_mv,0);
      return phist::rcp(mv_clone,true); 
    }

    static Teuchos::RCP<GhostMV > CloneCopy( const GhostMV& mv, const std::vector<int>& index )
    {
      ENTER_FCN(__FUNCTION__);    
      ghost_densemat_t* _mv = const_cast<GhostMV&>(mv).get();
      TEUCHOS_TEST_FOR_EXCEPTION(index.size() == 0,std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,GhostMV>::CloneCopy(mv,index): numvecs must be greater than zero.");

      TEUCHOS_TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::runtime_error,
          "Belos::MultiVecTraits<Scalar,GhostMV>::CloneCopy(mv,index): indices must be >= zero.");

      TEUCHOS_TEST_FOR_EXCEPTION( (size_t)*std::max_element(index.begin(),index.end()) >= GetNumberVecs(mv), std::runtime_error,
          "Belos::MultiVecTraits<Scalar,GhostMV>::CloneCopy(mv,index): indices must be < mv.traits.ncols.");

      bool contig=true;
      for (typename std::vector<int>::size_type j=1; j<index.size(); ++j) {
        if (index[j] != index[j-1]+1) {
          // not contiguous; short circuit
          contig=false;
          break;
        }
      }

      if (contig)
      {
        ghost_densemat_t *result = NULL;
        _mv->clone(_mv,&result,index.size(),index[0]);

        return phist::rcp(result,true);
      }
      else
      {
        ghost_densemat_t* result;
        ghost_densemat_traits_t vtraits = _mv->traits;
                vtraits.ncols=index.size();
        // copy the data even if the input vector is itself a view
        // (bitwise NAND operation to unset the view flag if set)
        vtraits.flags = (ghost_densemat_flags_t)((int)vtraits.flags & ~(int)GHOST_DENSEMAT_VIEW);
        ghost_densemat_create(&result,_mv->context,vtraits);
        // allocates memory
        Scalar zero=Teuchos::ScalarTraits<Scalar>::zero();
        result->fromScalar(result, &zero);
        // copy columns one by one
        for (int j=0;j<index.size();j++)
        {
          ghost_densemat_t *result_j;
          result->viewVec(result, &result_j, 1, j);
          result_j->fromVec(result_j,_mv,index[j]);
        }
        return phist::rcp(result,true);
      }
    }

    static Teuchos::RCP<GhostMV > 
    CloneCopy (const GhostMV& mv, 
	       const Teuchos::Range1D& index)
    {
      ENTER_FCN(__FUNCTION__);    
      ghost_densemat_t* _mv = const_cast<GhostMV&>(mv).get();
      const bool validRange = index.size() > 0 && 
	index.lbound() >= 0 && 
	index.ubound() < GetNumberVecs(mv);
      if (! validRange)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<Scalar, GhostMV<...> >::"
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
      return phist::rcp(_mv->clone(_mv,index.ubound()-index.lbound()+1,index.lbound),true);
    }


    static Teuchos::RCP<GhostMV > CloneViewNonConst( GhostMV& mv, const std::vector<int>& index )
    {
      ENTER_FCN(__FUNCTION__);    
      ghost_densemat_t* _mv = mv.get();
      TEUCHOS_TEST_FOR_EXCEPTION(index.size() == 0,std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,GhostMV>::CloneView(mv,index): numvecs must be greater than zero.");
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION( *std::min_element(index.begin(),index.end()) < 0, std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,GhostMV>::CloneView(mv,index): indices must be >= zero.");
      TEUCHOS_TEST_FOR_EXCEPTION( (size_t)*std::max_element(index.begin(),index.end()) >= mv.traits.ncols, std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,GhostMV>::CloneView(mv,index): indices must be < mv.traits.ncols.");
#endif

      bool constStride=true;                                                                                                      
      int stride=1;
      if (_mv->traits.flags&GHOST_DENSEMAT_SCATTERED) 
      {
        constStride=false;
      }
      else
      {
        if (index.size()>1) stride=index[1]-index[0];                                                                                                         
        for (typename std::vector<int>::size_type j=2; j<index.size(); ++j)                                                         
        {                                                                                                                         
          if (index[j] != index[j-1]+stride)                                                                                      
          {                                                                                                                     
            constStride=false;                                                                                                    
            break;                                                                                                                
          }
        }
      }
        
    ghost_densemat_t* result=NULL;

    if (constStride==false)
    {
#ifdef GHOST_HAVE_LONGIDX
      // ghost expects long ints here, while we get ints. So we copy them over:
      std::vector<ghost_idx_t> clone_index(index.size());
      for (int i=0;i<index.size();i++)
      {
        clone_index[i]=(ghost_idx_t)index[i];
      }
#else
      const std::vector<ghost_idx_t>& clone_index=index;
#endif
      _mv->viewScatteredVec(_mv,&result,(ghost_idx_t)index.size(),(ghost_idx_t*)&clone_index[0]);
    }
    else
    {
      // constant stride
      
      // stride k: first simply view the vector, then manually set pointers and stride
      _mv->viewVec(_mv,&result,index.size(),index[0]);
      if (stride!=1)
      {
        for (int i=0;i<index.size();i++)
        {
          result->val[i] = _mv->val[index[i]];
        }
        //TODO: this is quite a nasty hack to allow the ghost_gemm function to
        // work for constant stride access, it should be fixed somehow in ghost
        // so that the GPU stuff works as well etc.
        result->traits.nrowspadded = _mv->traits.nrowspadded*stride;
      }
    }
    return phist::rcp(result,true);
  }


    static Teuchos::RCP<GhostMV > 
    CloneViewNonConst (GhostMV& mv, 
		       const Teuchos::Range1D& index)
    {
      ENTER_FCN(__FUNCTION__);    
      ghost_densemat_t* _mv=mv.get();
      // NOTE (mfh 11 Jan 2011) We really should check for possible
      // overflow of int here.  However, the number of columns in a
      // multivector typically fits in an int.
      const int numCols = static_cast<int> (_mv->traits.ncols);
      const bool validRange = index.size() > 0 && 
	index.lbound() >= 0 && index.ubound() < numCols;
      if (! validRange)
	{
	  std::ostringstream os;
	  os << "Belos::MultiVecTraits<Scalar, GhostMV<...> >::"
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
      
      ghost_densemat_t* result=NULL;
      _mv->viewVec(_mv,&result,index.ubound()-index.lbound()+1, index.lbound());

      return phist::rcp(result,true);
    }


    static Teuchos::RCP<const GhostMV > CloneView(const GhostMV& mv, const std::vector<int>& index )
    {
      ENTER_FCN(__FUNCTION__);    
      return Teuchos::rcp_dynamic_cast<const GhostMV >
        (CloneViewNonConst(const_cast<GhostMV&>(mv),index));
    }

    static Teuchos::RCP<const GhostMV > 
    CloneView (const GhostMV& mv, 
	       const Teuchos::Range1D& index)
    {
      ENTER_FCN(__FUNCTION__);    
      return Teuchos::rcp_dynamic_cast<const GhostMV >
        (CloneViewNonConst(const_cast<GhostMV&>(mv),index));
    }

    static int GetVecLength( const GhostMV& mv )
    {
// TOOD - this is not available in ghost,
// but as it is only used for confirming 
// that vectors are compatible for certain
// operations in Belos, I think we can
// cheat here and simply return the local vec length instead.
    return mv.get()->traits.nrows;
  }

  static int GetNumberVecs( const GhostMV& mv )
  { 
    return mv.get()->traits.ncols;
  }


    static bool HasConstantStride( const GhostMV& mv )
    { return (mv.get()->traits.flags&GHOST_DENSEMAT_SCATTERED==false); }

    static void MvTimesMatAddMv( Scalar alpha, const GhostMV& A, 
                                 const Teuchos_sdMat_t& B, 
                                 Scalar beta, GhostMV& mv )
    {
      ENTER_FCN(__FUNCTION__);    
      // create view of Teuchos matrix as GhostMV
      ghost_densemat_t* Bghost=createGhostViewOfTeuchosSDM(B);
      ghost_densemat_t* _A = (ghost_densemat_t*)A.get();
      ghost_densemat_t* _mv = (ghost_densemat_t*)mv.get();
      // multiply
      const char* trans="N";
      ghost_gemm(_mv,_A,Bghost,(char*)trans,&alpha,&beta,GHOST_GEMM_NO_REDUCE);
      
      Bghost->destroy(Bghost);
    }

    // compute mv = alpha*A + beta*B. This function is abused in Belos by aliasing mv and A 
    // or B, e.g. A = 0*A+1*B instead of A=B. We therefore have to be a bit careful with
    // the memcopy we use if either alpha or beta are 0.
    static void MvAddMv( Scalar alpha, const GhostMV& A, Scalar beta, const GhostMV& B, GhostMV& mv )
    {
      ENTER_FCN(__FUNCTION__);
      ghost_densemat_t* _mv = mv.get();
      ghost_densemat_t* _A = (ghost_densemat_t*)A.get();
      ghost_densemat_t* _B = (ghost_densemat_t*)B.get();
      
      Scalar zero=st::zero();
      Scalar one=st::one();

      // mv = A
      if (alpha==zero)
      {
        if (beta==zero)
        {
          _mv->fromScalar(_mv,(void*)&zero);
        }
        else
        {
          // Belos allows aliasing here, putting views of A as MV etc. (cf. MVOPTester).
          // This is a problem because fromVec uses memcpy, which does not allow aliasing.
          bool mv_is_B = (_B->val[0] == _mv->val[0]) &&
                         (_B->val[_B->traits.ncols-1]==_mv->val[_mv->traits.ncols-1]) &&
                         (_B->traits.nrowspadded == _mv->traits.nrowspadded);
          // NOTE: we do not check for partial overlap, or if one of the two is a 'scattered 
          // view' or something of the kind. So this function may cause problems if used in
          // unexpected ways because of the memcpy here. 
          if (mv_is_B==false)
          {
            _mv->fromVec(_mv,_B,0);
          }
          if (beta!=one)
          {
            _mv->scale(_mv,&beta);
          }
        }
      }
      else // alpha!=0
      {
      // cf. comment on aliasing above.
        bool mv_is_A = (_A->val[0] == _mv->val[0]) &&
                         (_A->val[_A->traits.ncols-1]==_mv->val[_mv->traits.ncols-1]) &&
                         (_A->traits.nrowspadded == _mv->traits.nrowspadded);
          if (mv_is_A==false)
          {
            _mv->fromVec(_mv,_A,0);
          }
        if (alpha!=one)
        {
          _mv->scale(_mv,&alpha);
        }
        // now we have mv=alpha*A
        
        // mv = alpha*A + beta*B
        if (beta!=zero)
        {
          _mv->axpy(_mv,_B,(void*)&beta);
        }
      }
      //std::cerr << "alpha*A+beta*B="<<std::endl;
      //_mv->print(_mv);
    }

    static void MvScale ( GhostMV& mv, Scalar alpha )
    {
      ENTER_FCN(__FUNCTION__);    
      ghost_densemat_t* _mv = const_cast<GhostMV&>(mv).get();
      _mv->scale(_mv,(void*)&alpha); 
    }

    static void MvScale ( GhostMV& mv, const std::vector<Scalar>& alphas )
    {
      ENTER_FCN(__FUNCTION__);    
      ghost_densemat_t* _mv = const_cast<GhostMV&>(mv).get();
      void* val = (void*)(&alphas[0]);
      _mv->vscale(_mv,val);
    }

    // C=alpha*A*B
    static void MvTransMv( Scalar alpha, const GhostMV& A, const GhostMV& B, Teuchos_sdMat_t& C)
    {
      ENTER_FCN(__FUNCTION__);    
      ghost_densemat_t* Cghost=createGhostViewOfTeuchosSDM(C);

      Scalar beta = st::zero();
      const char* trans=phist::ScalarTraits<Scalar>::is_complex()? "C": "T";
      ghost_gemm(Cghost,  const_cast<ghost_densemat_t*>(A.get()),
                   const_cast<ghost_densemat_t*>(B.get()),
                   (char*)trans,
                   (void*)&alpha, (void*)&beta,
                   GHOST_GEMM_ALL_REDUCE);
      Cghost->destroy(Cghost);
    }

    static void MvDot( const GhostMV& A, const GhostMV& B, std::vector<Scalar> &dots)
    {
      ENTER_FCN(__FUNCTION__);    
      TEUCHOS_TEST_FOR_EXCEPTION(GetNumberVecs(A) != GetNumberVecs(B),std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,GhostMV>::MvDot(A,B,dots): A and B must have the same number of vectors.");

      TEUCHOS_TEST_FOR_EXCEPTION(dots.size() < (typename std::vector<int>::size_type)GetNumberVecs(A),std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,GhostMV>::MvDot(A,B,dots): dots must have room for all dot products.");

      Teuchos::ArrayView<Scalar> av(dots);
      ghost_dot((void*)&dots[0],const_cast<ghost_densemat_t*>(A.get()),const_cast<ghost_densemat_t*>(B.get()));
    }

    static void MvNorm(const GhostMV& mv, std::vector<magn_t> &normvec, NormType type=TwoNorm)
    {
      ENTER_FCN(__FUNCTION__);    
      ghost_densemat_t* _mv = const_cast<GhostMV&>(mv).get();
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION(normvec.size() < (typename std::vector<int>::size_type)mv.traits.ncols,std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,GhostMV>::MvNorm(mv,normvec): normvec must have room for all norms.");
#endif
      Teuchos::Array<Scalar> av(normvec.size());
      Teuchos::ArrayView<typename st::magn_t> nv(normvec);
      TEUCHOS_TEST_FOR_EXCEPTION(type != TwoNorm,std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,GhostMV>::MvNorm(mv,normvec): MvNorm only accepts TwoNorm up to now.");
      
      switch (type) {
        case OneNorm:
          break;
        case TwoNorm:
          ghost_dot(av.getRawPtr(), _mv, _mv);
          for (int i=0;i<av.size();i++)
          {
            nv[i]=st::real(st::sqrt(av[i]));
          }
          break;
        case InfNorm:
          break;
      }
    }

    static void SetBlock( const GhostMV& A, const std::vector<int>& index, GhostMV& mv )
    {
      ENTER_FCN(__FUNCTION__);    
#ifdef HAVE_TPETRA_DEBUG
      TEUCHOS_TEST_FOR_EXCEPTION((typename std::vector<int>::size_t)GetNumberVecs(A) < index.size(),std::invalid_argument,
          "Belos::MultiVecTraits<Scalar,GhostMV>::SetBlock(A,index,mv): index must be the same size as A.");
#endif
      // note the dual meaning of get() here: RCP.get() gives raw pointer to GhostMV,
      // GhostMV.get() gives raw pointer to ghost_densemat_t
      
      // view the columns that we want to set in mv:
      Teuchos::RCP<GhostMV> mvsub = CloneViewNonConst(mv,index);
      ghost_densemat_t* _mvsub = mvsub->get();
      Teuchos::RCP<GhostMV> Asub = Teuchos::null;
      ghost_densemat_t* _Asub = const_cast<ghost_densemat_t*>(A.get());
      if ((typename std::vector<int>::size_type)_Asub->traits.ncols > index.size()) {
        // this get is the GhostMV function to get an ghost_densemat_t*
        Asub = CloneViewNonConst(const_cast<GhostMV&>(A),Teuchos::Range1D(0,index.size()-1));
        _Asub= Asub->get();
      }
    _mvsub->fromVec(_mvsub,_Asub,0);
    return;
/*
    std::cout <<"MvCopy: indices: "<<std::endl;
    for (int i=0;i<index.size();i++)
      std::cout << index[i]<<" ";
    std::cout<<std::endl;
    
    std::cout << "complete target matrix of MvCopy():"<<std::endl;
    mv.get()->print(mv.get());

    std::cout << "target columns of MvCopy():"<<std::endl;
    _mvsub->print(_mvsub);

    std::cout << "input block of MvCopy() was:"<<std::endl;
    _Asub->print(_Asub);
*/
    }

    static void
    SetBlock (const GhostMV& A, 
	      const Teuchos::Range1D& index, 
	      GhostMV& mv)
    {
      ENTER_FCN(__FUNCTION__);        
      // We've already validated the static casts above.
      const int numColsA = GetNumberVecs(A);
      const int numColsMv = GetNumberVecs(mv);
      // 'index' indexes into mv; it's the index set of the target.
      const bool validIndex = index.lbound() >= 0 && index.ubound() < numColsMv;
      // We can't take more columns out of A than A has.
      const bool validSource = index.size() <= numColsA;

      if (! validIndex || ! validSource)
	{
	  std::ostringstream os;
	  os <<	"Belos::MultiVecTraits<Scalar, GhostMV<Scalar, ..."
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
      typedef Teuchos::RCP<GhostMV > MV_ptr;
      typedef Teuchos::RCP<const GhostMV > const_MV_ptr;

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

      mv_view->get()->fromVec(mv_view->get(),A_view->get(),0);
    }

    static void
    Assign (const GhostMV& A, 
	    GhostMV& mv)
    {
      ENTER_FCN(__FUNCTION__);        
      // Range1D bounds are signed; size_t is unsigned.
      // Assignment of GhostMV is a deep copy.

      const int numColsA = GetNumberVecs(A);
      const int numColsMv = GetNumberVecs(mv);
      if (numColsA > numColsMv)
      {
          std::ostringstream os;
          os <<	"Belos::MultiVecTraits<Scalar, GhostMV<Scalar, ..."
            "> >::Assign(A, mv): ";
          TEUCHOS_TEST_FOR_EXCEPTION(numColsA > numColsMv, std::invalid_argument,
                             os.str() << "Input multivector 'A' has " 
                             << numColsA << " columns, but output multivector "
                             "'mv' has only " << numColsMv << " columns.");
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
      }
      if (numColsA == numColsMv)
      {
        mv.get()->fromVec(mv.get(),(ghost_densemat_t*)A.get(),0);
      }
      else
      {
          ghost_densemat_t* mv_view = 
            CloneViewNonConst (mv, Teuchos::Range1D(0, numColsA-1))->get();
          // copy mv(:,1:numColsA)=A
          mv_view->fromVec(mv_view,(ghost_densemat_t*)A.get(),0);
      }
    }


    static void MvRandom( GhostMV& mv )
    { 
      mv.get()->fromRand(mv.get());
    }

    static void MvInit( GhostMV& mv, Scalar alpha = st::zero() )
    {
      mv.get()->fromScalar(mv.get(),(void*)&alpha);
    }

    static void MvPrint( const GhostMV& mv, std::ostream& os )
    {
      // TODO - the stream argument is ignored, ghost always prints to stdout
      ghost_densemat_t* _mv = const_cast<ghost_densemat_t*>(mv.get());
      _mv->print(_mv);
    }

  // private helper function
  static ghost_densemat_t* createGhostViewOfTeuchosSDM
        (const Teuchos_sdMat_t& M)
  {
    ENTER_FCN(__FUNCTION__);
      ghost_densemat_traits_t dmtraits;/*=new ghost_vtraits_t;*/
                dmtraits.flags = GHOST_DENSEMAT_DEFAULT;
                dmtraits.nrows=M.numRows();
                dmtraits.nrowshalo=M.numRows();
                dmtraits.nrowspadded=M.stride();
                dmtraits.ncols=M.numCols();
                dmtraits.datatype=st::ghost_dt;

      // The context and communicator are supposed to be irrelevant in an sdMat,
      // but it is not clear wether this is handled correctly everywhere i ghost.
      // For the moment we can afford to just put in MPI_COMM_WORLD at this point.
      MPI_Comm comm = MPI_COMM_WORLD;
      ghost_context_t* ctx=NULL;
      ghost_error_t gerr=ghost_context_create(&ctx,M.numRows(), M.numRows(), 
          GHOST_CONTEXT_DEFAULT, NULL, GHOST_SPARSEMAT_SRC_NONE, comm, 1.0);
      if (gerr!=GHOST_SUCCESS) PHIST_OUT(PHIST_ERROR,"GHOST error (%s) in file %s, line %d",
        phist_ghost_error2str(gerr),__FILE__,__LINE__);
      //TODO - check return values everywhere
      ghost_densemat_t* Mghost;
      ghost_densemat_create(&Mghost,ctx,dmtraits);
      Mghost->viewPlain(Mghost, (void*)M.values(),
                Mghost->traits.nrows, Mghost->traits.ncols,
                0,0,Mghost->traits.nrowspadded);

    return Mghost;
  }


#ifdef HAVE_BELOS_TSQR
    /// \typedef tsqr_adaptor_type
    /// \brief TsqrAdaptor specialization for ghost_densemat_t
    ///
    typedef ghost::TsqrAdaptor<Scalar> tsqr_adaptor_type;
#endif // HAVE_BELOS_TSQR
};

} // end of Belos namespace 

#endif 
// end of file BELOS_GHOST_ADAPTER_HPP
