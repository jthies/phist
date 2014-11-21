
namespace Belos {

using ::phist::ScalarTraits;

  template<>
  class MultiVecTraits<_ST_, ::phist::MultiVector< _ST_ > >
  {  
  public:

    typedef _ST_ Scalar;
    typedef ::phist::MultiVector< _ST_ > MV;
    typedef ::phist::ScalarTraits<Scalar> st;
    typedef typename st::magn_t magn_t;
    //! serial dense matrix from Teuchos, we need this for e.g. the BLAS interface.
    //! Note: the index type *must* be int here, not int64_t, so we decided to have
    //! phist local indices ints, even if ghost uses int64_t.
    typedef Teuchos::SerialDenseMatrix<int,Scalar> Teuchos_sdMat_t;

    static Teuchos::RCP<MV > Clone( const MV& mv, const int numvecs )
    {
      ENTER_FCN(__FUNCTION__);
      PHIST_TCHK_IERR(ierr_=-99,ierr_);
    }

    static Teuchos::RCP<MV > CloneCopy( const MV& mv )
    {
      ENTER_FCN(__FUNCTION__);
      PHIST_TCHK_IERR(ierr_=-99,ierr_);
    }

    static Teuchos::RCP<MV > CloneCopy( const MV& mv, const std::vector<int>& index )
    {
      ENTER_FCN(__FUNCTION__);
      PHIST_TCHK_IERR(ierr_=-99,ierr_);      
    }

    static Teuchos::RCP<MV > 
    CloneCopy (const MV& mv, 
	       const Teuchos::Range1D& index)
    {
      ENTER_FCN(__FUNCTION__);
      PHIST_TCHK_IERR(ierr_=-99,ierr_);      
    }


    static Teuchos::RCP<MV > CloneViewNonConst( MV& mv, 
    const std::vector<int>& index )
    {
      ENTER_FCN(__FUNCTION__);
      PHIST_TCHK_IERR(ierr_=-99,ierr_);      
    }

    static Teuchos::RCP<MV > 
    CloneViewNonConst (MV& mv, 
		       const Teuchos::Range1D& index)
    {
      ENTER_FCN(__FUNCTION__);
      PHIST_TCHK_IERR(ierr_=-99,ierr_);      
    }


    static Teuchos::RCP<const MV > CloneView(const MV& mv, const std::vector<int>& index )
    {
    }

    static Teuchos::RCP<const MV > 
    CloneView (const MV& mv, 
	       const Teuchos::Range1D& index)
    {
      ENTER_FCN(__FUNCTION__);    
      return Teuchos::rcp_dynamic_cast<const MV >
        (CloneViewNonConst(const_cast<MV&>(mv),index));
    }

    static int GetVecLength( const MV& mv )
    {
      ENTER_FCN(__FUNCTION__);
      PHIST_TCHK_IERR(ierr_=-99,ierr_);      
    }

  static int GetNumberVecs( const MV& mv )
  { 
      ENTER_FCN(__FUNCTION__);
      PHIST_TCHK_IERR(ierr_=-99,ierr_);      
  }


    static bool HasConstantStride( const MV& mv )
    { 
      ENTER_FCN(__FUNCTION__);
      PHIST_TCHK_IERR(ierr_=-99,ierr_);      
    }

    static void MvTimesMatAddMv( Scalar alpha, const MV& A, 
                                 const Teuchos_sdMat_t& B, 
                                 Scalar beta, MV& mv )
    {
      ENTER_FCN(__FUNCTION__);
      PHIST_TCHK_IERR(ierr_=-99,ierr_);      
    }

    // compute mv = alpha*A + beta*B. This function is abused in Belos by aliasing mv and A 
    // or B, e.g. A = 0*A+1*B instead of A=B. We therefore have to be a bit careful with
    // the memcopy we use if either alpha or beta are 0.
    static void MvAddMv( Scalar alpha, const MV& A, Scalar beta, const MV& B, MV& mv )
    {
      ENTER_FCN(__FUNCTION__);
      PHIST_TCHK_IERR(ierr_=-99,ierr_);      
    }

    static void MvScale ( MV& mv, Scalar alpha )
    {
      ENTER_FCN(__FUNCTION__);
      PHIST_TCHK_IERR(ierr_=-99,ierr_);      
    }

    static void MvScale ( MV& mv, const std::vector<Scalar>& alphas )
    {
      ENTER_FCN(__FUNCTION__);
      PHIST_TCHK_IERR(ierr_=-99,ierr_);      
    }

    // C=alpha*A*B
    static void MvTransMv( Scalar alpha, const MV& A, const MV& B, Teuchos_sdMat_t& C)
    {
      ENTER_FCN(__FUNCTION__);
      PHIST_TCHK_IERR(ierr_=-99,ierr_);      
    }

    static void MvDot( const MV& A, const MV& B, std::vector<Scalar> &dots)
    {
      ENTER_FCN(__FUNCTION__);    
      PHIST_TCHK_IERR(ierr_=-99,ierr_);      
    }

    static void MvNorm(const MV& mv, std::vector<magn_t> &normvec, NormType type=TwoNorm)
    {
      ENTER_FCN(__FUNCTION__);    
      PHIST_TCHK_IERR(ierr_=-99,ierr_);      
    }

    static void SetBlock( const MV& A, const std::vector<int>& index, MV& mv )
    {
      ENTER_FCN(__FUNCTION__);    
      PHIST_TCHK_IERR(ierr_=-99,ierr_);      
    }

    static void
    SetBlock (const MV& A, 
	      const Teuchos::Range1D& index, 
	      MV& mv)
    {
      ENTER_FCN(__FUNCTION__);        
      PHIST_TCHK_IERR(ierr_=-99,ierr_);      
    }

    static void
    Assign (const MV& A, 
	    MV& mv)
    {
      ENTER_FCN(__FUNCTION__); 
      PHIST_TCHK_IERR(ierr_=-99,ierr_);
    }


    static void MvRandom( MV& mv )
    { 
      ENTER_FCN(__FUNCTION__); 
      PHIST_TCHK_IERR(ierr_=-99,ierr_);
    }

    static void MvInit( MV& mv, Scalar alpha = st::zero() )
    {
      ENTER_FCN(__FUNCTION__); 
      PHIST_TCHK_IERR(ierr_=-99,ierr_);
    }

    static void MvPrint( const MV& mv, std::ostream& os )
    {
      ENTER_FCN(__FUNCTION__); 
      PHIST_TCHK_IERR(ierr_=-99,ierr_);
    }

  // private helper function
  static TYPE(mvec_ptr) createPhistViewOfTeuchosSDM
        (const Teuchos_sdMat_t& M)
  {
      ENTER_FCN(__FUNCTION__); 
      PHIST_TCHK_IERR(ierr_=-99,ierr_);
  }

private:

  static int ierr_;

};

} // end of Belos namespace 

