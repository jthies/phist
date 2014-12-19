#define PHIST_rcp phist::PREFIX(rcp)

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

#ifdef HAVE_BELOS_TSQR
    /// \typedef tsqr_adaptor_type
    /// \brief TsqrAdaptor specialization for ghost_densemat_t
    ///
    typedef phist::TsqrAdaptor<Scalar> tsqr_adaptor_type;
#endif // HAVE_BELOS_TSQR

    static Teuchos::RCP<MV> Clone( const MV& mv, const int numvecs )
    {
      PHIST_ENTER_FCN(__FUNCTION__);
      TYPE(const_mvec_ptr) mv_in;
      TYPE(mvec_ptr) mv_out=NULL;
      mv_in=mv.get();
      const_map_ptr_t map=NULL;
      PHIST_TCHK_IERR(SUBR(mvec_get_map)(mv_in,&map,&iflag_),iflag_);
      PHIST_TCHK_IERR(SUBR(mvec_create)(&mv_out,map,(int)numvecs,&iflag_),iflag_);
      return PHIST_rcp(mv_out,true);
    }

    static Teuchos::RCP<MV> CloneCopy( const MV& mv )
    {
      PHIST_ENTER_FCN(__FUNCTION__);
      int nvec=GetNumberVecs(mv);
      Teuchos::RCP<MV> v_out=Clone(mv,nvec);
      PHIST_TCHK_IERR(SUBR(mvec_get_block)(mv.get(),v_out->get(),0,nvec-1,&iflag_),iflag_);
      return v_out;
    }

    static Teuchos::RCP<MV> CloneCopy( const MV& mv, const std::vector<int>& index )
    {
      PHIST_ENTER_FCN(__FUNCTION__);

      if (is_contig(index))
      {
        Teuchos::Range1D r(index[0],*(index.end()-1));
        return CloneCopy(mv,r);
      }

      // not contiguous: clone column by column
      // (we might pick contiguous blocks, but 
      // for the moment it is OK as it like this)
      TYPE(const_mvec_ptr) v_in;
      TYPE(mvec_ptr) v_out, v_col=NULL;
      const_map_ptr_t map;
      int nvec_in, nvec_out;
      v_in=mv.get();
      PHIST_TCHK_IERR(SUBR(mvec_num_vectors)(v_in,&nvec_in,&iflag_),iflag_);
      PHIST_TCHK_IERR(SUBR(mvec_get_map)(v_in,&map,&iflag_),iflag_);

      nvec_out=index.size();
      PHIST_TCHK_IERR(SUBR(mvec_create)(&v_out,map,nvec_out,&iflag_),iflag_);
      
      for (int i=0;i<nvec_out;i++)
      {
        PHIST_TCHK_IERR(SUBR(mvec_view_block)(v_out,&v_col,i,i,&iflag_),iflag_);
        PHIST_TCHK_IERR(SUBR(mvec_get_block)(v_in,v_col,index[i],index[i],&iflag_),iflag_);
      }
      if (v_col!=NULL)
      {
        PHIST_TCHK_IERR(SUBR(mvec_delete)(v_col,&iflag_),iflag_);
      }
      
      return PHIST_rcp(v_out,true);
    }

    static Teuchos::RCP<MV>
    CloneCopy (const MV& mv, 
	       const Teuchos::Range1D& index)
    {
      PHIST_ENTER_FCN(__FUNCTION__);
      TYPE(const_mvec_ptr) v_in;
      TYPE(mvec_ptr) v_out;
      const_map_ptr_t map;
      int nvec_out=index.ubound()-index.lbound()+1;
      v_in=mv.get();
      PHIST_TCHK_IERR(SUBR(mvec_get_map)(v_in,&map,&iflag_),iflag_);

      PHIST_TCHK_IERR(SUBR(mvec_create)(&v_out,map,nvec_out,&iflag_),iflag_);
      
      PHIST_TCHK_IERR(SUBR(mvec_get_block)(v_in,v_out,index.lbound(),index.ubound(),&iflag_),iflag_);
      return PHIST_rcp(v_out,true);
    }


    static Teuchos::RCP<MV> CloneViewNonConst( MV& mv, 
    const std::vector<int>& index )
    {
      PHIST_ENTER_FCN(__FUNCTION__);
      if (is_contig(index))
      {
        Teuchos::Range1D r(index[0],*(index.end()-1));
        return CloneViewNonConst(mv,r);
      }
      // the phist interface does not supporte `scattered views'.
      // If it turns out that we need them here we could build the
      // functionality into the phist::MultiVector wrapper by storing
      // an array of contiguous views.
      PHIST_TCHK_IERR(iflag_=-99,iflag_);
      return Teuchos::null;
    }

    static Teuchos::RCP<MV> 
    CloneViewNonConst (MV& mv, 
		       const Teuchos::Range1D& index)
    {
      PHIST_ENTER_FCN(__FUNCTION__);
      TYPE(mvec_ptr) v_in, v_out=NULL;
      v_in=mv.get();

      PHIST_TCHK_IERR(SUBR(mvec_view_block)(v_in,&v_out,index.lbound(),index.ubound(),&iflag_),iflag_);
      return PHIST_rcp(v_out,true);
    }


    static Teuchos::RCP<const MV > CloneView(const MV& mv, const std::vector<int>& index )
    {
      PHIST_ENTER_FCN(__FUNCTION__);    
      return Teuchos::rcp_dynamic_cast<const MV >
        (CloneViewNonConst(const_cast<MV&>(mv),index));
    }

    static Teuchos::RCP<const MV > 
    CloneView (const MV& mv, 
	       const Teuchos::Range1D& index)
    {
      PHIST_ENTER_FCN(__FUNCTION__);    
      return Teuchos::rcp_dynamic_cast<const MV >
        (CloneViewNonConst(const_cast<MV&>(mv),index));
    }

    static int GetVecLength( const MV& mv )
    {
      PHIST_ENTER_FCN(__FUNCTION__);
      // I think this should return the global length, but
      // I don't see why that would ever be necessary (?)
      // Anyway, the phist interface gives access only to 
      // the local length.
      int local_length;
      PHIST_TCHK_IERR(SUBR(mvec_my_length)(mv.get(),&local_length,&iflag_),iflag_);
      return local_length;
    }

  static int GetNumberVecs( const MV& mv )
  { 
      PHIST_ENTER_FCN(__FUNCTION__);
      int nvec;
      PHIST_TCHK_IERR(SUBR(mvec_num_vectors)(mv.get(),&nvec,&iflag_),iflag_);
      return nvec;
  }


    static bool HasConstantStride( const MV& mv )
    { 
      PHIST_ENTER_FCN(__FUNCTION__);
      // the phist interface does not currently allow
      // 'scattered' vectors as they are called in ghost.
      return true;
    }

    static void MvTimesMatAddMv( Scalar alpha, const MV& A, 
                                 const Teuchos_sdMat_t& B, 
                                 Scalar beta, MV& mv )
    {
      PHIST_ENTER_FCN(__FUNCTION__);
      PHIST_TCHK_IERR(SUBR(mvec_times_sdMat)(alpha,A.get(),convertSDM(B),beta,mv.get(),&iflag_),iflag_);
    }

    // compute mv = alpha*A + beta*B. This function is abused in Belos by aliasing mv and A 
    // or B, e.g. A = 0*A+1*B instead of A=B. We therefore have to be a bit careful with
    // the memcopy we use if either alpha or beta are 0.
    static void MvAddMv( Scalar alpha, const MV& A, Scalar beta, const MV& B, MV& mv )
    {
      PHIST_ENTER_FCN(__FUNCTION__);
      PHIST_TCHK_IERR(SUBR(mvec_add_mvec)(alpha,A.get(),st::zero(),mv.get(),&iflag_),iflag_);
      PHIST_TCHK_IERR(SUBR(mvec_add_mvec)(beta,B.get(),st::one(),mv.get(),&iflag_),iflag_);
    }

    static void MvScale ( MV& mv, Scalar alpha )
    {
      PHIST_ENTER_FCN(__FUNCTION__);
      PHIST_TCHK_IERR(SUBR(mvec_scale)(mv.get(),alpha,&iflag_),iflag_);
    }

    static void MvScale ( MV& mv, const std::vector<Scalar>& alphas )
    {
      PHIST_ENTER_FCN(__FUNCTION__);
      PHIST_TCHK_IERR(SUBR(mvec_vscale)(mv.get(),&alphas[0],&iflag_),iflag_);
    }

    // C=alpha*A*B
    static void MvTransMv( Scalar alpha, const MV& A, const MV& B, Teuchos_sdMat_t& C)
    {
      PHIST_ENTER_FCN(__FUNCTION__);
      PHIST_TCHK_IERR(SUBR(mvecT_times_mvec)(alpha,A.get(),B.get(),st::zero(),convertSDM(C),&iflag_),iflag_);
    }

    static void MvDot( const MV& A, const MV& B, std::vector<Scalar> &dots)
    {
      PHIST_ENTER_FCN(__FUNCTION__);    
      PHIST_TCHK_IERR(SUBR(mvec_dot_mvec)(A.get(),B.get(),&dots[0],&iflag_),iflag_);
    }

    static void MvNorm(const MV& mv, std::vector<magn_t> &normvec, NormType type=TwoNorm)
    {
      PHIST_ENTER_FCN(__FUNCTION__);
      if (type!=TwoNorm)
      {
        PHIST_SOUT(PHIST_WARNING,"can only compute 2-norm in this interface\n"
                                 "(file %s, line %d)\n",__FILE__,__LINE__);
      }
      PHIST_TCHK_IERR(SUBR(mvec_norm2)(mv.get(),&normvec[0],&iflag_),iflag_);
    }

    static void SetBlock( const MV& A, const std::vector<int>& index, MV& mv )
    {
      PHIST_ENTER_FCN(__FUNCTION__);
      if (is_contig(index))
      {
      SetBlock(A,Teuchos::Range1D(index[0],*(index.end()-1)), mv);
      }
      else
      {
        PHIST_SOUT(PHIST_VERBOSE,"wARNING, extra copy operation due to missing implementation\n"
                              "         (file %s, line %d)\n",__FILE__,__LINE__);
      Teuchos::RCP<MV> mv_tmp=CloneCopy(A,index);
      Assign(*mv_tmp,mv);
      }
    }

    static void
    SetBlock (const MV& A, 
	      const Teuchos::Range1D& index, 
	      MV& mv)
    {
      PHIST_ENTER_FCN(__FUNCTION__);
      PHIST_CHK_IERR(SUBR(mvec_set_block)(mv.get(),A.get(),index.lbound(),index.ubound(),&iflag_),iflag_);
    }

    static void
    Assign (const MV& A, 
	    MV& mv)
    {
      PHIST_ENTER_FCN(__FUNCTION__);
      int nvecA;
      PHIST_TCHK_IERR(SUBR(mvec_num_vectors)(A.get(),&nvecA,&iflag_),iflag_);
      PHIST_TCHK_IERR(SUBR(mvec_get_block)(A.get(),mv.get(),0,nvecA-1,&iflag_),iflag_);
    }


    static void MvRandom( MV& mv )
    { 
      PHIST_ENTER_FCN(__FUNCTION__); 
      PHIST_TCHK_IERR(SUBR(mvec_random)(mv.get(),&iflag_),iflag_);
    }

    static void MvInit( MV& mv, Scalar alpha = st::zero() )
    {
      PHIST_ENTER_FCN(__FUNCTION__); 
      PHIST_TCHK_IERR(SUBR(mvec_put_value)(mv.get(),alpha,&iflag_),iflag_);
    }

    static void MvPrint( const MV& mv, std::ostream& os )
    {
      PHIST_ENTER_FCN(__FUNCTION__);
      //TODO: we could probably mess with
      //      stream buffers etc to redirect
      //      the output from mvec_print to 
      //      the given stream. Alternatively
      //      we could extract the raw data and
      //      print it manually, but as we can
      //      simply switch to a kernel lib that
      //      has its own interface (like Tpetra),
      //      this function is not that important.
      os << "MVPrint not implemented in file "<<__FILE__<<std::endl;
      return;
    }

  // private helper function
  static TYPE(const_sdMat_ptr) convertSDM
        (const Teuchos_sdMat_t& M)
  {
      PHIST_ENTER_FCN(__FUNCTION__);
      return (TYPE(const_sdMat_ptr))convertSDM((Teuchos_sdMat_t&)(M));
  }

  // private helper function
  static TYPE(sdMat_ptr) convertSDM
        (Teuchos_sdMat_t& M)
  {
      PHIST_ENTER_FCN(__FUNCTION__);
      TYPE(sdMat_ptr) M_out=0;
      // oh oh, memory leak
      comm_ptr_t comm=NULL;
      PHIST_TCHK_IERR(phist_comm_create(&comm,&iflag_),iflag_);
      PHIST_TCHK_IERR(SUBR(sdMat_create_view)(&M_out,comm,M.values(),
        M.stride(),M.numRows(),M.numCols(),&iflag_),iflag_);
    return M_out;
  }

private:

  static int iflag_;


private:

static bool is_contig(const std::vector<int>& index)
{
      bool contig=true;
      for (typename std::vector<int>::size_type j=1; j<index.size(); ++j) 
      {
        if (index[j] != index[j-1]+1) {
          contig=false;
          break;
        }
      }
      return contig;
}

};

} // end of Belos namespace 

#undef PHIST_rcp