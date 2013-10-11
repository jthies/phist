#include <algorithm>

#ifndef CLASSNAME
#error "file not included correctly"
#endif

/*! Test fixure. */
class CLASSNAME: public KernelTestWithVectors<_ST_,_N_,_NV_> 
  {

public:


  /*! Set up routine.
   */
  virtual void SetUp()
    {
    KernelTestWithVectors<_ST_,_N_,_NV_>::SetUp();
    if (typeImplemented_)
      {
      for (int j=0;j<nvec_;j++)
        for (int i=0;i<nloc_*stride_;i+=stride_)
          {
          vec1_vp_[j*lda_+i]=random_number();
          vec2_vp_[j*lda_+i]=st::one();
          }
      }
    }

  /*! Clean up.
   */
  virtual void TearDown() 
    {
    KernelTestWithVectors<_ST_,_N_,_NV_>::TearDown();
    }

};

  TEST_F(CLASSNAME, my_length)
    {
    if (typeImplemented_)
      {
      int nloc;
      SUBR(mvec_my_length)(vec1_,&nloc,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_EQ(nloc_, nloc); 
      }
    }

  TEST_F(CLASSNAME, num_vectors) 
    {
    if (typeImplemented_)
      {
      int nvec;
      SUBR(mvec_num_vectors)(vec1_,&nvec,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_EQ(nvec_, nvec);
      }
    }


  TEST_F(CLASSNAME, put_value) 
    {
    if (typeImplemented_)
      {
      ST val = (_ST_)42.0 + (ST)3.0*st::I();
      SUBR(mvec_put_value)(vec1_,val,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec1_vp_,nloc_,nvec_,lda_,stride_,val));
      }
    }

  TEST_F(CLASSNAME, dot_mvec)
    {
    if (typeImplemented_)
      {
      for (int j=0;j<nvec_;j++)
        for (int i=0;i<nloc_*stride_;i+=stride_)
          {
        vec2_vp_[j*lda_+i]=mt::one()/st::conj(vec1_vp_[j*lda_+i]);
        }
      _ST_* dots = new ST[nvec_];
      SUBR(mvec_dot_mvec)(vec1_,vec2_,dots,&ierr_);
      ASSERT_EQ(0,ierr_);
            
      _ST_ val = st::one() * (ST)nglob_;
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(dots,nvec_,1,nvec_,1,val));
      delete [] dots;
      }
    
    }

  TEST_F(CLASSNAME, random)
    {
    if (typeImplemented_)
      {
      SUBR(mvec_random)(vec1_,&ierr_);
      ASSERT_EQ(0,ierr_);
      MT absval[nloc_*nvec_];
      int k=0;
      for (int j=0;j<nvec_;j++)
        {
        for (int i=0;i<nloc_*stride_;i+=stride_)
          {
          absval[k++]=st::abs(vec1_vp_[j*lda_+i]);
          }
        }
      MT minval=1.0;
      std::sort(absval,absval+k);            
      for (int j=1;j<k;j++)
        {
        minval=std::min(minval,mt::abs(absval[j]-absval[j-1]));
        }
      // force assertion failure if two 'random' numbers are the same
      if (minval<mt::eps()) ASSERT_EQ(0,-1);
      }
    }

  // X = 1*Y + 0*X = Y
  TEST_F(CLASSNAME, copy_by_axpy)
    {
    if (typeImplemented_)
      {
      ST alpha = st::one();
      ST beta  = st::zero();
      SUBR(mvec_add_mvec)(alpha,vec1_,beta,vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);
            
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nloc_,nvec_,lda_,stride_));
      }
    
    }

  // X = 0*Y + a*X = a*X
  TEST_F(CLASSNAME, scale_by_axpy)
    {
    if (typeImplemented_)
      {
      ST alpha = st::zero();
      ST beta  = st::rand();
      PHIST_OUT(9,"axpy, alpha=%f+%f i, beta=%f+%f i",st::real(alpha),
        st::imag(alpha),st::real(beta),st::imag(beta));
      PrintVector(std::cerr,"before scale",
        vec2_vp_,nloc_,lda_,stride_,mpi_comm_);
      SUBR(mvec_add_mvec)(alpha,vec1_,beta,vec2_,&ierr_);
      PrintVector(std::cerr,"after scale",
        vec2_vp_,nloc_,lda_,stride_,mpi_comm_);
      ASSERT_EQ(0,ierr_);
            
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec2_vp_,nloc_,nvec_,lda_,stride_,beta));
      }
    
    }


// TODO - missing tests

// only test the Belos interface for ghost, we didn't write
// the interfaces for Epetra or Tpetra so it is not our problem.
#ifdef PHIST_KERNEL_LIB_GHOST
#ifdef DO_BELOS_TESTS
  // runs all tests from the Belos MvTraits tester
  TEST_F(CLASSNAME, belos_iface)
    {
    if (typeImplemented_)
      {
      Teuchos::RCP<Belos::OutputManager<ST> > MyOM
        = Teuchos::rcp( new Belos::OutputManager<ST>() );
      MyOM->setVerbosity( Belos::Warnings|Belos::Debug);

      ghost_vec_t* v = (ghost_vec_t*)vec1_;
      Teuchos::RCP<phist::GhostMV> ivec = phist::rcp(v,false);

      // test the multivector and its adapter
      bool berr=Belos::TestMultiVecTraits<ST,phist::GhostMV>(MyOM,ivec);
      ASSERT_EQ(true,berr);
      }
    
    }
    
#endif
#endif
