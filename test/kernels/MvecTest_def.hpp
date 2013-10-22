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

      // check that the random function does not change the pointer
      ST* ptr;
      lidx_t lda_new;
      SUBR(mvec_extract_view)(vec1_,&ptr,&lda_new,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_EQ(lda_,lda_new);
      ASSERT_EQ(vec1_vp_,ptr);
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
      
      // check that the random function does not change the pointer
      ST* ptr;
      lidx_t lda_new;
      SUBR(mvec_extract_view)(vec1_,&ptr,&lda_new,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_EQ(lda_,lda_new);
      ASSERT_EQ(vec1_vp_,ptr);
            
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
      ASSERT_EQ(true,minval>mt::eps()); 
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


  // view certain columns, manipulate them and check it changes the
  // correct locations in the original one
  TEST_F(CLASSNAME, view_block)
    {
    if (typeImplemented_)
      {
      int jmin=std::min(2,nvec_-1);
      int jmax=std::min(5,nvec_-1);
      TYPE(mvec_ptr) v1_view=NULL;
      SUBR(mvec_view_block)(vec1_,&v1_view,jmin,jmax,&ierr_);
      ASSERT_EQ(0,ierr_);
      
      _MT_ norms_V1[nvec_];
      _MT_ norms_V1view[jmax-jmin+1];
      
      SUBR(mvec_norm2)(vec1_,norms_V1,&ierr_);
      ASSERT_EQ(0,ierr_);      

      SUBR(mvec_norm2)(v1_view,norms_V1view,&ierr_);
      ASSERT_EQ(0,ierr_);
      
      // compare elements one-by-one because our ArraysEqual expects ST rather than MT as 
      // type here.
      for (int j=jmin;j<=jmax;j++)
        {
        ASSERT_REAL_EQ(norms_V1[j],norms_V1view[j-jmin]);
        }
      // set all the viewed entries to a certain value and check that the original vector is 
      // changed.
      _ST_ val = random_number();
      SUBR(mvec_put_value)(v1_view,val,&ierr_);
      ASSERT_EQ(0,ierr_);
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec1_vp_+jmin*lda_,nloc_,jmax-jmin+1,lda_,stride_,val));
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
