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
#ifdef PHIST_KERNEL_LIB_FORTRAN
          vec1_vp_[j+i*lda_]=random_number();
          vec2_vp_[j+i*lda_]=st::one();
#else
          vec1_vp_[j*lda_+i]=random_number();
          vec2_vp_[j*lda_+i]=st::one();
#endif
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
      ST val = (_ST_)42.0 + (ST)3.0*st::cmplx_I();
      SUBR(mvec_put_value)(vec1_,val,&ierr_);
      ASSERT_EQ(0,ierr_);
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec1_vp_,nvec_,nloc_,lda_,stride_,val));
#else
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec1_vp_,nloc_,nvec_,lda_,stride_,val));
#endif

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
#ifdef PHIST_KERNEL_LIB_FORTRAN
        vec2_vp_[j+i*lda_]=mt::one()/st::conj(vec1_vp_[j+i*lda_]);
#else
        vec2_vp_[j*lda_+i]=mt::one()/st::conj(vec1_vp_[j*lda_+i]);
#endif
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
#ifdef PHIST_KERNEL_LIB_FORTRAN
          absval[k++]=st::abs(vec1_vp_[j+i*lda_]);
#else
          absval[k++]=st::abs(vec1_vp_[j*lda_+i]);
#endif
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

  // 2-norm, nrm2=sqrt(v'v)
  TEST_F(CLASSNAME, norm2)
    {
    if (typeImplemented_)
      {
      int ilower;     
      phist_map_get_ilower(map_,&ilower,&ierr_);
      ASSERT_EQ(0,ierr_);
      for (int j=0;j<nvec_;j++)
        {
        for (int i=0;i<nloc_*stride_;i+=stride_)
          {
#ifdef PHIST_KERNEL_LIB_FORTRAN
          vec1_vp_[j+i*lda_]=ilower+i;
#else
          vec1_vp_[j*lda_+i]=ilower+i;
#endif
          }
        }
      MT expect = 0.0;
      for (int i=0;i<_N_;i++)
        {
        expect+=(MT)(i*i);
        }
      expect=mt::sqrt(expect);

      MT nrm2[nvec_];
      SUBR(mvec_norm2)(vec1_,nrm2,&ierr_);
      ASSERT_EQ(0,ierr_);
      for (int i=0;i<nvec_;i++)
        {
        ASSERT_REAL_EQ(expect,nrm2[i]);
        }
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

#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nvec_,nloc_,lda_,stride_));
#else
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nloc_,nvec_,lda_,stride_));
#endif
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
            
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec2_vp_,nvec_,nloc_,lda_,stride_,beta));
#else
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec2_vp_,nloc_,nvec_,lda_,stride_,beta));
#endif
      }
    
    }

  // X = Y*diag(a_1,...,a_nvec) + a*X
  TEST_F(CLASSNAME, random_vadd)
  {
    if( typeImplemented_ )
    {
      ST beta = st::rand();
      ST alpha[_NV_];
      for(int i = 0; i < nvec_; i++)
        alpha[i] = st::rand();

      SUBR(mvec_vadd_mvec)(alpha,vec1_,beta,vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);

      // calculate solution by hand
      for(int i = 0; i < nloc_; i++)
        for(int j = 0; j < nvec_; j++)
#ifdef PHIST_KERNEL_LIB_FORTRAN
          vec1_vp_[j+i*lda_] = alpha[j]*vec1_vp_[j+i*lda_]+beta;
#else
          vec1_vp_[j*lda_+i] = alpha[j]*vec1_vp_[j*lda_+i]+beta;
#endif
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nvec_,nloc_,lda_,stride_));
#else
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nloc_,nvec_,lda_,stride_));
#endif
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
      // create a view of the view
      TYPE(mvec_ptr) v1_vv=NULL;
      SUBR(mvec_view_block)(v1_view,&v1_vv,0,jmax-jmin,&ierr_);
      ASSERT_EQ(0,ierr_);
      
      // now this should delete the original view and create a new one,
      // all vectors must remain valid:
      SUBR(mvec_view_block)(vec1_,&v1_view,jmin,jmax,&ierr_);
      
      _MT_ norms_V1[nvec_];
      _MT_ norms_V1view[nvec_];
      _MT_ norms_V1vv[nvec_];
      
      SUBR(mvec_norm2)(vec1_,norms_V1,&ierr_);
      ASSERT_EQ(0,ierr_);      

      SUBR(mvec_norm2)(v1_view,norms_V1view,&ierr_);
      ASSERT_EQ(0,ierr_);

      SUBR(mvec_norm2)(v1_vv,norms_V1vv,&ierr_);
      ASSERT_EQ(0,ierr_);
      
      // compare elements one-by-one because our ArraysEqual expects ST rather than MT as 
      // type here.
      for (int j=jmin;j<=jmax;j++)
        {
        ASSERT_REAL_EQ(norms_V1[j],norms_V1view[j-jmin]);
        ASSERT_REAL_EQ(norms_V1[j],norms_V1vv[j-jmin]);
        }
      // set all the viewed entries to a certain value and check that the original vector is 
      // changed.
      _ST_ val = random_number();
      SUBR(mvec_put_value)(v1_view,val,&ierr_);
      ASSERT_EQ(0,ierr_);
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec1_vp_+jmin,jmax-jmin+1,nloc_,lda_,stride_,val));
#else
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec1_vp_+jmin*lda_,nloc_,jmax-jmin+1,lda_,stride_,val));
#endif

      // new norms after changing columns
      SUBR(mvec_norm2)(vec1_,norms_V1,&ierr_);
      ASSERT_EQ(0,ierr_);      
      
      // delete the views without harming v1
      SUBR(mvec_delete)(v1_view,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_delete)(v1_vv,&ierr_);
      ASSERT_EQ(0,ierr_);
      // check that v1 still works
      SUBR(mvec_norm2)(vec1_,norms_V1vv,&ierr_);
      ASSERT_EQ(0,ierr_);
      for (int j=0;j<nvec_;j++)
        {
        ASSERT_REAL_EQ(norms_V1[j],norms_V1vv[j]);
        }
      }
    }


  // view some columns as a new mvec, compare ||V|| calculations
  // and check that modifying the view vector modifies the original ones
  TEST_F(CLASSNAME, view_scattered)
  {
    if( typeImplemented_ )
    {
      const int nc=4;
      int idx[nc]={3,0,2,5};
      for (int i=0;i<nc;i++)
      {
        idx[i]=std::min(idx[i],nvec_-1);
      }
      MT norms0[nvec_];
      MT norms1[nc];
      SUBR(mvec_norm2)(vec1_,norms0,&ierr_);
      ASSERT_EQ(0,ierr_);

      TYPE(mvec_ptr) V=NULL;
      SUBR(mvec_view_scattered)(vec1_,&V,idx,nc,&ierr_);
      ASSERT_EQ(0,ierr_);
      SUBR(mvec_norm2)(V,norms1,&ierr_);
      ASSERT_EQ(0,ierr_);

      for (int i=0;i<nc;i++)
      {
        ASSERT_REAL_EQ(norms0[idx[i]],norms1[i]);
        ASSERT_TRUE(norms0[idx[i]]!=mt::zero());
      }
        
      // now randomize the view and check again
      SUBR(mvec_random)(V,&ierr_);
      ASSERT_EQ(0,ierr_);

      SUBR(mvec_norm2)(vec1_,norms0,&ierr_);
      ASSERT_EQ(0,ierr_);

      SUBR(mvec_norm2)(V,norms1,&ierr_);
      ASSERT_EQ(0,ierr_);

      for (int i=0;i<nc;i++)
      {
        ASSERT_REAL_EQ(norms0[idx[i]],norms1[i]);
      }
    }
  }

  // copy in and out columns
  TEST_F(CLASSNAME, get_set_block)
    {
    if (typeImplemented_)
      {
      int jmin=std::min(2,nvec_-1);
      int jmax=std::min(5,nvec_-1);
      TYPE(mvec_ptr) v1_copy=NULL;
      SUBR(mvec_create)(&v1_copy,map_,jmax-jmin+1,&ierr_);
      ASSERT_EQ(0,ierr_);
      
      SUBR(mvec_get_block)(vec1_,v1_copy,jmin,jmax,&ierr_);
      ASSERT_EQ(0,ierr_);

      _MT_ norms_V1[nvec_];
      _MT_ norms_V1copy[nvec_];
      
      SUBR(mvec_norm2)(vec1_,norms_V1,&ierr_);
      ASSERT_EQ(0,ierr_);      

      SUBR(mvec_norm2)(v1_copy,norms_V1copy,&ierr_);
      ASSERT_EQ(0,ierr_);

      // compare elements one-by-one because our ArraysEqual expects ST rather than MT as 
      // type here.
      for (int j=jmin;j<=jmax;j++)
        {
        ASSERT_REAL_EQ(norms_V1[j],norms_V1copy[j-jmin]);
        }
      // set all the viewed entries to a certain value and check that the original vector is 
      // changed.
      _ST_ val = random_number();
      SUBR(mvec_put_value)(v1_copy,val,&ierr_);
      ASSERT_EQ(0,ierr_);
      
      // check that the norms of v1 are unchanged
      SUBR(mvec_norm2)(vec1_,norms_V1copy,&ierr_);
      ASSERT_EQ(0,ierr_);

      for (int j=0;j<nvec_;j++)
        {
        ASSERT_REAL_EQ(norms_V1[j],norms_V1copy[j]);
        }

      // compute the new norms
      // check that the norms of v1 are unchanged
      SUBR(mvec_norm2)(v1_copy,norms_V1copy,&ierr_);
      ASSERT_EQ(0,ierr_);
      
      // now set the block as the corresponding columns of v2
      SUBR(mvec_set_block)(vec2_,v1_copy,jmin,jmax,&ierr_);
      ASSERT_EQ(0,ierr_);

      // compute the new norms
      // check that the norms of v1 are unchanged
      SUBR(mvec_norm2)(vec2_,norms_V1,&ierr_);
      ASSERT_EQ(0,ierr_);

      for (int j=jmin;j<=jmax;j++)
        {
        ASSERT_REAL_EQ(norms_V1[j],norms_V1copy[j-jmin]);
        }
      }
    }

  TEST_F(CLASSNAME, scale)
  {
    if( typeImplemented_ )
    {
      _ST_ scale = st::rand();

      SUBR(mvec_scale)(vec2_,scale,&ierr_);
      ASSERT_EQ(0,ierr_);
#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec2_vp_,nvec_,nloc_,lda_,stride_,scale));
#else
      ASSERT_REAL_EQ(mt::one(),ArrayEqual(vec2_vp_,nloc_,nvec_,lda_,stride_,scale));
#endif

      SUBR(mvec_add_mvec)(st::one(),vec1_,st::zero(),vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);

      scale = st::rand();
      SUBR(mvec_scale)(vec1_,scale,&ierr_);
      ASSERT_EQ(0,ierr_);
      // apply scale to vec2_ by hand
      for(int i = 0; i < nloc_; i++)
        for(int j = 0; j < nvec_; j++)
#ifdef PHIST_KERNEL_LIB_FORTRAN
          vec2_vp_[j+i*lda_] *= scale;
#else
          vec2_vp_[j*lda_+i] *= scale;
#endif

#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nvec_,nloc_,lda_,stride_));
#else
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nloc_,nvec_,lda_,stride_));
#endif
    }
  }

  TEST_F(CLASSNAME, vscale)
  {
    if( typeImplemented_ )
    {
      _ST_ scale[_NV_];
      for(int i = 0; i < _NV_; i++)
        scale[i] = st::rand();

      SUBR(mvec_add_mvec)(st::one(),vec1_,st::zero(),vec2_,&ierr_);
      ASSERT_EQ(0,ierr_);

      SUBR(mvec_vscale)(vec1_,scale,&ierr_);
      ASSERT_EQ(0,ierr_);
      // apply scale to vec2_ by hand
      for(int i = 0; i < nloc_; i++)
        for(int j = 0; j < nvec_; j++)
#ifdef PHIST_KERNEL_LIB_FORTRAN
          vec2_vp_[j+i*lda_] *= scale[j];
#else
          vec2_vp_[j*lda_+i] *= scale[j];
#endif

#ifdef PHIST_KERNEL_LIB_FORTRAN
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nvec_,nloc_,lda_,stride_));
#else
      ASSERT_REAL_EQ(mt::one(),ArraysEqual(vec1_vp_,vec2_vp_,nloc_,nvec_,lda_,stride_));
#endif
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
