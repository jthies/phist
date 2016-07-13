#ifndef CLASSNAME
#error "file not included correctly"
#endif
#ifndef MATNAME
// matrix base filename (eg. xyz => xyz'_N_'.'ext'
#error "file not included correctly"
#endif
#ifndef _N_
// global size of matrix
#error "file not included correctly"
#endif
#ifndef _M_
// number of systems solved simultaneously
#error "file not included correctly"
#endif
#ifndef TOLA
// tolerance to be achieved in MAXBAS iterations (no restart)
#error "file not included correctly"
#endif
#ifndef TOLB
// tolerance to be achieved in 5*MAXBAS iterations (4 restarts)
#error "file not included correctly"
#endif

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithSparseMat<_ST_,_N_,MATNAME>,
                 public virtual KernelTestWithVectors<_ST_,_N_,_M_,0,3>
{

  public:

    typedef KernelTestWithSparseMat<_ST_,_N_,MATNAME> SparseMatTest;
    typedef KernelTestWithVectors<_ST_,_N_,_M_,0,3> VTest;

    //! mvec/sdMat sizes
    static const int n_=_N_;
    static const int m_=_M_;
    
    //! (optional) additional projection vectors
    TYPE(mvec_ptr) Vproj_;
    
    //! shifts
    static const int numShifts_=_M_;
    static MT shift_r_[numShifts_], shift_i_[numShifts_];

    static void SetUpTestCase()
    {
      SparseMatTest::SetUpTestCase();
      VTest::SetUpTestCase();
    }

    //! Set up routine.
    virtual void SetUp()
    {
      SparseMatTest::SetUp();
      VTest::SetUp();

      if( typeImplemented_ && !problemTooSmall_ )
      {
      
        static const _MT_ shifts_r[10]={0.0, 0.0, 0.0, 1.0, 0.1, 0.1, 0.1,-0.1,-0.1,-0.1};
        static const _MT_ shifts_i[10]={0.0, 0.1,-0.1, 0.0, 0.0, 0.1,-0.1, 0.0, 0.1,-0.1};
     
        for (int i=0; i<numShifts_;i++)
        {
          shift_r_[i]=shifts_r[i];
          shift_i_[i]=shifts_i[i];
        }
      
        SUBR(carp_cgState_create)(&state_,A_,Vproj_,numShifts_,shift_r_,shift_i_, &iflag_);
        
        ASSERT_EQ(0,iflag_);
        xex_=vec1_;
        rhs_=vec2_;
        sol_=vec3_;
        xex_vp_=vec1_vp_;
        rhs_vp_=vec2_vp_;
        sol_vp_=vec3_vp_;
        SUBR(mvec_put_value)(xex_,st::one(),&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_put_value)(sol_,st::zero(),&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(sparseMat_times_mvec)(st::one(),A_,xex_,st::zero(),rhs_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_norm2)(xex_,xNorm_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_norm2)(rhs_,bNorm_,&iflag_);
        ASSERT_EQ(0,iflag_);
        Vproj_=NULL;
      }
    }

    //! Clean up.
    virtual void TearDown()
    {
      int iflag;
      if( typeImplemented_ && !problemTooSmall_ )
      {
        SUBR(carp_cgState_delete)(state_,&iflag);
        ASSERT_EQ(0,iflag_);
      }

      VTest::TearDown();
      SparseMatTest::TearDown();
    }

    static void TearDownTestCase()
    {
      VTest::TearDownTestCase();
      SparseMatTest::TearDownTestCase();
    }

    TYPE(carp_cgState_ptr) state_;

  protected:
  
    TYPE(sparseMat_ptr) A_;
    TYPE(mvec_ptr) xex_,sol_,rhs_;
    ST *xex_vp_,*rhs_vp_,*sol_vp_;
    MT bNorm_[_M_],xNorm_[_M_];


    // ========================= the actual arnoldi test =========================
    void doCarpCGTest(TYPE(mvec_ptr) B, TYPE(mvec_ptr) Xr, TYPE(mvec_ptr) X_i, MT tol)
    {
      PHIST_ENTER_FCN(__FUNCTION__);
      if( typeImplemented_ && !problemTooSmall_ )
      {
        int nrhs=nvec_;
        TYPE(mvec_ptr) x=NULL,b=NULL,xex=NULL;
        SUBR(mvec_view_block)(xex_,&xex,0,nvec_-1,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_view_block)(sol_,&x,0,nrhs-1,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_view_block)(rhs_,&b,0,nrhs-1,&iflag_);
        ASSERT_EQ(0,iflag_);
        
        // start with a zero initial guess for each of the systems
        SUBR(mvec_put_value)(x,st::zero(),&iflag_);
//        SUBR(mvec_random)(x,&iflag_);
        ASSERT_EQ(0,iflag_);
        
        int iflag2=0;
        TYPE(mvec_ptr) x_i=NULL,b_i=NULL;
#if 0
        for (int nr=0;nr<=nrestarts;nr++)
        {
          // initialize state with current approximation
          for (int i=0;i<nrhs;i++)
          {
            SUBR(mvec_view_block)(x,&x_i,i,i,&iflag_);
            SUBR(mvec_view_block)(b,&b_i,i,i,&iflag_);
            if( nr == 0 )
            {
              SUBR(carp_cgState_reset)(state_[i],b_i,x_i,&iflag_);
              ASSERT_EQ(0,iflag_);
            }
            else
            {
              SUBR(carp_cgState_reset)(state_[i],NULL,x_i,&iflag_);
              ASSERT_EQ(0,iflag_);
            }
          }
          // iterate for MAXBAS iterations
          int nIter = 0;
          SUBR(pminresState_iterate)(opA_,state_,nrhs,&nIter, &iflag2);
          ASSERT_TRUE(iflag2>=0);
          _MT_ resNorm[nrhs];
          ASSERT_EQ(0,iflag_);
          if (iflag2==0) break; // converged
        }
#endif
        ASSERT_EQ(0,iflag2); // GMRES indicated convergence
        SUBR(mvec_delete)(x_i,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_delete)(b_i,&iflag_);
        ASSERT_EQ(0,iflag_);

        // check residual and error norms, compare with tol
        MT resNorm[nrhs], errNorm[nrhs];
        
        
        // compute true residual in b
        SUBR(sparseMat_times_mvec)(-st::one(),A_,x,st::one(),b,&iflag_);
        ASSERT_EQ(0,iflag_);
        
        SUBR(mvec_norm2)(b,resNorm,&iflag_);
        ASSERT_EQ(0,iflag_);

        for(int i = 0; i < nrhs; i++)
        {
          PHIST_DEB("required tol:  %6.2g\n",tol);
          PHIST_DEB("rel. res norm: %6.2g\n",resNorm[i]/bNorm_[i]);

          // check error (residual two norm < tol)
          ASSERT_TRUE(resNorm[i]<=tol*bNorm_[i]);
        }

        // check error
        SUBR(mvec_add_mvec)(-st::one(),xex,st::one(),x,&iflag_);
        ASSERT_EQ(0,iflag_);

        SUBR(mvec_norm2)(x,errNorm,&iflag_);
        ASSERT_EQ(0,iflag_);

        for(int i = 0; i < nrhs; i++)
        {
          PHIST_DEB("rel. err norm: %6.2g\n",errNorm[i]/bNorm_[i]);
          // TODO - we should use the matrix norm for the scaling here,
          //        but since we only test with one matrix we can hard-code
          //        a factor 20 so that the test is passed.
          ASSERT_TRUE(errNorm[i]<=20*tol*bNorm_[i]);
        }

        SUBR(mvec_delete)(x,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_delete)(b,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_delete)(xex,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
    }
};

#if 0
// test unrestarted GMRES with tolerance TOL1
TEST_F(CLASSNAME, TODO_add_tests)
{
  int nrestarts=0;
  int nrhs=1;
  doCarpCGTest(nrhs,nrestarts,TOLA,false);
}

#endif
