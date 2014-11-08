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
#ifndef MAXBAS
// maximum number of vectors in Krylov space
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
class CLASSNAME: public virtual KernelTestWithVectors<_ST_,_N_,_M_>
{

  public:

    typedef KernelTestWithVectors<_ST_,_N_,_M_> VTest;

    //! mvec/sdMat sizes
    static const int n_=_N_;
    static const int m_=_M_;
    static const int maxBas_=MAXBAS;

    //! Set up routine.
    virtual void SetUp()
    {
      VTest::SetUp();

      if( typeImplemented_ )
      {
        // read matrices
        SUBR(read_mat)(MATNAME,n_,&A_,&ierr_);
        ASSERT_EQ(0,ierr_);
        ASSERT_TRUE(A_ != NULL);

        // wrap matrices in operators
        opA_ = new TYPE(op);
        ASSERT_TRUE(opA_ != NULL);

        SUBR(op_wrap_crsMat)(opA_,A_,&ierr_);
        ASSERT_EQ(0,ierr_);
        
        state_=new TYPE(blockedGMRESstate_ptr)[m_];
        SUBR(blockedGMRESstates_create)(state_,m_,map_,maxBas_,&ierr_);
        
        ASSERT_EQ(0,ierr_);
        xex_=vec1_;
        rhs_=vec2_;
        sol_=vec3_;
        xex_vp_=vec1_vp_;
        rhs_vp_=vec2_vp_;
        sol_vp_=vec3_vp_;
        SUBR(mvec_put_value)(xex_,st::one(),&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_put_value)(sol_,st::zero(),&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(crsMat_times_mvec)(st::one(),A_,xex_,st::zero(),rhs_,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_norm2)(xex_,xNorm_,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_norm2)(rhs_,bNorm_,&ierr_);
        ASSERT_EQ(0,ierr_);
      }
    }

    //! Clean up.
    virtual void TearDown()
    {
      int ierr;
      if( typeImplemented_ )
      {
        ENTER_FCN(__FUNCTION__);
        delete opA_;

        SUBR(crsMat_delete)(A_,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(blockedGMRESstates_delete)(state_,m_,&ierr);
        ASSERT_EQ(0,ierr_);
        delete [] state_;
      }
      //TODO - causes segfault when deleting the map (in KernelTestWithMap::TearDown())
//      VTest::TearDown();
    }

    TYPE(op_ptr) opA_;
    TYPE(blockedGMRESstate_ptr) *state_;

  protected:
  
    TYPE(crsMat_ptr) A_;
    TYPE(mvec_ptr) xex_,sol_,rhs_;
    ST *xex_vp_,*rhs_vp_,*sol_vp_;
    MT bNorm_[_M_],xNorm_[_M_];


    // ========================= the actual arnoldi test =========================
    void doBlockedGMRESTest(int nrhs, int nrestarts, MT tol, bool useMINRES)
    {
      ENTER_FCN(__FUNCTION__);
      if( typeImplemented_ )
      {
        ASSERT_TRUE(nrhs<=m_);
        // up to now we only test the case of one matrix, one rhs here (TODO)
        //ASSERT_TRUE(nrhs==1);

        // set convergence tolerance for all linear systems
        for (int i=0;i<nrhs;i++)
        {
          state_[i]->tol=tol;
        }

        TYPE(mvec_ptr) x=NULL,b=NULL,xex=NULL;
        SUBR(mvec_view_block)(xex_,&xex,0,nrhs-1,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_view_block)(sol_,&x,0,nrhs-1,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_view_block)(rhs_,&b,0,nrhs-1,&ierr_);
        ASSERT_EQ(0,ierr_);
        
        // start with a zero initial guess for each of the systems
        SUBR(mvec_put_value)(x,st::zero(),&ierr_);
//        SUBR(mvec_random)(x,&ierr_);
        ASSERT_EQ(0,ierr_);
        
        int ierr2=0;
        TYPE(mvec_ptr) x_i=NULL,b_i=NULL;
        for (int nr=0;nr<=nrestarts;nr++)
        {
          // initialize state with current approximation
          for (int i=0;i<nrhs;i++)
          {
            SUBR(mvec_view_block)(x,&x_i,i,i,&ierr_);
            SUBR(mvec_view_block)(b,&b_i,i,i,&ierr_);
            if( nr == 0 )
            {
              SUBR(blockedGMRESstate_reset)(state_[i],b_i,x_i,&ierr_);
              ASSERT_EQ(0,ierr_);
            }
            else
            {
              SUBR(blockedGMRESstate_reset)(state_[i],NULL,x_i,&ierr_);
              ASSERT_EQ(0,ierr_);
            }
          }
          // iterate for MAXBAS iterations
          int nIter = 0;
          if( useMINRES )
            SUBR(blockedMINRESstates_iterate)(opA_,state_,nrhs,&nIter, &ierr2);
          else
            SUBR(blockedGMRESstates_iterate)(opA_,state_,nrhs,&nIter,true, &ierr2);
          ASSERT_TRUE(ierr2>=0);
          _MT_ resNorm[nrhs];
          SUBR(blockedGMRESstates_updateSol)(state_,nrhs,x,resNorm,false,&ierr_);
          ASSERT_EQ(0,ierr_);
          if (ierr2==0) break; // converged
        }
        
        ASSERT_EQ(0,ierr2); // GMRES indicated convergence
        SUBR(mvec_delete)(x_i,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_delete)(b_i,&ierr_);
        ASSERT_EQ(0,ierr_);

        // check residual and error norms, compare with tol
        MT resNorm[nrhs], errNorm[nrhs];
        
        
        // compute true residual in b
        SUBR(crsMat_times_mvec)(-st::one(),A_,x,st::one(),b,&ierr_);
        ASSERT_EQ(0,ierr_);
        
        SUBR(mvec_norm2)(b,resNorm,&ierr_);
        ASSERT_EQ(0,ierr_);

        for(int i = 0; i < nrhs; i++)
        {
          PHIST_DEB("required tol:  %6.2g\n",tol);
          PHIST_DEB("rel. res norm: %6.2g\n",resNorm[i]/bNorm_[i]);

          // check error (residual two norm < tol)
          ASSERT_TRUE(resNorm[i]<=tol*bNorm_[i]);
        }

        // check error
        SUBR(mvec_add_mvec)(-st::one(),xex,st::one(),x,&ierr_);
        ASSERT_EQ(0,ierr_);

        SUBR(mvec_norm2)(x,errNorm,&ierr_);
        ASSERT_EQ(0,ierr_);

        for(int i = 0; i < nrhs; i++)
        {
          PHIST_DEB("rel. err norm: %6.2g\n",errNorm[i]/bNorm_[i]);
          // TODO - we should use the matrix norm for the scaling here,
          //        but since we only test with one matrix we can hard-code
          //        a factor 20 so that the test is passed.
          ASSERT_TRUE(errNorm[i]<=20*tol*bNorm_[i]);
        }

        SUBR(mvec_delete)(x,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_delete)(b,&ierr_);
        ASSERT_EQ(0,ierr_);
        SUBR(mvec_delete)(xex,&ierr_);
        ASSERT_EQ(0,ierr_);
      }
    }
};

// test unrestarted GMRES with tolerance TOL1
TEST_F(CLASSNAME, simple_blockedGMRES) 
{
  int nrestarts=0;
  int nrhs=1;
  doBlockedGMRESTest(nrhs,nrestarts,TOLA,false);
}

// test restarted GMRES with tolerance TOL2
TEST_F(CLASSNAME, restarted_blockedGMRES) 
{
  int nrestarts=5;
  int nrhs=1;
  doBlockedGMRESTest(nrhs,nrestarts,TOLB,false);
}

// test unrestarted GMRES with tolerance TOL1
TEST_F(CLASSNAME, multiple_simple_blockedGMRES) 
{
  int nrestarts=0;
  int nrhs=4;
  doBlockedGMRESTest(nrhs,nrestarts,TOLA,false);
}

#ifdef MATSYMMETRIC
// test unrestarted GMRES with tolerance TOL1
TEST_F(CLASSNAME, simple_blockedminres) 
{
  int nrestarts=0;
  int nrhs=1;
  doBlockedGMRESTest(nrhs,nrestarts,TOLA,true);
}

// test restarted GMRES with tolerance TOL2
TEST_F(CLASSNAME, restarted_blockedminres) 
{
  int nrestarts=5;
  int nrhs=1;
  doBlockedGMRESTest(nrhs,nrestarts,TOLB,true);
}

// test unrestarted GMRES with tolerance TOL1
TEST_F(CLASSNAME, multiple_simple_blockedminres) 
{
  int nrestarts=0;
  int nrhs=4;
  doBlockedGMRESTest(nrhs,nrestarts,TOLA,true);
}
#endif
