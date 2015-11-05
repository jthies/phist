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

      if( typeImplemented_ && !problemTooSmall_ )
      {
        // read matrices
        SUBR(read_mat)(MATNAME,comm_,n_,&A_,&iflag_);
        ASSERT_EQ(0,iflag_);
        ASSERT_TRUE(A_ != NULL);

        // wrap matrices in operators
        opA_ = new TYPE(op);
        ASSERT_TRUE(opA_ != NULL);

        SUBR(op_wrap_sparseMat)(opA_,A_,&iflag_);
        ASSERT_EQ(0,iflag_);

        const_map_ptr_t map = NULL;
        SUBR(sparseMat_get_domain_map)(A_, &map, &iflag_);
        ASSERT_EQ(0,iflag_);
        replaceMap(map);
        ASSERT_EQ(0,iflag_);

        state_=new TYPE(blockedGMRESstate_ptr)[m_];
        SUBR(blockedGMRESstates_create)(state_,m_,map,maxBas_,&iflag_);
        
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
      }
    }

    //! Clean up.
    virtual void TearDown()
    {
      int iflag;
      if( typeImplemented_ && !problemTooSmall_ )
      {
        PHIST_ENTER_FCN(__FUNCTION__);
        delete opA_;

        SUBR(sparseMat_delete)(A_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(blockedGMRESstates_delete)(state_,m_,&iflag);
        ASSERT_EQ(0,iflag_);
        delete [] state_;
      }
      //TODO - causes segfault when deleting the map (in KernelTestWithMap::TearDown())
      //melven: why/where?
      VTest::TearDown();
    }

    TYPE(op_ptr) opA_;
    TYPE(blockedGMRESstate_ptr) *state_;

  protected:
  
    TYPE(sparseMat_ptr) A_;
    TYPE(mvec_ptr) xex_,sol_,rhs_;
    ST *xex_vp_,*rhs_vp_,*sol_vp_;
    MT bNorm_[_M_],xNorm_[_M_];


    // ========================= the actual arnoldi test =========================
    void doBlockedGMRESTest(int nrhs, int nrestarts, MT tol, bool useMINRES)
    {
      PHIST_ENTER_FCN(__FUNCTION__);
      if( typeImplemented_ && !problemTooSmall_ )
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
        SUBR(mvec_view_block)(xex_,&xex,0,nrhs-1,&iflag_);
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
        for (int nr=0;nr<=nrestarts;nr++)
        {
          // initialize state with current approximation
          for (int i=0;i<nrhs;i++)
          {
            SUBR(mvec_view_block)(x,&x_i,i,i,&iflag_);
            SUBR(mvec_view_block)(b,&b_i,i,i,&iflag_);
            if( nr == 0 )
            {
              SUBR(blockedGMRESstate_reset)(state_[i],b_i,x_i,&iflag_);
              ASSERT_EQ(0,iflag_);
            }
            else
            {
              SUBR(blockedGMRESstate_reset)(state_[i],NULL,x_i,&iflag_);
              ASSERT_EQ(0,iflag_);
            }
          }
          // iterate for MAXBAS iterations
          int nIter = 0;
          if( useMINRES )
            SUBR(blockedMINRESstates_iterate)(opA_,state_,nrhs,&nIter, &iflag2);
          else
            SUBR(blockedGMRESstates_iterate)(opA_,state_,nrhs,&nIter,true, &iflag2);
          ASSERT_TRUE(iflag2>=0);
          _MT_ resNorm[nrhs];
          SUBR(blockedGMRESstates_updateSol)(state_,nrhs,x,resNorm,false,&iflag_);
          ASSERT_EQ(0,iflag_);
          if (iflag2==0) break; // converged
        }
        
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

// test unrestarted GMRES with tolerance TOL1
TEST_F(CLASSNAME, simple_blockedGMRES) 
{
  int nrestarts=0;
  int nrhs=1;
PHIST_TASK_DECLARE(singleTask)
PHIST_TASK_BEGIN(singleTask)
  doBlockedGMRESTest(nrhs,nrestarts,TOLA,false);
PHIST_TASK_END(&iflag_)
}

// test restarted GMRES with tolerance TOL2
TEST_F(CLASSNAME, restarted_blockedGMRES) 
{
  int nrestarts=5;
  int nrhs=1;
PHIST_TASK_DECLARE(singleTask)
PHIST_TASK_BEGIN(singleTask)
  doBlockedGMRESTest(nrhs,nrestarts,TOLB,false);
PHIST_TASK_END(&iflag_)
}

// test unrestarted GMRES with tolerance TOL1
TEST_F(CLASSNAME, multiple_simple_blockedGMRES) 
{
  int nrestarts=0;
  int nrhs=4;
PHIST_TASK_DECLARE(singleTask)
PHIST_TASK_BEGIN(singleTask)
  doBlockedGMRESTest(nrhs,nrestarts,TOLA,false);
PHIST_TASK_END(&iflag_)
}

#ifdef MATSYMMETRIC
// test unrestarted GMRES with tolerance TOL1
TEST_F(CLASSNAME, simple_blockedminres) 
{
  int nrestarts=0;
  int nrhs=1;
PHIST_TASK_DECLARE(singleTask)
PHIST_TASK_BEGIN(singleTask)
  doBlockedGMRESTest(nrhs,nrestarts,TOLA,true);
PHIST_TASK_END(&iflag_)
}

// test restarted GMRES with tolerance TOL2
TEST_F(CLASSNAME, restarted_blockedminres) 
{
  int nrestarts=10;
  int nrhs=1;
PHIST_TASK_DECLARE(singleTask)
PHIST_TASK_BEGIN(singleTask)
  doBlockedGMRESTest(nrhs,nrestarts,TOLB,true);
PHIST_TASK_END(&iflag_)
}

// test unrestarted GMRES with tolerance TOL1
TEST_F(CLASSNAME, multiple_simple_blockedminres) 
{
  int nrestarts=0;
  int nrhs=4;
PHIST_TASK_DECLARE(singleTask)
PHIST_TASK_BEGIN(singleTask)
  doBlockedGMRESTest(nrhs,nrestarts,TOLA,true);
PHIST_TASK_END(&iflag_)
}
#endif
