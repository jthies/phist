/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "../tools/TestHelpers.h"
#ifndef CLASSNAME
#error "file not included correctly."
#endif

/*! Test fixure. */
class CLASSNAME: public virtual KernelTestWithSparseMat<_ST_,_N_,_N_,MATNAME>,
                 public virtual KernelTestWithSparseMat<_ST_,_N_,_N_,PRECNAME>,
                 public virtual KernelTestWithVectors<_ST_,_N_,_NV_,0,4>,
                 public virtual KernelTestWithVectors<_ST_,_N_,_NVP_,0,2>,
                 public virtual KernelTestWithSdMats<_ST_,_NVP_,_NV_>,
                 public JadaTestWithOpts
{

  public:
    typedef KernelTestWithSparseMat<_ST_,_N_,_N_,MATNAME> ATest;
    typedef KernelTestWithSparseMat<_ST_,_N_,_N_,PRECNAME> PTest;
    typedef KernelTestWithVectors<_ST_,_N_,_NV_,0,4> VTest;
    typedef KernelTestWithVectors<_ST_,_N_,_NVP_,0,2> QTest;
    typedef KernelTestWithSdMats<_ST_,_NVP_,_NV_> MTest;

    static void SetUpTestCase()
    {
      int sparseMatCreateFlag=getSparseMatCreateFlag(_N_,_NV_);
      ATest::SetUpTestCase(sparseMatCreateFlag);
      PTest::SetUpTestCase(sparseMatCreateFlag,ATest::context_);
      VTest::SetUpTestCase();
      QTest::SetUpTestCase();
      MTest::SetUpTestCase();
    }

    /*! Set up routine.
    */
    virtual void SetUp()
    {
      ATest::SetUp();
      PTest::SetUp();
      VTest::SetUp();
      QTest::SetUp();
      MTest::SetUp();
      JadaTestWithOpts::SetUp();

      jadaOpts_.innerSolvType=_SOLVTYPE_;
      jadaOpts_.innerSolvMaxIters=_MAXINNER_;
      jadaOpts_.innerSolvMaxBas=_MAXINNER_;
      jadaOpts_.blockSize=_NV_;
      jadaOpts_.minBas=12;
      jadaOpts_.maxBas=jadaOpts_.minBas+10*_NV_;
      jadaOpts_.convTol=1e-6;

      if( typeImplemented_ && !problemTooSmall_ )
      {
        // disable the test because TSQR will not work.
        // This situation occurs if we have a small matrix (_N_=25, say)
        // and many Q vectors (e.g. 10) with multiple MPI procs.
        int globalTooSmall = _N_ <= std::min(_NVP_,_NV_);
#ifdef PHIST_HAVE_MPI
        int localTooSmall = nloc_ <= std::min(_NVP_,_NV_);
        iflag_ = MPI_Allreduce(&localTooSmall, &globalTooSmall, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        ASSERT_EQ(0,iflag_);
#endif
        problemTooSmall_ = globalTooSmall != 0;
      }

      if (typeImplemented_ && !problemTooSmall_)
      {
        opA_ = new TYPE(linearOp);
        SUBR(linearOp_wrap_sparseMat)(opA_, ATest::A_, &iflag_);
        ASSERT_EQ(0,iflag_);

        opP_ = new TYPE(linearOp);
        SUBR(linearOp_wrap_sparseMat)(opP_, PTest::A_, &iflag_);
        ASSERT_EQ(0,iflag_);
        
        v1_=VTest::vec1_;
        v2_=VTest::vec2_;
        v3_=VTest::vec3_;
        v4_=VTest::vec4_;
        q_=QTest::vec1_;
        sigma_ = new _ST_[_NV_];
        noSigma_ = new _ST_[std::max(_NV_,_NVP_)];
        negSigma_ = new _ST_[_NV_];
        for(int i = 0; i < _NV_; i++)
        {
          // there are hopefully no eigenvalues in this region so the matrix doesn't get nearly singular
          sigma_[i] = (_ST_)30*st::one() + (_ST_)5*st::prand();
          negSigma_[i] = -sigma_[i];
          noSigma_[i]  = st::zero();
        }
        for (int i=0; i<_NVP_; i++) noSigma_[i]=st::zero();

        // create random orthogonal Q
        SUBR(mvec_random)(q_,&iflag_);
        ASSERT_EQ(0,iflag_);
        TYPE(sdMat_ptr) Rtmp;
        SUBR(sdMat_create)(&Rtmp,_NVP_,_NVP_,comm_,&iflag_);
        ASSERT_EQ(0,iflag_);
        int rankQ = 0;
        SUBR(orthog)(NULL,q_,NULL,Rtmp,NULL,4,&rankQ,&iflag_);
        ASSERT_GE(iflag_,0);
        SUBR(sdMat_delete)(Rtmp,&iflag_);
        ASSERT_EQ(0,iflag_);

        jdOp_ = new TYPE(linearOp);
        SUBR(jadaOp_create)(opA_,NULL,q_,NULL,negSigma_,_NV_,jdOp_,&iflag_);
        ASSERT_EQ(0,iflag_);

        //note: we also want to apply this operator to Q, hence the max(_NV_,_NVP_) argument for the number of shifts
        jdPrec_ = new TYPE(linearOp);
        SUBR(jadaPrec_create)(opP_,q_,q_,noSigma_,std::max(_NV_,_NVP_),jdPrec_,1,&iflag_);
        ASSERT_EQ(0,iflag_);

        // setup system to solve, exact x and A*x
        SUBR(mvec_random)(v2_,&iflag_);
        ASSERT_EQ(0,iflag_);
        // x needs to be in q^orth
        SUBR(mvecT_times_mvec)(st::one(),q_,v2_,st::zero(),mat1_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_times_sdMat)(-st::one(),q_,mat1_,st::one(),v2_,&iflag_);
        ASSERT_EQ(0,iflag_);
        _MT_ tmp[_NV_];
        SUBR(mvec_normalize)(v2_, tmp, &iflag_);
        ASSERT_EQ(0,iflag_);
        jdOp_->apply(st::one(),jdOp_->A,v2_,st::zero(),v3_,&iflag_);
        ASSERT_EQ(0,iflag_);
        // as we are only interested in the direction of vec3_, scale it to one
        SUBR(mvec_normalize)(v3_, tmp, &iflag_);
        ASSERT_EQ(0,iflag_);

        jadaOpts_.preconType=phist_USER_PRECON;
        jadaOpts_.preconOp=opP_;
      }
    }

    /*! Clean up.
    */
    virtual void TearDown() 
    {
      if (typeImplemented_ && !problemTooSmall_)
      {
        SUBR(jadaOp_delete)(jdOp_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(jadaOp_delete)(jdPrec_,&iflag_);
        ASSERT_EQ(0,iflag_);
        if( jdOp_ != NULL )
          delete jdOp_;
        jdOp_ = NULL;
        if( jdPrec_ != NULL )
          delete jdPrec_;
        jdOp_ = NULL;
        if( opA_ != NULL )
          delete opA_;
        opA_ = NULL;
        if( opP_ != NULL )
          delete opP_;
        opP_ = NULL;
        if( sigma_ != NULL )
          delete[] sigma_;
        sigma_ = NULL;
        if( negSigma_ != NULL )
          delete[] negSigma_;
        negSigma_ = NULL;
        if( noSigma_ != NULL )
          delete[] noSigma_;
        noSigma_ = NULL;
      }
      MTest::TearDown();
      VTest::TearDown();
      QTest::TearDown();
      ATest::TearDown();
      PTest::TearDown();
    }

    static void TearDownTestCase()
    {
      MTest::TearDownTestCase();
      VTest::TearDownTestCase();
      QTest::TearDownTestCase();
      ATest::TearDownTestCase();
      PTest::TearDownTestCase();
    }

    void checkResiduals(_MT_ tol[_NV_])
    {
      // the calculated approximations should be stored in vec2, the JD-residual in vec3_

      // the solution should be scaled to one
      _MT_ solutionNorm[_NV_];
      SUBR(mvec_norm2)(v2_, solutionNorm, &iflag_);
      ASSERT_EQ(0,iflag_);
      for(int i = 0; i < _NV_; i++)
      {
        ASSERT_NEAR(mt::one(), solutionNorm[i], 100*VTest::releps());
      }

      // we cannot directly compare the vectors because the solution is scaled
      jdOp_->apply(st::one(),jdOp_->A,v2_,st::zero(),v1_,&iflag_);
      ASSERT_EQ(0,iflag_);
      // rather check that v1_ and vec3_ point in the same direction
      _MT_ tmp[_NV_];
      SUBR(mvec_normalize)(v1_, tmp, &iflag_);
      ASSERT_EQ(0,iflag_);
      _ST_ dot[_NV_];
      SUBR(mvec_dot_mvec)(v3_, v1_, dot, &iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_vadd_mvec)(dot, v3_, -st::one(), v1_, &iflag_);
      ASSERT_EQ(0,iflag_);
      _MT_ resNorm[_NV_];
      SUBR(mvec_norm2)(v1_, resNorm, &iflag_);
      ASSERT_EQ(0,iflag_);
      // check that the resnorm is actually near the required norm
      PHIST_SOUT(PHIST_INFO, "required tol.:");
      for(int i = 0; i < _NV_; i++)
      {
        PHIST_SOUT(PHIST_INFO, "\t%8.4e", tol[i]);
      }
      PHIST_SOUT(PHIST_INFO, "\nachieved tol.:");
      for(int i = 0; i < _NV_; i++)
      {
        PHIST_SOUT(PHIST_INFO, "\t%8.4e", resNorm[i]);
      }
      PHIST_SOUT(PHIST_INFO, "\n");
      for(int i = 0; i < _NV_; i++)
      {
        ASSERT_LT(resNorm[i], 20*tol[i]);
      }
    }
    
    // solve _NV_ identical linear systems iteratively using MINRES or GMRES with or without a left preconditioner.
    // The first column of rhs and sol is used for all _NV_ systems, *num_iter contains max_iter on input and actual
    // num_iter on output. iflag=0 is returned if the solver converged within the given *num_iter to a residual<tol.
    void run_krylov_solver(TYPE(linearOp_ptr) A_op, TYPE(linearOp_ptr) leftPrec, TYPE(linearOp_ptr) rightPrec,
                           TYPE(mvec_ptr) rhs, TYPE(mvec_ptr) sol, 
                           phist_ElinSolv method, int *num_iter, _MT_ tol, int* iflag)
    {
      if (method==phist_GMRES||method==phist_MINRES)
      {
        // use a jadaOp for the matrix operator A without any projection space. This allows us to add
        // a left preconditioner.
        TYPE(linearOp_ptr) jdOp = new TYPE(linearOp);
        std::vector<_ST_> noSigma(_NV_,st::zero());
        SUBR(jadaOp_create)(A_op,NULL,NULL,NULL,&noSigma[0],_NV_,jdOp,&iflag_);
        ASSERT_EQ(0,iflag_);
        
        TYPE(mvec_ptr) Prhs=rhs;
        MvecOwner<_ST_> _Prhs;        
        if (leftPrec!=NULL)
        {
          SUBR(jadaOp_set_leftPrecond)(jdOp,leftPrec,&iflag_);
          ASSERT_EQ(0,iflag_);
          // left-precondition the rhs
          Prhs=NULL;
          SUBR(mvec_clone_shape)(&Prhs,rhs,&iflag_);
          _Prhs.set(Prhs);
          ASSERT_EQ(0,iflag_);
          SUBR(linearOp_apply)(st::one(),leftPrec,rhs,st::zero(),Prhs,&iflag_);
          ASSERT_EQ(0,iflag_);
        }
        
        int block_size=_NV_, max_blocks=*num_iter;
        TYPE(blockedGMRESstate_ptr) states[block_size];
        PHIST_CHK_IERR(SUBR(blockedGMRESstates_create)(states, block_size, map_, max_blocks, iflag), *iflag);

        for (int i=0; i<block_size; i++) states[i]->tol=tol;

        // copy the first rhs and sol column to all the others
        TYPE(mvec_ptr) rhs0=NULL,Prhs0=NULL,sol0=NULL;

        int nv;
        SUBR(mvec_num_vectors)(rhs,&nv,iflag);
        ASSERT_EQ(0,*iflag);
        ASSERT_EQ(_NV_,nv);
        SUBR(mvec_num_vectors)(sol,&nv,iflag);
        ASSERT_EQ(0,*iflag);
        ASSERT_EQ(_NV_,nv);

        SUBR(mvec_view_block)(Prhs,&Prhs0,0,0,iflag);
        SUBR(mvec_view_block)(rhs,&rhs0,0,0,iflag);
        ASSERT_EQ(0,*iflag);
        SUBR(mvec_view_block)(sol,&sol0,0,0,iflag);
        ASSERT_EQ(0,*iflag);
        
        // make sure these views are eventually deleted
        MvecOwner<_ST_> _rhs0(rhs0), _Prhs0(Prhs0),_sol0(sol0);
        
        for (int i=1; i<nv; i++)
        {
          SUBR(mvec_set_block)(rhs,rhs0,i,i,iflag);
          ASSERT_EQ(0,*iflag);
          SUBR(mvec_set_block)(Prhs,Prhs0,i,i,iflag);
          ASSERT_EQ(0,*iflag);
          SUBR(mvec_set_block)(sol,sol0,i,i,iflag);
          ASSERT_EQ(0,*iflag);
        }
    
        // first setup block_size systems  to work on,
    
        for (int i=0; i<block_size; i++)
        {
          // reset selected state object with 0 initial guess
          PHIST_CHK_IERR(SUBR(blockedGMRESstate_reset)(states[i], Prhs0, sol0, iflag), *iflag);
        }
        if (method==phist_GMRES)
        {
          int useIMGS=1;
          PHIST_CHK_NEG_IERR(SUBR(blockedGMRESstates_iterate)(jdOp, rightPrec,states, block_size, num_iter, useIMGS, iflag), *iflag);
        }
        else if (method==phist_MINRES)
        {
          PHIST_CHK_NEG_IERR(SUBR(blockedMINRESstates_iterate)(jdOp, rightPrec,states, block_size, num_iter, iflag), *iflag);
        }
        int min_total_iter=9999999,max_total_iter=-1;
        int num_converged=0;
        for (int i=0; i<block_size; i++)
        {
          min_total_iter=std::min(min_total_iter,states[i]->totalIter);
          max_total_iter=std::max(max_total_iter,states[i]->totalIter);
          if (states[i]->status==0)
          {
            num_converged++;
          }
        }
        // check that all of the identical linear systems converged
        ASSERT_EQ(block_size,num_converged);
        // in the same number of iterations
        ASSERT_EQ(min_total_iter,max_total_iter);
        _MT_ res_norms[block_size];
        // check that the result is correct by comparing with the actual residual
        PHIST_CHK_IERR(SUBR(blockedGMRESstates_updateSol)(states,block_size,rightPrec,sol,res_norms,false,iflag),*iflag);
     
        PHIST_CHK_IERR(SUBR(blockedGMRESstates_delete)(states,block_size,iflag),*iflag);
        PHIST_CHK_IERR(SUBR(jadaOp_delete)(jdOp,iflag),*iflag);
        delete jdOp;
      }
      else
      {
        PHIST_CHK_IERR(*iflag=PHIST_NOT_IMPLEMENTED,*iflag);
      }
    }

    TYPE(linearOp_ptr) opA_ = NULL, opP_ = NULL;
    TYPE(linearOp_ptr) jdOp_ = NULL, jdPrec_ = NULL;
    TYPE(mvec_ptr) q_=NULL, v1_=NULL,v2_=NULL,v3_=NULL,v4_=NULL;
    _ST_ *sigma_ = NULL, *negSigma_=NULL, *noSigma_=NULL;
};

  TEST_F(CLASSNAME, wrap_Ainv_as_USER_PRECON)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
        // now we have our approximate inverse as a linearOp, we can wrap it once more to get
        // a valid phist preconditioner. This object adds an 'update' function etc, which should
        // return with an error. It also takes care of checking for apply_shifted and falling back
        // to apply if it is NULL.
        TYPE(linearOp) userPrec;
        // we explicitly disable the 'apply_shifted' function of our approximate inverse
        opP_->apply_shifted=NULL;
        SUBR(precon_create)(&userPrec,ATest::A_,sigma_[0],NULL,NULL,NULL,
                "user_defined",NULL,opP_,&iflag_);
        ASSERT_EQ(0,iflag_);
        // aplying the preconditioner should be the same as applying the original Ainv matrix
        _ST_ alpha=st::prand(),beta=st::prand();
        SUBR(mvec_random)(v1_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_random)(v2_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_add_mvec)(st::one(),v2_,st::zero(),v3_,&iflag_);
        ASSERT_EQ(0,iflag_);
        SUBR(mvec_add_mvec)(st::one(),v2_,st::zero(),v4_,&iflag_);
        ASSERT_EQ(0,iflag_);
        
        // reference solution using plain spMVM interface
        SUBR(sparseMat_times_mvec)(alpha,PTest::A_,v1_,beta,v2_,&iflag_);
        ASSERT_EQ(0,iflag_);

        // sanity check if our linearOp wrapper works
        SUBR(linearOp_apply)(alpha,opP_,v1_,beta,v3_,&iflag_);
        ASSERT_EQ(0,iflag_);
        ASSERT_NEAR(mt::one(),MvecsEqual(v2_,v3_),10*VTest::releps());

        // yet another wrapper ontop: userPrec.A_ is opP_ wrapped as an internal precon struct with a phist_Eprecon to identify it as USER_DEFINED
        SUBR(linearOp_apply)(alpha,&userPrec,v1_,beta,v4_,&iflag_);
        ASSERT_EQ(0,iflag_);
        ASSERT_NEAR(mt::one(),MvecsEqual(v2_,v4_),10*VTest::releps());
    
        SUBR(precon_delete)(&userPrec,&iflag_);
        ASSERT_EQ(0,iflag_);
    }
  }

  // test if q is orthonormal (prerequisite for other tests, although it does not have anything to do with preconditioning)
  TEST_F(CLASSNAME, q_is_orthonormal)
  {
    _ST_* q_vp=NULL;
    phist_lidx ldq;
    SUBR(mvec_extract_view)(q_,&q_vp,&ldq,&iflag_);
    ASSERT_EQ(0,iflag_);
    EXPECT_NEAR(mt::one(),QTest::ColsAreNormalized(q_vp,nloc_,ldq,QTest::stride_,mpi_comm_),(MT)100.*QTest::releps(q_));
    EXPECT_NEAR(mt::one(),QTest::ColsAreOrthogonal(q_vp,nloc_,ldq,QTest::stride_,mpi_comm_),(MT)100.*QTest::releps(q_));
  }

  // apply wrapped projected preconditioner and check PJD\Q = (I-(P\Q)inv(Q'P\Q)Q')P\Q = 0
  TEST_F(CLASSNAME, apply_jadaPrec)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(mvec_ptr) pq=NULL;
      SUBR(mvec_clone_shape)(&pq,q_,&iflag_);
      ASSERT_EQ(0,iflag_);
      MvecOwner<_ST_> _pq(pq);
      SUBR(linearOp_apply)(st::one(),jdPrec_,q_,st::zero(),pq,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      // should give 0 if the projection works
      ASSERT_NEAR(mt::one(),MvecEqual(pq,st::zero()),VTest::releps());
    }
  }


  TEST_F(CLASSNAME, precon_improves_krylov_solver)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      SUBR(mvec_put_value)(v1_,st::one(),&iflag_);
      ASSERT_EQ(0,iflag_);
      _MT_ nrm0[_NV_];
      SUBR(mvec_normalize)(v1_,nrm0,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_put_value)(v2_,st::zero(),&iflag_);
      ASSERT_EQ(0,iflag_);

      // compute initial residual norm. vec1 is b, vec2 x0, vec3 will be b-Ax0
      SUBR(mvec_add_mvec)(st::one(),v1_,st::zero(),v3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sparseMat_times_mvec)(-st::one(),ATest::A_,v2_,st::one(),v3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      _MT_ normR0[_NV_];
      SUBR(mvec_norm2)(v3_,normR0,&iflag_);
            ASSERT_EQ(0,iflag_);

      // copy X0 to vec3 and vec4
      SUBR(mvec_add_mvec)(st::one(),v2_,st::zero(),v3_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(mvec_add_mvec)(st::one(),v2_,st::zero(),v4_,&iflag_);
      ASSERT_EQ(0,iflag_);
      
      int max_iter=_MAXBAS_;
      _MT_ tol=0.1;
      int num_iter0=max_iter, num_iter1=max_iter, num_iter2=max_iter;
      
      // run GMRES or MINRES without preconditioning
      run_krylov_solver(opA_, NULL, NULL,
                        v1_, v2_,
                        _SOLVTYPE_, &num_iter0, tol, &iflag_);
      ASSERT_EQ(0,iflag_);

      // with left preconditioning
      run_krylov_solver(opA_, opP_, NULL,
                        v1_, v3_, 
                        _SOLVTYPE_, &num_iter1, tol, &iflag_);
      ASSERT_EQ(0,iflag_);
      
      // with right preconditioning
      run_krylov_solver(opA_, NULL, opP_,
                        v1_, v4_, 
                        _SOLVTYPE_, &num_iter2, tol, &iflag_);
      ASSERT_EQ(0,iflag_);
    
      // using a preconditioner should reduce the number of iteations somewhat
      ASSERT_LT(num_iter1,(int)(0.8*num_iter0));
      // left or right preconditioning should do roughly the same
      ASSERT_LT(std::abs(num_iter1-num_iter2),2);
      // and the results should be almost the same, too
//      VTest::PrintVector(PHIST_VERBOSE, "x without preconditioning", 
//                VTest::vec2_vp_, VTest::nloc_, VTest::lda_, VTest::stride_,mpi_comm_);
//      VTest::PrintVector(PHIST_VERBOSE, "x with left preconditioning", 
//                VTest::vec3_vp_, nloc_, VTest::lda_, VTest::stride_,mpi_comm_);
//      VTest::PrintVector(PHIST_VERBOSE, "x with right preconditioning", 
//                VTest::vec4_vp_, nloc_, VTest::lda_, VTest::stride_,mpi_comm_);
                
      EXPECT_NEAR(mt::one(),MvecsEqual(v2_,v3_),tol*normR0[0]);
      EXPECT_NEAR(mt::one(),MvecsEqual(v2_,v4_),tol*normR0[0]);
      EXPECT_NEAR(mt::one(),MvecsEqual(v3_,v4_),tol*normR0[0]);
    }
  }


  TEST_F(CLASSNAME, DISABLED_jadaCorrectionSolver)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      TYPE(jadaCorrectionSolver_ptr) solver = NULL;
      SUBR(jadaCorrectionSolver_create)(&solver, jadaOpts_, map_, &iflag_);
      ASSERT_EQ(0, iflag_);

      _MT_ tol[_NV_];
      for(int i = 0; i < _NV_; i++)
      {
        // create some random tolerance
        tol[i] = (_MT_)exp((_MT_)-8*mt::one() + (_MT_)4*mt::prand());
      }

      SUBR(mvec_put_value)(v2_, st::zero(), &iflag_);
      ASSERT_EQ(0, iflag_);

      // run
      SUBR(jadaCorrectionSolver_run)(solver, opA_, NULL, q_, NULL, sigma_,v3_, NULL, tol, 200, v2_, 1, 0, 0, &iflag_);
      ASSERT_EQ(0, iflag_);

      // check all solutions
      checkResiduals(tol);

      SUBR(jadaCorrectionSolver_delete)(solver, &iflag_);
      ASSERT_EQ(0, iflag_);
    }
  }

  

  TEST_F(CLASSNAME, DISABLED_subspacejada)
  {
    if( typeImplemented_ && !problemTooSmall_ )
    {
      int nEig1=3,nEig2=nEig1;
      jadaOpts_.numEigs=nEig1;
      jadaOpts_.which=_WHICH_;
      int blockDim=_NV_;
      int nIter1=100,nIter2=2*nIter1;
      int nQ=nEig1+blockDim-1;
      TYPE(mvec_ptr) Q1=NULL,Q2=NULL;
      PHISTTEST_MVEC_CREATE(&Q1,map_,nQ,&iflag_);
      ASSERT_EQ(0,iflag_);
      PHISTTEST_MVEC_CREATE(&Q2,map_,nQ,&iflag_);
      ASSERT_EQ(0,iflag_);
      TYPE(sdMat_ptr) R1=NULL,R2=NULL;
      SUBR(sdMat_create)(&R1,nQ,nQ,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      SUBR(sdMat_create)(&R2,nQ,nQ,comm_,&iflag_);
      ASSERT_EQ(0,iflag_);
      MvecOwner<_ST_> _Q1(Q1),_Q2(Q2);
      SdMatOwner<_ST_> _R1(R1),_R2(R2);
      
      _CT_ ev1[nQ],ev2[nQ];
      _MT_ resNorm1[nQ],resNorm2[nQ];
      
      // set starting vector, note that this test is only compled with _NVP_==1
      jadaOpts_.v0=q_;
      
      // use fixed number of inner iterations because preconditioning distorts the GMRES
      // convergence tolerance
      jadaOpts_.innerSolvBaseTol=0.0;
      
      // run with preconditioner
      SUBR(subspacejada)(opA_, NULL, jadaOpts_,
          Q1, R1, ev1, resNorm1, &nEig1, &nIter1, &iflag_);
      EXPECT_EQ(0,iflag_);

      // run without preconditioner
      jadaOpts_.preconOp=NULL;
      SUBR(subspacejada)(opA_, NULL, jadaOpts_,
          Q2, R2, ev2, resNorm2, &nEig2, &nIter2, &iflag_);
      EXPECT_EQ(0,iflag_);
      
      // the preconditioner should at least not hurt the overall convergence:
      ASSERT_LE(nIter1,nIter2);
    }
  }
