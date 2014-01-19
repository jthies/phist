#ifndef CLASSNAME
#error "file not included correctly"
#endif

/*! Test fixure. */
class CLASSNAME: public KernelTestWithSdMats<_ST_,_N_,_N_>
{

  public:

    typedef KernelTestWithSdMats<_ST_,_N_,_N_> MTest;

    //! mvec/sdMat sizes
    static const int n_=_N_;
    static const int m_=_N_;

    // different test scenarios
    std::vector<int> nselect_;
    std::vector<int> nsort_;

    std::complex<MT> ev_[_N_];

    // for testing PartialSchurDecompReorder
    std::vector<int> permutation;
    std::vector<MT> resNorm;

    /*! Set up routine.
    */
    virtual void SetUp()
    {
      MTest::SetUp();

      if( typeImplemented_ )
      {
        // a few test cases
        nselect_.push_back(0);                  nsort_.push_back(0);
        nselect_.push_back(1);                  nsort_.push_back(0);
        nselect_.push_back(1);                  nsort_.push_back(1);
        nselect_.push_back(_N_);                nsort_.push_back(_N_);
        nselect_.push_back(std::min(5,_N_));    nsort_.push_back(0);
        nselect_.push_back(std::min(7,_N_));    nsort_.push_back(3);
        nselect_.push_back(std::min(20,_N_));   nsort_.push_back(7);
        nselect_.push_back(std::min(8,_N_));    nsort_.push_back(8);

        // create a diagonal matrix with some interesting features for the diag_* tests
        SUBR(sdMat_put_value)(mat3_,st::zero(),&ierr_);
        ASSERT_EQ(0,ierr_);
        ST *diag = new ST[nrows_];
        for (int i=0;i<nrows_;i++)
        {
          diag[i]=st::prand();
        }
        if (nrows_==10)
        {
          diag[6]=diag[1];
          diag[9]=diag[1];
          diag[0]=(ST)st::real(diag[0]);
          diag[4]=-diag[0];
        }
        if (nrows_==50)
        {
          for(int i = 0; i < nrows_; i++)
          {
            diag[i] = i*(i+1) % 20;
          }
        }
        for (int i=0;i<nrows_;i++)
        {
          mat3_vp_[i*m_lda_+i]=diag[i];
        }
        delete [] diag;
      }
    }

    /*! Clean up.
    */
    virtual void TearDown()
    {
      KernelTestWithType<_ST_>::TearDown();
    }

    void DoSchurDecompTest(TYPE(const_sdMat_ptr) A_clone, eigSort_t which, _MT_ tol, bool onlyDoReorderTest)
    {
      if (!typeImplemented_) return;
      ASSERT_EQ(nsort_.size(),nselect_.size());
      for (int c=0; c < nselect_.size(); c++)
      {
        int nselect = nselect_[c];
        int nsort = nsort_[c];
        if( onlyDoReorderTest && nsort == 0 )
          continue;
        PHIST_OUT(PHIST_INFO,"==================================================\n");
        if( which == LM ) {
          PHIST_OUT(PHIST_INFO,"CASE LM nselect %d, nsort %d, tol %e\n",nselect,nsort, tol);
        } else if( which == SM ) {
          PHIST_OUT(PHIST_INFO,"CASE SM nselect %d, nsort %d, tol %e\n",nselect,nsort, tol);
        } else if( which == LR ) {
          PHIST_OUT(PHIST_INFO,"CASE LR nselect %d, nsort %d, tol %e\n",nselect,nsort, tol);
        } else if( which == SR ) {
          PHIST_OUT(PHIST_INFO,"CASE SR nselect %d, nsort %d, tol %e\n",nselect,nsort, tol);
        }
        PHIST_OUT(PHIST_INFO,"==================================================\n");

        SUBR(sdMat_add_sdMat)(st::one(),A_clone,st::zero(),mat1_,&this->ierr_);
        ASSERT_EQ(0,this->ierr_);

        SUBR(sdMat_random)(mat2_,&this->ierr_);
        ASSERT_EQ(0,this->ierr_);

#if PHIST_OUTLEV>=PHIST_DEBUG
        PHIST_DEB("input matrix to Schur-decomp:\n");
        SUBR(sdMat_print)(mat1_,&ierr_);
#endif
        SUBR(SchurDecomp)(mat1_vp_,m_lda_,mat2_vp_,m_lda_,n_,nselect,nsort,which,tol,ev_,&this->ierr_);
        PHIST_DEB("resulting T:\n");
#if PHIST_OUTLEV>=PHIST_DEBUG
        SUBR(sdMat_print)(mat1_,&ierr_);
        ASSERT_EQ(0,ierr_);
#endif


        PHIST_OUT(PHIST_INFO,"Checking results...\n");
        CheckSchurDecomp(A_clone, which, nselect, nsort, tol);
        if( HasFatalFailure() )
          return;

        // don't do this always so the results don't get too messed up
        if( onlyDoReorderTest )
        {
          std::vector<_MT_> resNorm(nsort);
          std::vector<_MT_> resNormOrig(nsort);
          std::vector<int> permutation(nsort);
          for(int i = 0; i < nsort; i++)
            resNormOrig[i] = resNorm[i] = exp(mt::prand());
          SUBR(ReorderPartialSchurDecomp)(mat1_vp_,m_lda_,mat2_vp_,m_lda_,n_,nsort,which,tol,
              &resNorm[0],ev_,&permutation[0],&ierr_);
          ASSERT_EQ(0,ierr_);
          PHIST_DEB("resulting T:\n");
#if PHIST_OUTLEV>=PHIST_DEBUG
          SUBR(sdMat_print)(mat1_,&ierr_);
          ASSERT_EQ(0,ierr_);
#endif


          // should still be a valid schur decomposition!
          PHIST_OUT(PHIST_INFO,"Checking reordering results...\n");
          CheckSchurDecomp(A_clone, which, nsort, nsort, tol);

          // check if the permutation of resNorm is correct
          PHIST_OUT(PHIST_INFO,"resNorm array:\nold\t\tnew\t\tperm\n");
          for(int i = 0; i < nsort; i++)
          {
            PHIST_OUT(PHIST_INFO,"%8.4e\t%8.4e\t%d\n",resNormOrig[i],resNorm[i],permutation[i]);
            ASSERT_TRUE( resNorm.at(i) == resNormOrig.at(permutation[i]) );
          }

          // test if the reordering itself was successful
          for(int i = 1; i < nsort; i++)
          {
            if (which==LM || which==SM)
            {
              if( mt::abs(ct::abs(ev_[i])-ct::abs(ev_[i-1])) > tol )
                continue;
            }
            else if( which==LR || which==SR )
            {
              if( mt::abs(ct::real(ev_[i])-ct::real(ev_[i-1])) > tol )
                continue;
            }

            ASSERT_TRUE( resNorm[i] > resNorm[i-1] );
          }
        }
      }//cases
    }//DoSchurDecompTest


    void CheckSchurDecomp(TYPE(const_sdMat_ptr) A_clone, eigSort_t which, int nselect, int nsort, _MT_ tol)
    {
      PHIST_DEB("check AS=ST\n");
      SUBR(sdMat_times_sdMat)(st::one(),A_clone,mat2_,st::zero(),mat4_,&ierr_);
      ASSERT_EQ(0,ierr_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      //SUBR(sdMat_print)(mat4_,&ierr_);
      //ASSERT_EQ(0,ierr_);
#endif
      SUBR(sdMat_times_sdMat)(-st::one(),mat2_,mat1_,st::one(),mat4_,&ierr_);
      ASSERT_EQ(0,ierr_);
#if PHIST_OUTLEV>=PHIST_DEBUG
      //SUBR(sdMat_print)(mat4_,&ierr_);
      //ASSERT_EQ(0,ierr_);
#endif
      ASSERT_NEAR(mt::one(),ArrayEqual(mat4_vp_,nrows_,ncols_,m_lda_,1,st::zero()),1000*mt::eps());

      PHIST_OUT(PHIST_INFO,"eigenvalue array:\n");
      for (int i=0;i<n_;i++)
      {
        // test the traits class on the way:
        ASSERT_REAL_EQ(ct::abs(ev_[i]),std::abs(ev_[i]));
        ASSERT_REAL_EQ(ct::real(ev_[i]),std::real(ev_[i]));
        ASSERT_REAL_EQ(ct::imag(ev_[i]),std::imag(ev_[i]));
        PHIST_OUT(PHIST_INFO,"%8.4f%+8.4fi\tabs=%8.4f\n",ct::real(ev_[i]),ct::imag(ev_[i]),ct::abs(ev_[i]));
      }


      // check that the eigenvalues on the diagonal of T have the same ordering as those in 
      // ev_
      for (int i=0;i<n_;i++)
      {
        ASSERT_REAL_EQ(ct::real(ev_[i]), st::real(mat1_vp_[i*m_lda_+i]));
#ifdef IS_COMPLEX
        ASSERT_REAL_EQ(ct::imag(ev_[i]), st::imag(mat1_vp_[i*m_lda_+i]));
#endif      
      }
      // check that the first nsort eigenvalues are sorted correctly
      for (int i=1;i<nsort;i++)
      {
        if (which==LM)
        {
          ASSERT_LE(ct::abs(ev_[i]), ct::abs(ev_[i-1])+tol);
        }
        else if (which==SM)
        {
          ASSERT_GE(ct::abs(ev_[i]), ct::abs(ev_[i-1])-tol);
        }      
        else if (which==LR)
        {
          ASSERT_LE(ct::real(ev_[i]), ct::real(ev_[i-1])+tol);
        }
        else if (which==SR)
        {
          ASSERT_GE(ct::real(ev_[i]), ct::real(ev_[i-1])-tol);
        }
      }
      // check that the first nselect are the largest/smallest etc globally
      if( nselect > 0 )
      {
        MT val;
        if (which==LM||which==SM) val=ct::abs(ev_[0]);
        else if (which==LR||which==SR) val=ct::real(ev_[0]);
        for (int i=1;i<nselect;i++)
        {
          if (which==LM)
          {
            val=std::min(ct::abs(ev_[i]),val);
          }
          else if (which==SM)
          {
            val=std::max(ct::abs(ev_[i]),val);
          }      
          else if (which==LR)
          {
            val=std::min(ct::real(ev_[i]),val);
          }
          else if (which==SR)
          {
            val=std::max(ct::real(ev_[i]),val);
          }
        }
        for (int i=nselect;i<n_;i++)
        {
          if (which==LM)
          {
            ASSERT_LE(ct::abs(ev_[i]), val+tol);
          }
          else if (which==SM)
          {
            ASSERT_GE(ct::abs(ev_[i]), val-tol);
          }
          if (which==LR)
          {
            ASSERT_LE(ct::real(ev_[i]), val+tol);
          }
          else if (which==SR)
          {
            ASSERT_GE(ct::real(ev_[i]), val-tol);
          }
        }
      }

      // check the T matrix is upper triangular (with 2x2 blocks for complex eigs in the 
      // real case)
#ifdef IS_COMPLEX
      for (int j=0;j<_N_;j++)
      {
        for (int i=j+1;i<_N_;i++)
        {
          ASSERT_REAL_EQ(mt::one(),mt::one()+st::abs(mat1_vp_[j*m_lda_+i]));
        }//i
      }//j
#else
      for (int j=0;j<_N_;j++)
      {
        PHIST_DEB("j=%d, ev[j]=%8.4f %+8.4fi\n",j,ct::real(ev_[j]),ct::imag(ev_[j]));
        if (ct::imag(ev_[j])<=mt::zero() && j+1<_N_)
        {
          ASSERT_REAL_EQ(mt::one(),mt::one()+st::abs(mat1_vp_[j*m_lda_+j+1]));
        }
        for (int i=j+2;i<_N_;i++)
        {
          ASSERT_REAL_EQ(mt::one(),mt::one()+st::abs(mat1_vp_[j*m_lda_+i]));
        }//i
      }//j
#endif
    }//DoSchurDecompTest
};

TEST_F(CLASSNAME, rand_LM) 
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->ierr_);
    DoSchurDecompTest(mat3_,LM,mt::zero(),false);
  }
}

TEST_F(CLASSNAME, rand_SM) 
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->ierr_);
    DoSchurDecompTest(mat3_,SM,mt::zero(),false);
  }
}

TEST_F(CLASSNAME, rand_LR) 
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->ierr_);
    DoSchurDecompTest(mat3_,LR,mt::zero(),false);
  }
}

TEST_F(CLASSNAME, rand_SR) 
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->ierr_);
    DoSchurDecompTest(mat3_,SR,mt::zero(),false);
  }
}

TEST_F(CLASSNAME, diag_LM) 
{
  DoSchurDecompTest(mat3_,LM,mt::zero(),false);
}

TEST_F(CLASSNAME, diag_SM) 
{
  DoSchurDecompTest(mat3_,SM,mt::zero(),false);
}

TEST_F(CLASSNAME, diag_LR) 
{
  DoSchurDecompTest(mat3_,LR,mt::zero(),false);
}

TEST_F(CLASSNAME, diag_SR) 
{
  DoSchurDecompTest(mat3_,SR,mt::zero(),false);
}


TEST_F(CLASSNAME, rand_LM_tol)
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->ierr_);
    DoSchurDecompTest(mat3_,LM,(_MT_)0.3,false);
  }
}

TEST_F(CLASSNAME, rand_SM_tol)
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->ierr_);
    DoSchurDecompTest(mat3_,SM,(_MT_)0.7,false);
  }
}

TEST_F(CLASSNAME, rand_LR_tol)
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->ierr_);
    DoSchurDecompTest(mat3_,LR,(_MT_)0.2,false);
  }
}

TEST_F(CLASSNAME, rand_SR_tol)
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->ierr_);
    DoSchurDecompTest(mat3_,SR,(_MT_)0.5,false);
  }
}

TEST_F(CLASSNAME, diag_LM_tol)
{
  DoSchurDecompTest(mat3_,LM,(_MT_)0.3,false);
}

TEST_F(CLASSNAME, diag_SM_tol)
{
  DoSchurDecompTest(mat3_,SM,(_MT_)0.2,false);
}

TEST_F(CLASSNAME, diag_LR_tol)
{
  DoSchurDecompTest(mat3_,LR,(_MT_)0.1,false);
}

TEST_F(CLASSNAME, diag_SR_tol)
{
  DoSchurDecompTest(mat3_,SR,(_MT_)0.03,false);
}


#ifdef IS_COMPLEX
TEST_F(CLASSNAME, rand_LM_tol_reorder)
#else
TEST_F(CLASSNAME, DISABLED_rand_LM_tol_reorder)
#endif
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->ierr_);
    DoSchurDecompTest(mat3_,LM,(_MT_)0.3,true);
  }
}

#ifdef IS_COMPLEX
TEST_F(CLASSNAME, rand_SM_tol_reorder)
#else
TEST_F(CLASSNAME, DISABLED_rand_SM_tol_reorder)
#endif
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->ierr_);
    DoSchurDecompTest(mat3_,SM,(_MT_)0.7,true);
  }
}

#ifdef IS_COMPLEX
TEST_F(CLASSNAME, rand_LR_tol_reorder)
#else
TEST_F(CLASSNAME, DISABLED_rand_LR_tol_reorder)
#endif
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->ierr_);
    DoSchurDecompTest(mat3_,LR,(_MT_)0.2,true);
  }
}

#ifdef IS_COMPLEX
TEST_F(CLASSNAME, rand_SR_tol_reorder)
#else
TEST_F(CLASSNAME, DISABLED_rand_SR_tol_reorder)
#endif
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->ierr_);
    DoSchurDecompTest(mat3_,SR,(_MT_)0.5, true);
  }
}

TEST_F(CLASSNAME, diag_LM_tol_reorder)
{
  DoSchurDecompTest(mat3_,LM,(_MT_)0.3,true);
}

TEST_F(CLASSNAME, diag_SM_tol_reorder)
{
  DoSchurDecompTest(mat3_,SM,(_MT_)0.2,true);
}

TEST_F(CLASSNAME, diag_LR_tol_reorder)
{
  DoSchurDecompTest(mat3_,LR,(_MT_)0.1,true);
}

TEST_F(CLASSNAME, diag_SR_tol_reorder)
{
  DoSchurDecompTest(mat3_,SR,(_MT_)0.03,true);
}

