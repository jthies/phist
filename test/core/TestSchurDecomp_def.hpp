/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
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
        SUBR(sdMat_put_value)(mat3_,st::zero(),&iflag_);
        PHIST_CHK_IERR(SUBR(sdMat_from_device)(mat3_,&iflag_),iflag_);
        ASSERT_EQ(0,iflag_);
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
        PHIST_CHK_IERR(SUBR(sdMat_to_device)(mat3_,&iflag_),iflag_);
        delete [] diag;
      }
    }

    /*! Clean up.
    */
    virtual void TearDown()
    {
      MTest::TearDown();
    }

    void DoSchurDecompTest(phist_EeigSort which, _MT_ tol, bool onlyDoReorderTest)
    {
      if( !typeImplemented_ )
        return;
      ASSERT_EQ(nsort_.size(),nselect_.size());
      for (int c=0; c < nselect_.size(); c++)
      {
        int nselect = nselect_[c];
        int nsort = nsort_[c];
        if( onlyDoReorderTest && nsort == 0 )
          continue;
        PHIST_SOUT(PHIST_VERBOSE,"==================================================\n");
        if( which == phist_LM ) {
          PHIST_SOUT(PHIST_VERBOSE,"CASE phist_LM nselect %d, nsort %d, tol %e\n",nselect,nsort, tol);
        } else if( which == phist_SM ) {
          PHIST_SOUT(PHIST_VERBOSE,"CASE phist_SM nselect %d, nsort %d, tol %e\n",nselect,nsort, tol);
        } else if( which == phist_LR ) {
          PHIST_SOUT(PHIST_VERBOSE,"CASE phist_LR nselect %d, nsort %d, tol %e\n",nselect,nsort, tol);
        } else if( which == phist_SR ) {
          PHIST_SOUT(PHIST_VERBOSE,"CASE phist_SR nselect %d, nsort %d, tol %e\n",nselect,nsort, tol);
        }
        PHIST_SOUT(PHIST_VERBOSE,"==================================================\n");

        SUBR(sdMat_add_sdMat)(st::one(),mat3_,st::zero(),mat1_,&this->iflag_);
        ASSERT_EQ(0,this->iflag_);

        SUBR(sdMat_random)(mat2_,&this->iflag_);
        ASSERT_EQ(0,this->iflag_);

        SUBR(sdMat_sync_values)(mat2_, comm_, &iflag_);
        ASSERT_EQ(0,iflag_);
      
        PrintSdMat(PHIST_DEBUG,"input matrix to Schur-decomp",
                mat1_vp_, m_lda_, 1,mpi_comm_);
        
        PHIST_CHK_IERR(SUBR(sdMat_from_device)(mat1_,&iflag_),iflag_);
        SUBR(SchurDecomp)(mat1_vp_,m_lda_,mat2_vp_,m_lda_,n_,nselect,nsort,which,tol,ev_,&this->iflag_);
        PHIST_CHK_IERR(SUBR(sdMat_to_device)(mat1_,&iflag_),iflag_);
        PHIST_CHK_IERR(SUBR(sdMat_to_device)(mat2_,&iflag_),iflag_);

        PrintSdMat(PHIST_DEBUG,"resulting T",
                mat1_vp_, m_lda_, 1,mpi_comm_);


        PHIST_SOUT(PHIST_VERBOSE,"Checking results...\n");
        CheckSchurDecomp(which, nselect, nsort, tol);
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
              &resNorm[0],ev_,&permutation[0],&iflag_);
          ASSERT_EQ(0,iflag_);
          SUBR(sdMat_to_device)(mat1_,&iflag_);
          ASSERT_EQ(0,iflag_);
          SUBR(sdMat_to_device)(mat2_,&iflag_);
          ASSERT_EQ(0,iflag_);
          
        PrintSdMat(PHIST_DEBUG,"resulting T",
                mat1_vp_, m_lda_, 1,mpi_comm_);


          // should still be a valid schur decomposition!
          PHIST_SOUT(PHIST_VERBOSE,"Checking reordering results...\n");
          CheckSchurDecomp(which, nsort, nsort, tol);

          // check if the permutation of resNorm is correct
          PHIST_SOUT(PHIST_DEBUG,"resNorm array:\nold\t\tnew\t\tperm\n");
          for(int i = 0; i < nsort; i++)
          {
            PHIST_SOUT(PHIST_DEBUG,"%8.4e\t%8.4e\t%d\n",resNormOrig[i],resNorm[i],permutation[i]);
            ASSERT_TRUE( resNorm.at(i) == resNormOrig.at(permutation[i]) );
          }

          // test if the reordering itself was successful
          for(int i = 1; i < nsort; i++)
          {
            if (which==phist_LM || which==phist_SM)
            {
              if( mt::abs(ct::abs(ev_[i])-ct::abs(ev_[i-1])) > tol )
                continue;
            }
            else if( which==phist_LR || which==phist_SR )
            {
              if( mt::abs(ct::real(ev_[i])-ct::real(ev_[i-1])) > tol )
                continue;
            }

            ASSERT_TRUE( resNorm[i] > resNorm[i-1] );
          }
        }
      }//cases
    }//DoSchurDecompTest


    void CheckSchurDecomp(phist_EeigSort which, int nselect, int nsort, _MT_ tol)
    {
      PHIST_DEB("check AS=ST\n");
      SUBR(sdMat_times_sdMat)(st::one(),mat3_,mat2_,st::zero(),mat4_,&iflag_);
      ASSERT_EQ(0,iflag_);

      SUBR(sdMat_times_sdMat)(-st::one(),mat2_,mat1_,st::one(),mat4_,&iflag_);
      ASSERT_EQ(0,iflag_);

      PHIST_CHK_IERR(SUBR(sdMat_from_device)(mat4_,&iflag_),iflag_);
      ASSERT_NEAR(mt::one(),ArrayEqual(mat4_vp_,nrows_,ncols_,m_lda_,1,st::zero()),1000*mt::eps());

      PHIST_SOUT(PHIST_DEBUG,"eigenvalue array:\n");
      for (int i=0;i<n_;i++)
      {
        // test the traits class on the way:
        ASSERT_REAL_EQ(ct::abs(ev_[i]),std::abs(ev_[i]));
        ASSERT_REAL_EQ(ct::real(ev_[i]),std::real(ev_[i]));
        ASSERT_REAL_EQ(ct::imag(ev_[i]),std::imag(ev_[i]));
        PHIST_SOUT(PHIST_DEBUG,"%8.4f%+8.4fi\tabs=%8.4f\n",ct::real(ev_[i]),ct::imag(ev_[i]),ct::abs(ev_[i]));
      }


      // check that the eigenvalues on the diagonal of T have the same ordering as those in 
      // ev_
      PHIST_CHK_IERR(SUBR(sdMat_from_device)(mat1_,&iflag_),iflag_);
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
        if (which==phist_LM)
        {
          ASSERT_LE(ct::abs(ev_[i]), ct::abs(ev_[i-1])+tol);
        }
        else if (which==phist_SM)
        {
          ASSERT_GE(ct::abs(ev_[i]), ct::abs(ev_[i-1])-tol);
        }      
        else if (which==phist_LR)
        {
          ASSERT_LE(ct::real(ev_[i]), ct::real(ev_[i-1])+tol);
        }
        else if (which==phist_SR)
        {
          ASSERT_GE(ct::real(ev_[i]), ct::real(ev_[i-1])-tol);
        }
      }
      // check that the first nselect are the largest/smallest etc globally
      if( nselect > 0 )
      {
        MT val=mt::zero();
        if (which==phist_LM||which==phist_SM) val=ct::abs(ev_[0]);
        else if (which==phist_LR||which==phist_SR) val=ct::real(ev_[0]);
        for (int i=1;i<nselect;i++)
        {
          if (which==phist_LM)
          {
            val=std::min(ct::abs(ev_[i]),val);
          }
          else if (which==phist_SM)
          {
            val=std::max(ct::abs(ev_[i]),val);
          }      
          else if (which==phist_LR)
          {
            val=std::min(ct::real(ev_[i]),val);
          }
          else if (which==phist_SR)
          {
            val=std::max(ct::real(ev_[i]),val);
          }
        }
        for (int i=nselect;i<n_;i++)
        {
          if (which==phist_LM)
          {
            ASSERT_LE(ct::abs(ev_[i]), val+tol);
          }
          else if (which==phist_SM)
          {
            ASSERT_GE(ct::abs(ev_[i]), val-tol);
          }
          if (which==phist_LR)
          {
            ASSERT_LE(ct::real(ev_[i]), val+tol);
          }
          else if (which==phist_SR)
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
        else
        {
          ASSERT_REAL_EQ(mt::one(),mt::one());
        }
        for (int i=j+2;i<_N_;i++)
        {
          ASSERT_REAL_EQ(mt::one(),mt::one()+st::abs(mat1_vp_[j*m_lda_+i]));
        }//i
      }//j
#endif

      // check that we have the same result on all processes
      int stride = 1;
      ASSERT_REAL_EQ(mt::one(), ArrayParallelReplicated(mat3_vp_,nrows_,ncols_,m_lda_,stride,mflag_)); // <- input matrix
      ASSERT_REAL_EQ(mt::one(), ArrayParallelReplicated(mat1_vp_,nrows_,ncols_,m_lda_,stride,mflag_));
      ASSERT_REAL_EQ(mt::one(), ArrayParallelReplicated(mat2_vp_,nrows_,ncols_,m_lda_,stride,mflag_));
    }//CheckSchurDecomp
};

TEST_F(CLASSNAME, rand_LM) 
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->iflag_);
    SUBR(sdMat_sync_values)(mat3_, comm_, &iflag_);
    DoSchurDecompTest(phist_LM,mt::zero(),false);
  }
}

TEST_F(CLASSNAME, rand_SM) 
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->iflag_);
    SUBR(sdMat_sync_values)(mat3_, comm_, &iflag_);
    DoSchurDecompTest(phist_SM,mt::zero(),false);
  }
}

TEST_F(CLASSNAME, rand_LR) 
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->iflag_);
    SUBR(sdMat_sync_values)(mat3_, comm_, &iflag_);
    DoSchurDecompTest(phist_LR,mt::zero(),false);
  }
}

TEST_F(CLASSNAME, rand_SR) 
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->iflag_);
    SUBR(sdMat_sync_values)(mat3_, comm_, &iflag_);
    DoSchurDecompTest(phist_SR,mt::zero(),false);
  }
}

TEST_F(CLASSNAME, diag_LM) 
{
  DoSchurDecompTest(phist_LM,mt::zero(),false);
}

TEST_F(CLASSNAME, diag_SM) 
{
  DoSchurDecompTest(phist_SM,mt::zero(),false);
}

TEST_F(CLASSNAME, diag_LR) 
{
  DoSchurDecompTest(phist_LR,mt::zero(),false);
}

TEST_F(CLASSNAME, diag_SR) 
{
  DoSchurDecompTest(phist_SR,mt::zero(),false);
}


TEST_F(CLASSNAME, rand_LM_tol)
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->iflag_);
    SUBR(sdMat_sync_values)(mat3_, comm_, &iflag_);
    DoSchurDecompTest(phist_LM,(_MT_)0.3,false);
  }
}

TEST_F(CLASSNAME, rand_SM_tol)
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->iflag_);
    SUBR(sdMat_sync_values)(mat3_, comm_, &iflag_);
    DoSchurDecompTest(phist_SM,(_MT_)0.7,false);
  }
}

TEST_F(CLASSNAME, rand_LR_tol)
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->iflag_);
    SUBR(sdMat_sync_values)(mat3_, comm_, &iflag_);
    DoSchurDecompTest(phist_LR,(_MT_)0.2,false);
  }
}

TEST_F(CLASSNAME, rand_SR_tol)
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->iflag_);
    SUBR(sdMat_sync_values)(mat3_, comm_, &iflag_);
    DoSchurDecompTest(phist_SR,(_MT_)0.5,false);
  }
}

TEST_F(CLASSNAME, diag_LM_tol)
{
  DoSchurDecompTest(phist_LM,(_MT_)0.3,false);
}

TEST_F(CLASSNAME, diag_SM_tol)
{
  DoSchurDecompTest(phist_SM,(_MT_)0.2,false);
}

TEST_F(CLASSNAME, diag_LR_tol)
{
  DoSchurDecompTest(phist_LR,(_MT_)0.1,false);
}

TEST_F(CLASSNAME, diag_SR_tol)
{
  DoSchurDecompTest(phist_SR,(_MT_)0.03,false);
}


#ifdef IS_COMPLEX
TEST_F(CLASSNAME, rand_LM_tol_reorder)
#else
TEST_F(CLASSNAME, DISABLED_rand_LM_tol_reorder)
#endif
{
  if( typeImplemented_ )
  {
    SUBR(sdMat_random)(mat3_,&this->iflag_);
    SUBR(sdMat_sync_values)(mat3_, comm_, &iflag_);
    DoSchurDecompTest(phist_LM,(_MT_)0.3,true);
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
    SUBR(sdMat_random)(mat3_,&this->iflag_);
    SUBR(sdMat_sync_values)(mat3_, comm_, &iflag_);
    DoSchurDecompTest(phist_SM,(_MT_)0.7,true);
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
    SUBR(sdMat_random)(mat3_,&this->iflag_);
    SUBR(sdMat_sync_values)(mat3_, comm_, &iflag_);
    DoSchurDecompTest(phist_LR,(_MT_)0.2,true);
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
    SUBR(sdMat_random)(mat3_,&this->iflag_);
    SUBR(sdMat_sync_values)(mat3_, comm_, &iflag_);
    DoSchurDecompTest(phist_SR,(_MT_)0.5, true);
  }
}

TEST_F(CLASSNAME, diag_LM_tol_reorder)
{
  DoSchurDecompTest(phist_LM,(_MT_)0.3,true);
}

TEST_F(CLASSNAME, diag_SM_tol_reorder)
{
  DoSchurDecompTest(phist_SM,(_MT_)0.2,true);
}

TEST_F(CLASSNAME, diag_LR_tol_reorder)
{
  DoSchurDecompTest(phist_LR,(_MT_)0.1,true);
}

TEST_F(CLASSNAME, diag_SR_tol_reorder)
{
  DoSchurDecompTest(phist_SR,(_MT_)0.03,true);
}

