#ifndef CLASSNAME
#error "file not included correctly"
#endif
#ifndef MATNAME
#error "file not included correctly"
#endif
#ifndef _N_
#error "file not included correctly"
#endif

/*! Test fixure. */
class CLASSNAME: public KernelTestWithVectors<_ST_,_N_,_M_>
{

  public:

    typedef KernelTestWithVectors<_ST_,_N_,_M_> VTest;

    //! mvec/sdMat sizes
    static const int n_=_N_;
    static const int m_=_M_;

    //! Set up routine.
    virtual void SetUp()
    {
      int ierr;
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

        SUBR(op_wrap_crsMat)(opA_,A_,&ierr);
        ASSERT_EQ(0,ierr);
      }
    }

    //! Clean up.
    virtual void TearDown()
    {
      int ierr;
      if( typeImplemented_ )
      {

        delete opA_;

        SUBR(crsMat_delete)(A_,&ierr_);
        ASSERT_EQ(0,ierr_);
      }

      VTest::TearDown();
    }

    TYPE(op_ptr) opA_;

  protected:
  
    TYPE(crsMat_ptr) A_;

/*
    // ========================= the actual arnoldi test =========================
    void doArnoldiTest(TYPE(const_op_ptr) opA)
    {
      int ierr;
      if( typeImplemented_ )
      {
        // run simple_arnoldi
        SUBR(simple_arnoldi)(opA,NULL,v0_,V_,NULL,H_,m_,&ierr);
        ASSERT_EQ(0,ierr);

        // check orthogonality of V_
        ASSERT_REAL_EQ(mt::one(),VTest::ColsAreNormalized(V_vp_,nloc_,ldaV_,stride_,mpi_comm_));
        ASSERT_REAL_EQ(mt::one(),VTest::ColsAreOrthogonal(V_vp_,nloc_,ldaV_,stride_,mpi_comm_));

        // calculate A*V(:,1:m)
        opA->apply(st::one(),opA->A,Vm_,st::zero(),AV_,&ierr);
        ASSERT_EQ(0,ierr);
        // calculate V(:,1:m+1)*H(1:m+1,1:m)
        SUBR(mvec_times_sdMat)(st::one(),V_,H_,st::zero(),VH_,&ierr);
        ASSERT_EQ(0,ierr);

        // calculate AV_' := AV_ - VH_
        SUBR(mvec_add_mvec)(-st::one(),VH_,st::one(),AV_,&ierr);
        ASSERT_EQ(0,ierr);
        _MT_ vnorm[_M_];
        SUBR(mvec_norm2)(AV_,vnorm,&ierr);
        ASSERT_EQ(0,ierr);

        // check AV_' = AV_ - VH_ == 0
        for(int i = 0; i < _M_; i++)
          ASSERT_NEAR(mt::zero(),vnorm[i],100*mt::eps());
      }
    }
    */
};


TEST_F(CLASSNAME, some_test) 
{
/*
  int ierr;
  SUBR(mvec_put_value)(v0_,st::one(),&ierr);
  ASSERT_EQ(0,ierr);
  doArnoldiTest(opAeye_);
*/
}

