/*! Test fixture. */
template<int _Nglob, int _Nvec>
class KernelTestWithVectors<_ST_,_Nglob,_Nvec> : 
        public virtual KernelTestWithType< _ST_ >,
        public virtual KernelTestWithMap<_Nglob> 
  {

public:

  /*! Set up routine.
   */
virtual void SetUp()
  {
  KernelTestWithType< _ST_ >::SetUp();
  if (this->typeImplemented_)
    {
    KernelTestWithMap<_Nglob>::SetUp();
    this->nvec_=_Nvec;
    _SUBR_(mvec_create)(&vec1_,this->map_,this->nvec_,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);
    _SUBR_(mvec_extract_view)(vec1_,&vec1_vp_,&this->lda_,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);
    _SUBR_(mvec_create)(&vec2_,this->map_,this->nvec_,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);
    _SUBR_(mvec_extract_view)(vec2_,&vec2_vp_,&this->lda_,&this->ierr_);
        ASSERT_EQ(0,this->ierr_);
    stride_=1;
    }
  }

  /*! Clean up.
   */
virtual void TearDown() 
  {
  if (this->typeImplemented_)
    {
    _SUBR_(mvec_delete)(vec1_,&this->ierr_);
    _SUBR_(mvec_delete)(vec2_,&this->ierr_);
    KernelTestWithMap<_Nglob>::TearDown();
    }
  KernelTestWithType< _ST_ >::TearDown();
  }
  
  _TYPE_(mvec_ptr) vec1_, vec2_;
  _ST_ *vec1_vp_, *vec2_vp_;
  lidx_t nvec_,lda_,stride_;
  };
