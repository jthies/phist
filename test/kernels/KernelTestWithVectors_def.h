
/*! Test fixture. */
template<int _Nglob, int _Nvec>
class KernelTestWithVectors<_ST_,_Nglob,_Nvec> : 
        public KernelTestWithType< _ST_ >,
        public KernelTestWithMap<_Nglob> 
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
    ASSERT_EQ(this->ierr_, 0);
    _SUBR_(mvec_extract_view)(vec1_,&vec1_vp_,&this->lda_,&this->ierr_);
    ASSERT_EQ(this->ierr_, 0);
    _SUBR_(mvec_create)(&vec2_,this->map_,this->nvec_,&this->ierr_);
    ASSERT_EQ(this->ierr_, 0);
    _SUBR_(mvec_extract_view)(vec2_,&vec2_vp_,&this->lda_,&this->ierr_);
        ASSERT_EQ(this->ierr_, 0);
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
  int nvec_,lda_;
  };


