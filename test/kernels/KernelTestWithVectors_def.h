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
  
  static _MT_ ColsAreNormalized(const _ST_* vec_vp, int nloc, int lda, int stride)
    {
    _MT_ res=1.0;
    // see if all columns in vec2 have 2-norm 1
    _ST_ *norms = new _ST_[nvec_];
    for (int j=0;j<nvec_;j++)
      {
      _ST_ sum=zero();
      for (int i=0;i<stride*nloc;i+=stride)
        {
        _ST_ val=vec_vp[j*lda+i];
        sum+=val*_CONJ_(val); 
        }
      norms[j]=std::sqrt(sum);
      }
    res=ArrayEqual(norms,nvec_,1,nvec_,1,one());
    delete [] norms;
    return res;
    }


  // check if vectors are mutually orthogonal after QR factorization
  static _MT_ ColsAreOrthogonal(_ST_* vec_vp, int nloc, int lda, int stride) 
    {
    _MT_ res=1.0;
    int nsums=(nloc*nvec_-nloc)/2;
      _ST_ sums[nsums];
      int k=0;
      for (int j1=0;j1<nvec_;j1++)
      for (int j2=j1+1;j2<nvec_;j2++)
        {
        _ST_ sum=zero();
        for (int i=0;i<stride*nloc;i+=stride)
          {
          _ST_ val1=vec_vp[j1*lda+i];
          _ST_ val2=vec_vp[j2*lda+i];
          sum+=val1*_CONJ_(val2); 
          }
        sums[k++]=sum;
        }
      res=ArrayEqual(sums,nsums,1,nvec_,1,zero());
      return res;
      }

  _TYPE_(mvec_ptr) vec1_, vec2_;
  _ST_ *vec1_vp_, *vec2_vp_;
  static const int nvec_=_Nvec;
  lidx_t lda_, stride_;
  };

template<int n, int nvec>
const int KernelTestWithVectors<_ST_,n,nvec>::nvec_;

