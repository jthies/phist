
/** 
 */
template<>
class KernelTestWithType< _ST_ >
{
public:

bool typeImplemented_;
_ST_ zero_, one_, cmplx_i_;
bool verbose_;

/** Set up method.
 * Fills internal data vector with values 1.0, 2.0 and 3.0.
 */
virtual void SetUp() {
int ierr;
_SUBR_(type_avail)(&ierr);
typeImplemented_=(ierr==0);
zero_=_ZERO_;
one_=_ONE_;
cmplx_i_=_CMPLX_I_;
verbose_=true;

if (verbose_)
  {
  std::cout << "data type: ";
#ifdef _IS_COMPLEX_
  std::cout << " complex ";
#else
  std::cout << " real ";
#endif
#ifdef _IS_DOUBLE_
  std::cout << "double" << std::endl;
#else
  std::cout << "float" << std::endl;
#endif
  }
}

virtual void TearDown() {
}

_MT_ ArrayEqual(const _ST_* array, int n, int m, int lda, int stride, _ST_ value)
  {
  _MT_ maxval=(_MT_)0.0;
  _MT_ scal= std::abs(value);
  if (scal==(_MT_)0.0) scal=(_MT_)1.0;
  for (int i=0;i<n*stride;i+=stride)
    {
    for (int j=0; j<m;j++)
      {
      maxval=std::max(std::abs(array[j*lda+i]-value)/scal,maxval);
      }
    }
  return (_MT_)1.0+maxval;
  }


_MT_ ArraysEqual(const _ST_* arr1,const _ST_* arr2, int n, int m, int lda, int stride)
  {
  _MT_ maxval=(_MT_)0.0;
  for (int i=0;i<n*stride;i+=stride)
    {
    for (int j=0;j<m;j++)
      {
      _MT_ m = std::abs(arr1[j*lda+i]-arr2[j*lda+i]);
      _MT_ p = (std::abs(arr1[j*lda+i])+std::abs(arr2[j*lda+i]))*(_MT_)0.5;
      if (p==(_MT_)0.0) p=(_MT_)1.0;
      maxval=std::max(m/p,maxval);
      }
    }
  return (_MT_)1.0+maxval;
  }

inline _ST_ zero() const {return zero_;}

inline _ST_ one() const {return one_;}

inline _ST_ I() const 
  {
  return cmplx_i_;
  }

inline _ST_ random_number() const 
  {
  return (_MT_)(2*std::rand()-RAND_MAX)/(_MT_)RAND_MAX + 
         (_MT_)(2*std::rand()-RAND_MAX)/(_MT_)RAND_MAX * I();
  }

inline _MT_ eps() const {return std::numeric_limits< _MT_ >::epsilon(); }
};

