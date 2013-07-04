
/** 
 */
template<>
class KernelTestWithType< _ST_ >
{
public:

bool typeImplemented_;
static const bool verbose_=true;

/** Set up method.
 * Fills internal data vector with values 1.0, 2.0 and 3.0.
 */
virtual void SetUp() {
int ierr;
_SUBR_(type_avail)(&ierr);
typeImplemented_=(ierr==0);

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

static _MT_ ArrayEqual(const _ST_* array, int n, int m, int lda, int stride, _ST_ value)
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


static _MT_ ArraysEqual(const _ST_* arr1,const _ST_* arr2, int n, int m, int lda, int stride)
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

inline static _ST_ zero() {return _ZERO_;}

inline static _ST_ one() {return _ONE_;}

static inline _ST_ I()
  {
  return _CMPLX_I_;
  }

static inline _ST_ random_number() 
  {
  return (_MT_)(2*std::rand()-RAND_MAX)/(_MT_)RAND_MAX + 
         (_MT_)(2*std::rand()-RAND_MAX)/(_MT_)RAND_MAX * I();
  }

static inline _MT_ eps() {return std::numeric_limits< _MT_ >::epsilon(); }
};

