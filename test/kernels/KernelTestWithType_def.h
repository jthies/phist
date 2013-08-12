
/** 
 */
template<>
class KernelTestWithType< _ST_ >
{
public:

// including these typedefs allows us to write
// things like st::sqrt in derived classes
#include "phist_std_typedefs.hpp"

bool typeImplemented_;
static const bool verbose_=true;

/** Set up method.
 * Fills internal data vector with values 1.0, 2.0 and 3.0.
 */
virtual void SetUp() {

int ierr;
SUBR(type_avail)(&ierr);
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

static _MT_ ArrayEqual(const _ST_* array, int n, int m, lidx_t lda, lidx_t stride, _ST_ value)
  {
  MT maxval=mt::zero();
  MT scal= st::abs(value);
  if (scal==mt::zero()) scal=mt::one();
  for (int i=0;i<n*stride;i+=stride)
    {
    for (int j=0; j<m;j++)
      {
      maxval=std::max(st::abs(array[j*lda+i]-value)/scal,maxval);
      }
    }
  return (MT)1.0+maxval;
  }


static _MT_ ArraysEqual(const _ST_* arr1,const _ST_* arr2, int n, int m, lidx_t lda, lidx_t stride)
  {
  _MT_ maxval=mt::zero();
  for (int i=0;i<n*stride;i+=stride)
    {
    for (int j=0;j<m;j++)
      {
      MT m = st::abs(arr1[j*lda+i]-arr2[j*lda+i]);
      MT p = (st::abs(arr1[j*lda+i])+st::abs(arr2[j*lda+i]))*(MT)0.5;
      if (p==mt::zero()) p=mt::one();
      maxval=std::max(m/p,maxval);
      }
    }
  return mt::one()+maxval;
  }

static inline _ST_ random_number() 
  {
  return (MT)(2*std::rand()-RAND_MAX)/(MT)RAND_MAX + 
         (MT)(2*std::rand()-RAND_MAX)/(MT)RAND_MAX * st::I();
  }

};
