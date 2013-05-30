#include <cstdlib>
#include "gtest/gtest.h"

using namespace ::testing;

/** 
 */
template<>
class KernelTestWithType< _ST_ >
{
public:

bool typeImplemented_;
_ST_ zero_, one_;
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
verbose_=true;
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

inline _ST_ zero() {return zero_;}

inline _ST_ one() {return one_;}

_ST_ random_number() {return (_MT_)std::rand()/(_MT_)RAND_MAX + (_MT_)std::rand()/(_MT_)RAND_MAX * _Complex_I;}

};

