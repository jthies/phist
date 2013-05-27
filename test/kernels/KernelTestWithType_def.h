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

/** Set up method.
 * Fills internal data vector with values 1.0, 2.0 and 3.0.
 */
virtual void SetUp() {
int ierr;
_SUBR_(type_avail)(&ierr);
typeImplemented_=(ierr==0);
zero_=_ZERO_;
one_=_ONE_;
}

virtual void TearDown() {
}

AssertionResult ArrayEqual(const _ST_* array, int n, _ST_ value)
  {
  _MT_ sum=(_MT_)0.0;
  for (int i=0;i<n;i++)
    {
    sum+=std::abs(array[i]-value);
    }
  return (sum==(_MT_)0.0 ? AssertionSuccess() : AssertionFailure());
  }


AssertionResult ArraysEqual(const _ST_* arr1,const _ST_* arr2, int n)
  {
  _MT_ sum=(_MT_)0.0;
  for (int i=0;i<n;i++)
    {
    sum+=std::abs(arr1[i]-arr2[i]);
    }
  return (sum==(_MT_)0.0 ? AssertionSuccess() : AssertionFailure());
  }

inline _ST_ zero() {return zero_;}

inline _ST_ one() {return one_;}

_ST_ random_number() {return (_MT_)std::rand()/(_MT_)RAND_MAX + (_MT_)std::rand()/(_MT_)RAND_MAX * _Complex_I;}


};

