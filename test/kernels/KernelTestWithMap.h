#ifndef PHIST_KERNEL_TEST_WITH_MAP_H
#define PHIST_KERNEL_TEST_WITH_MAP_H

#include "KernelTest.h"

/** 
 */
template<int _Nglob>
class KernelTestWithMap: public virtual KernelTest
{
public:

/** Set up method.
 * Fills internal data vector with values 1.0, 2.0 and 3.0.
 */
virtual void SetUp() {
KernelTest::SetUp();
nloc_=nglob_; // TODO - parallel testing not implemented
phist_map_create(&map_,comm_,nglob_,&ierr_);
ASSERT_EQ(0,ierr_);
}

static const int nglob_=_Nglob;
int nloc_;
map_ptr_t map_;


};

template<int n>
const int KernelTestWithMap<n>::nglob_;

#endif
