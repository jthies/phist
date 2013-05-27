#ifndef PHIST_KERNEL_TEST_WITH_MAP_H
#define PHIST_KERNEL_TEST_WITH_MAP_H

#include "KernelTest.h"

/** 
 */
template<int _Nglob>
class KernelTestWithMap: public KernelTest
{
public:

/** Set up method.
 * Fills internal data vector with values 1.0, 2.0 and 3.0.
 */
virtual void SetUp() {
KernelTest::SetUp();
nglob_=_Nglob;
nloc_=nglob_; // TODO - parallel testing not implemented
phist_map_create(&map_,comm_,nglob_,&ierr_);
ASSERT_EQ(ierr_,0);
}

int nglob_, nloc_;
map_ptr_t map_;


};

#endif
