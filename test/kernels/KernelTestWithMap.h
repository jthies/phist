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
 */
virtual void SetUp() {
KernelTest::SetUp();
phist_map_create(&map_,comm_,nglob_,&ierr_);
ASSERT_EQ(0,ierr_);
phist_map_get_local_length(map_,&nloc_,&ierr_);
ASSERT_EQ(0,ierr_);
}

/** clean up the global test environment */
virtual void TearDown()
  {
  phist_map_delete(map_,&ierr_);
ASSERT_EQ(0,ierr_);
  KernelTest::TearDown();
  }

static const gidx_t nglob_=_Nglob;
lidx_t nloc_;
map_ptr_t map_;


};

template<int n>
const gidx_t KernelTestWithMap<n>::nglob_;

#endif
