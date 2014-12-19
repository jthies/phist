#ifndef PHIST_KERNEL_TEST_WITH_MAP_H
#define PHIST_KERNEL_TEST_WITH_MAP_H

#include "KernelTest.h"

/** 
 */
template<gidx_t _Nglob>
class KernelTestWithMap: public virtual KernelTest
{
public:

/** Set up method.
 */
virtual void SetUp()
{
  KernelTest::SetUp();
  map_ptr_t map;
  phist_map_create(&map,comm_,nglob_,&iflag_);
  ASSERT_EQ(0,iflag_);
  map_ = map;
  deleteMap_ = true;
  phist_map_get_local_length(map_,&nloc_,&iflag_);
  ASSERT_EQ(0,iflag_);
}

/** clean up the global test environment */
virtual void TearDown()
{
  if( deleteMap_ )
  {
    phist_map_delete((map_ptr_t)map_,&iflag_);
    ASSERT_EQ(0,iflag_);
  }
  KernelTest::TearDown();
}

/** replace the map, e.g. when we need to use the map of a sparseMat */
virtual void replaceMap(const_map_ptr_t map)
{
  if( deleteMap_ )
  {
    phist_map_delete((map_ptr_t)map_,&iflag_);
    ASSERT_EQ(0,iflag_);
    deleteMap_ = false;
  }
  map_ = map;
  phist_map_get_local_length(map_,&nloc_,&iflag_);
  ASSERT_EQ(0,iflag_);

}

static const gidx_t nglob_=_Nglob;
lidx_t nloc_;
const_map_ptr_t map_;
bool deleteMap_;
};

template<gidx_t n>
const gidx_t KernelTestWithMap<n>::nglob_;

#endif
