#ifndef PHIST_KERNEL_TEST_WITH_MAP_H
#define PHIST_KERNEL_TEST_WITH_MAP_H

#include "KernelTest.h"

/** 
 */
template<phist_gidx _Nglob>
class KernelTestWithMap: public virtual KernelTest
{
public:

static void SetUpTestCase()
{
  KernelTest::SetUpTestCase();
  problemTooSmall_ = mpi_size_ > nglob_;
  // this function may be called by multiple base classes, e.g.
  // by KernelTestWithVectors and KernelTestWithSparseMat, so  
  // do not re-build the map if it has been called already
  if( defaultMap_ == NULL)
  {
    phist_map_ptr map;
    iflag_=PHIST_SPARSEMAT_QUIET;
    phist_map_create(&map,comm_,nglob_,&iflag_);
    EXPECT_EQ(0,iflag_);
    defaultMap_=map;
    map_ = map;
    phist_map_get_local_length(map_,&nloc_,&iflag_);
    EXPECT_EQ(0,iflag_);
    phist_context_create(&defaultContext_,map_,map_,map_,&iflag_);
    EXPECT_EQ(0,iflag_);
  }
}

static void SetUpTestCaseWithMap(phist_const_map_ptr map)
{
  EXPECT_TRUE(map != NULL);
  // prevent getting called multiple times
  EXPECT_EQ(NULL,map_);

  SetUpTestCase(); // creates defaultMap_, determines problemTooSmall_ etc

  map_ = map;
  phist_map_get_local_length(map_,&nloc_,&iflag_);
  EXPECT_EQ(0,iflag_);
  
  // check if the given map is compatible with the default linear map. It may be permuted and repartitioned,
  // so allow all values >=0
  phist_maps_compatible(map_,defaultMap_,&iflag_);
  EXPECT_TRUE(iflag_>=0);
}

static void TearDownTestCase()
{
  if (defaultMap_!=NULL)
  {
    phist_context_delete(defaultContext_,&iflag_);
    EXPECT_EQ(0,iflag_);
    phist_map_delete(defaultMap_,&iflag_);
    EXPECT_EQ(0,iflag_);
    defaultMap_=NULL;
  }
  map_ = NULL; // this either belongs to another object or points to defaultMap_
}

static const phist_gidx nglob_=_Nglob;
static phist_lidx nloc_;
static phist_const_map_ptr map_;
static phist_map_ptr defaultMap_;
static phist_context_ptr defaultContext_;
static bool problemTooSmall_;
};

template<phist_gidx _Nglob>
const phist_gidx KernelTestWithMap<_Nglob>::nglob_;

template<phist_gidx _Nglob>
phist_lidx KernelTestWithMap<_Nglob>::nloc_ = 0;

template<phist_gidx _Nglob>
phist_const_map_ptr KernelTestWithMap<_Nglob>::map_ = NULL;

template<phist_gidx _Nglob>
phist_map_ptr KernelTestWithMap<_Nglob>::defaultMap_ = NULL;

template<phist_gidx _Nglob>
phist_context_ptr KernelTestWithMap<_Nglob>::defaultContext_ = NULL;

template<phist_gidx _Nglob>
bool KernelTestWithMap<_Nglob>::problemTooSmall_ = false;

#endif
