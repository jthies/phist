#ifndef PHIST_KERNEL_TEST_WITH_MAP_H
#define PHIST_KERNEL_TEST_WITH_MAP_H

#include "KernelTest.h"

/** 
 */
template<phist_gidx _Nglob>
class KernelTestWithMap: public virtual KernelTest
{
public:

static void SetUpTestCaseWithMap(phist_const_map_ptr map)
{
  EXPECT_TRUE(map != NULL);
  // prevent getting called multiple times
  EXPECT_EQ(NULL,map_);
  EXPECT_EQ(false,staticDeleteMap_);

  KernelTest::SetUpTestCase();

  problemTooSmall_ = mpi_size_ > nglob_;
  map_ = map;
  staticDeleteMap_ = false;
  phist_map_get_local_length(map_,&nloc_,&iflag_);
  EXPECT_EQ(0,iflag_);
}

static void SetUpTestCase()
{
  if( map_ == NULL )
  {
    // prevent getting called multiple times
    EXPECT_EQ(false,staticDeleteMap_);

    KernelTest::SetUpTestCase();

    phist_map_ptr map;
    problemTooSmall_ = mpi_size_ > nglob_;
    iflag_=PHIST_SPARSEMAT_QUIET;
    phist_map_create(&map,comm_,nglob_,&iflag_);
    EXPECT_EQ(0,iflag_);
    map_ = map;
    staticDeleteMap_ = true;
    phist_map_get_local_length(map_,&nloc_,&iflag_);
    EXPECT_EQ(0,iflag_);
  }
}

static void TearDownTestCase()
{
  if( staticDeleteMap_ )
  {
    EXPECT_TRUE(map_ != NULL);
    phist_map_delete((phist_map_ptr)map_,&iflag_);
    EXPECT_EQ(0,iflag_);

  }
  staticDeleteMap_ = false;
  map_ = NULL;
}

static const phist_gidx nglob_=_Nglob;
static phist_lidx nloc_;
static phist_const_map_ptr map_;
static bool staticDeleteMap_;
static bool problemTooSmall_;
};

template<phist_gidx _Nglob>
const phist_gidx KernelTestWithMap<_Nglob>::nglob_;

template<phist_gidx _Nglob>
phist_lidx KernelTestWithMap<_Nglob>::nloc_ = 0;

template<phist_gidx _Nglob>
phist_const_map_ptr KernelTestWithMap<_Nglob>::map_ = NULL;

template<phist_gidx _Nglob>
bool KernelTestWithMap<_Nglob>::staticDeleteMap_ = false;

template<phist_gidx _Nglob>
bool KernelTestWithMap<_Nglob>::problemTooSmall_ = false;

#endif
