#ifndef PHIST_KERNEL_TEST_WITH_MAP_H
#define PHIST_KERNEL_TEST_WITH_MAP_H

#include "KernelTest.h"

/** 
 */
template<gidx_t _Nglob>
class KernelTestWithMap: public virtual KernelTest
{
public:

static void SetUpTestCaseWithMap(const_map_ptr_t map)
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

    map_ptr_t map;
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
    phist_map_delete((map_ptr_t)map_,&iflag_);
    EXPECT_EQ(0,iflag_);

  }
  staticDeleteMap_ = false;
  map_ = NULL;
}

static const gidx_t nglob_=_Nglob;
static lidx_t nloc_;
static const_map_ptr_t map_;
static bool staticDeleteMap_;
static bool problemTooSmall_;
};

template<gidx_t _Nglob>
const gidx_t KernelTestWithMap<_Nglob>::nglob_;

template<gidx_t _Nglob>
lidx_t KernelTestWithMap<_Nglob>::nloc_ = 0;

template<gidx_t _Nglob>
const_map_ptr_t KernelTestWithMap<_Nglob>::map_ = NULL;

template<gidx_t _Nglob>
bool KernelTestWithMap<_Nglob>::staticDeleteMap_ = false;

template<gidx_t _Nglob>
bool KernelTestWithMap<_Nglob>::problemTooSmall_ = false;

#endif
