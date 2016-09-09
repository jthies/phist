#include "../tools/TestHelpers.h"

using namespace phist::testing;

/*! Test fixture. */
template<phist_gidx _Nglob>
class KernelTestWithMassMat<_ST_, _Nglob> : public virtual TestWithType< _ST_ >,
                                            public virtual KernelTest
{

  public:

    static void SetUpTestCase(phist_const_map_ptr map)
    {
      TestWithType<_ST_>::SetUpTestCase();
      KernelTest::SetUpTestCase();

      // preset problemTooSmall_ to abort appropriately
      bool problemTooSmall = mpi_size_ > _Nglob;
                        
      if( typeImplemented_ && !problemTooSmall )
      {
        iflag_ = 0;
      // create B_ as a tridiagonal hpd matrix
      ghost_gidx gnrows=_Nglob;
      // initialize rowFunc
      iflag_=phist::testing::PHIST_TG_PREFIX(hpd_tridiag)(-1,NULL,&gnrows,NULL,NULL);
      ASSERT_EQ(0,iflag_);
      SUBR(sparseMat_create_fromRowFuncAndMap)(&B_,map,3,&phist::testing::PHIST_TG_PREFIX(hpd_tridiag),NULL,&iflag_);
      ASSERT_EQ(0,iflag_);

        ASSERT_TRUE(B_ != NULL);
      }
    }

    virtual void SetUp()
    {
    }

    virtual void TearDown() 
    {
    }

    static void TearDownTestCase()
    {
      if( B_ != NULL )
      {
        SUBR(sparseMat_delete)(B_,&iflag_);
        ASSERT_EQ(0,iflag_);
        B_ = NULL;
      }
    }

    static TYPE(sparseMat_ptr) B_;
};


template<phist_gidx _Nglob>
TYPE(sparseMat_ptr) KernelTestWithMassMat<_ST_, _Nglob> :: B_;
