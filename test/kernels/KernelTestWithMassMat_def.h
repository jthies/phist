#include "../tools/TestHelpers.h"

using namespace phist::testing;

/*! Test fixture. */
template<phist_gidx _Nglob>
class KernelTestWithMassMat<_ST_, _Nglob> : public virtual TestWithType< _ST_ >,
                                            public virtual KernelTestWithMap<_Nglob>,
                                            public virtual KernelTest
{

  public:

    // if the ctx argument is NULL, the sparse matrix B will be constructed
    // without a context and its domain map will be used to replace the map_
    // of the base class (cf. KernelTestWithMap::SetUpTestCaseWithMap).
    // Otherwise (if the context comes, sy, from a KernelTestWithSparseMat),
    // the mass matrix will be made compatible to that class' sparseMat A_.
    static void SetUpTestCase(phist_const_context_ptr ctx=NULL)
    {
      TestWithType<_ST_>::SetUpTestCase();
      KernelTest::SetUpTestCase();
      
      my_map_=false;

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
        if (ctx==NULL)
        {
          SUBR(sparseMat_create_fromRowFunc)(&B_,comm_,gnrows,gnrows,3,&phist::testing::PHIST_TG_PREFIX(hpd_tridiag),NULL,&iflag_);
          ASSERT_EQ(0,iflag_);
          phist_const_map_ptr domain_map=NULL;
          SUBR(sparseMat_get_domain_map)(B_,&domain_map,&iflag_);
          ASSERT_EQ(0,iflag_);
          KernelTestWithMap<_Nglob>::SetUpTestCaseWithMap(domain_map);
          my_map_=true;
        }
        else
        {
          SUBR(sparseMat_create_fromRowFuncAndContext)(&B_,ctx,3,&phist::testing::PHIST_TG_PREFIX(hpd_tridiag),NULL,&iflag_);
          ASSERT_EQ(0,iflag_);
        }

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
        if (my_map_)
        {
          KernelTestWithMap<_Nglob>::TearDownTestCase();
        }
      }
    }

    static TYPE(sparseMat_ptr) B_;
    static bool my_map_;
};


template<phist_gidx _Nglob>
TYPE(sparseMat_ptr) KernelTestWithMassMat<_ST_, _Nglob> :: B_;

template<phist_gidx _Nglob>
bool KernelTestWithMassMat<_ST_, _Nglob> :: my_map_;
