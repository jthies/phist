#include "../tools/TestHelpers.h"

using namespace phist::testing;

extern "C" {

// prototype for a useful function from the driver_utils, we can't include
// the header here because it can only be included once after a phist_gen_X header.
void SUBR(create_matrix)(TYPE(sparseMat_ptr)* mat, phist_const_comm_ptr comm,
        const char* problem, int* iflag);

} //extern "C"

/*! Test fixture. */
template<phist_gidx _Nglob, MATNAME_ENUM _MatName, int _multipleDefinitionCounter>
class KernelTestWithSparseMat<_ST_, _Nglob, _MatName, _multipleDefinitionCounter> : public virtual TestWithType< _ST_ >,
                                                                                    public virtual KernelTest,
                                                                                    public virtual KernelTestWithMap<_Nglob> 
{

  public:

    static void SetUpTestCase(int sparseMatCreateFlag)
    {
      TestWithType<_ST_>::SetUpTestCase();
      KernelTest::SetUpTestCase();

      // preset problemTooSmall_ to abort appropriately
      bool problemTooSmall = mpi_size_ > _Nglob;
      
      iflag_ = sparseMatCreateFlag;

      if( typeImplemented_ && !problemTooSmall )
      {
        // read matrix first
        if (_MatName==MATNAME_IDFUNC)
        {
          // this is not handled by create_matrix, so it gets an extra treatment
          SUBR(sparseMat_create_fromRowFunc)(&A_,comm_,_Nglob,_Nglob,1,&SUBR(idfunc),NULL,&iflag_);
        }
        else if (MatNameEnumIsMatFunc(_MatName)==false)
        {
          SUBR(read_mat)(MatNameEnumToStr(_MatName),comm_,_Nglob,&A_,&iflag_);
        }
        else
        {
          SUBR(create_matrix)(&A_,comm_,MatNameEnumToStr(_MatName),&iflag_);
        }
        ASSERT_EQ(0,iflag_);
        ASSERT_TRUE(A_ != NULL);

        phist_const_map_ptr domain_map = NULL;
        SUBR(sparseMat_get_domain_map)(A_,&domain_map,&iflag_);
        ASSERT_EQ(0,iflag_);

        phist_const_map_ptr range_map = NULL;
        SUBR(sparseMat_get_domain_map)(A_,&range_map,&iflag_);
        ASSERT_EQ(0,iflag_);

        // check the size of the map
        phist_gidx map_nglob = 0;
        phist_map_get_global_length(domain_map,&map_nglob,&iflag_);
        ASSERT_EQ(0,iflag_);
        ASSERT_EQ(_Nglob,map_nglob);
        
        // check that the range- and domain map are the same. Our tests dassume
        // the matrix is square and symmetrically permuted.
        phist_maps_compatible(range_map,domain_map,&iflag_);
        ASSERT_EQ(0,iflag_);
        
        SUBR(sparseMat_get_context)(A_,&context_,&iflag_);
        ASSERT_EQ(0,iflag_);

        // now setup the map
        KernelTestWithMap<_Nglob>::SetUpTestCaseWithMap(domain_map);
      }
    }

    virtual void SetUp()
    {
      KernelTest::SetUp();
      KernelTestWithMap<_Nglob>::SetUp();

      if( this->typeImplemented_ && !this->problemTooSmall_ )
      {
        phist_const_map_ptr map = NULL;
        SUBR(sparseMat_get_domain_map)(A_,&map,&this->iflag_);
        ASSERT_EQ(0,this->iflag_);
        // make sure that the domain map of the matrix is the base map of this test class
        phist_maps_compatible(map,KernelTestWithMap<_Nglob>::map_,&iflag_);
        ASSERT_EQ(0,iflag_);
      }
    }

    virtual void TearDown() 
    {
      KernelTestWithMap<_Nglob>::TearDown();
      KernelTest::TearDown();
    }

    static void TearDownTestCase()
    {
      // respect teardown order!
      KernelTestWithMap<_Nglob>::TearDownTestCase();

      if( A_ != NULL )
      {
        SUBR(sparseMat_delete)(A_,&iflag_);
        ASSERT_EQ(0,iflag_);
        A_ = NULL;
      }

      KernelTest::TearDownTestCase();
    }

    static TYPE(sparseMat_ptr) A_;
    static phist_const_context_ptr context_;
};


template<phist_gidx _Nglob,MATNAME_ENUM _MatName,int _multipleDefinitionCounter>
TYPE(sparseMat_ptr) KernelTestWithSparseMat<_ST_,_Nglob,_MatName,_multipleDefinitionCounter> :: A_;
template<phist_gidx _Nglob,MATNAME_ENUM _MatName,int _multipleDefinitionCounter>
phist_const_context_ptr KernelTestWithSparseMat<_ST_,_Nglob,_MatName,_multipleDefinitionCounter> :: context_;
