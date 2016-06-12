#include "../tools/TestHelpers.h"

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

    static void SetUpTestCase()
    {
      TestWithType<_ST_>::SetUpTestCase();
      KernelTest::SetUpTestCase();

      // preset problemTooSmall_ to abort appropriately
      bool problemTooSmall = mpi_size_ > _Nglob;

      if( typeImplemented_ && !problemTooSmall )
      {
        // read matrix first
        if (MatNameEnumIsMatFunc(_MatName))
        {
          SUBR(sparseMat_create_fromRowFunc)(&A_,comm_,_Nglob,_Nglob,1,&SUBR(idfunc),NULL,&iflag_);
        }
        else
        {
          SUBR(read_mat)(MatNameEnumToStr(_MatName),comm_,_Nglob,&A_,&iflag_);
        }
        if( iflag_ != 0 )
        {
          SUBR(create_matrix)(&A_,comm_,MatNameEnumToStr(_MatName),&iflag_);
        }
        ASSERT_EQ(0,iflag_);
        ASSERT_TRUE(A_ != NULL);

        phist_const_map_ptr map = NULL;
        SUBR(sparseMat_get_domain_map)(A_,&map,&iflag_);
        ASSERT_EQ(0,iflag_);

        // check the size of the map
        phist_gidx map_nglob = 0;
        phist_map_get_global_length(map,&map_nglob,&iflag_);
        ASSERT_EQ(0,iflag_);
        ASSERT_EQ(_Nglob,map_nglob);

        // now setup the map
        KernelTestWithMap<_Nglob>::SetUpTestCaseWithMap(map);
      }
    }

    virtual void SetUp()
    {
      KernelTest::SetUp();
      KernelTestWithMap<_Nglob>::SetUp();

// with GHOST the map gets recreqted on the fly, so we cannot easily the identity of two maps
#ifndef PHIST_KERNEL_LIB_GHOST
      if( this->typeImplemented_ && !this->problemTooSmall_ )
      {
        phist_const_map_ptr map = NULL;
        SUBR(sparseMat_get_domain_map)(A_,&map,&this->iflag_);
        ASSERT_EQ(0,this->iflag_);
        ASSERT_EQ(map, this->map_);
      }
#endif
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
};


template<phist_gidx _Nglob, MATNAME_ENUM _MatName, int _multipleDefinitionCounter>
TYPE(sparseMat_ptr) KernelTestWithSparseMat<_ST_, _Nglob, _MatName, _multipleDefinitionCounter> :: A_;
