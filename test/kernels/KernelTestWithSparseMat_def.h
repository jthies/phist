/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "../tools/TestHelpers.h"

using namespace phist::testing;

extern "C" {

// prototype for a useful function from the driver_utils, we can't include
// the header here because it can only be included once after a phist_gen_X header.
void SUBR(create_matrix)(TYPE(sparseMat_ptr)* mat, phist_const_comm_ptr comm,
        const char* problem, int* iflag);

} //extern "C"

namespace phist {
namespace testing {
int PHIST_TG_PREFIX(idfunc)(ghost_gidx row, ghost_lidx *len, ghost_gidx* cols, void* vval, void *arg);
}
}

/*! Test fixture. */
template<phist_gidx _Nglob, phist_gidx _Mglob, MATNAME_ENUM _MatName, int _multipleDefinitionCounter>
class KernelTestWithSparseMat<_ST_, _Nglob, _Mglob, _MatName, _multipleDefinitionCounter> : 
                                                                                    public virtual TestWithType< _ST_ >,
                                                                                    public virtual KernelTest,
                                                                                    public virtual KernelTestWithMap<_Nglob>
{

  public:

    static void SetUpTestCase(int sparseMatCreateFlag, phist_const_context_ptr ctx=NULL)
    {
      TestWithType<_ST_>::SetUpTestCase();
      KernelTest::SetUpTestCase();

      // preset problemTooSmall_ to abort appropriately
      bool problemTooSmall = mpi_size_ > _Nglob;
      
      iflag_ = sparseMatCreateFlag;

      if( typeImplemented_ && !problemTooSmall )
      {
        iflag_=0;
        // read matrix first
        if (_MatName==MATNAME_IDFUNC)
        {
          // this is not handled by create_matrix, so it gets an extra treatment
          phist_gidx ini_dim[2];
          ini_dim[0]=_Nglob;
          ini_dim[1]=_Mglob;
          PHIST_TG_PREFIX(idfunc)(-1,NULL,ini_dim,NULL,NULL);
          SUBR(sparseMat_create_fromRowFunc)(&A_,comm_,_Nglob,_Mglob,1,&PHIST_TG_PREFIX(idfunc),NULL,&iflag_);
          // "uninitialize"
          PHIST_TG_PREFIX(idfunc)(-2,NULL,NULL,NULL,NULL);
        }
        else if (_MatName>=MATNAME_hpd_tridiag && _MatName<=MATNAME_nhid_tridiag_ainv)
        {
          phist_sparseMat_rowFunc rowfunc=NULL;
          if (_MatName==MATNAME_hpd_tridiag) rowfunc=phist::testing::PHIST_TG_PREFIX(hpd_tridiag);
          else if (_MatName==MATNAME_lapl_tridiag) rowfunc=phist::testing::PHIST_TG_PREFIX(lapl_tridiag);
          else if (_MatName==MATNAME_hid_tridiag) rowfunc=phist::testing::PHIST_TG_PREFIX(hid_tridiag);
          else if (_MatName==MATNAME_nhpd_tridiag) rowfunc=phist::testing::PHIST_TG_PREFIX(nhpd_tridiag);
          else if (_MatName==MATNAME_nhid_tridiag) rowfunc=phist::testing::PHIST_TG_PREFIX(nhid_tridiag);
          else if (_MatName==MATNAME_hpd_tridiag_ainv) rowfunc=phist::testing::PHIST_TG_PREFIX(hpd_tridiag_ainv);
          else if (_MatName==MATNAME_lapl_tridiag_ainv) rowfunc=phist::testing::PHIST_TG_PREFIX(lapl_tridiag_ainv);
          else if (_MatName==MATNAME_hid_tridiag_ainv) rowfunc=phist::testing::PHIST_TG_PREFIX(hid_tridiag_ainv);
          else if (_MatName==MATNAME_nhpd_tridiag_ainv) rowfunc=phist::testing::PHIST_TG_PREFIX(nhpd_tridiag_ainv);
          else if (_MatName==MATNAME_nhid_tridiag_ainv) rowfunc=phist::testing::PHIST_TG_PREFIX(nhid_tridiag_ainv);
          else
          {
            iflag_=-99;
          }
          ASSERT_EQ(0,iflag_);
          // create B_ as a tridiagonal hpd matrix
          ghost_gidx gnrows=_Nglob;
          // initialize rowFunc
          iflag_=rowfunc(-1,NULL,&gnrows,NULL,NULL);
          ASSERT_EQ(0,iflag_);
          if (ctx==NULL)
          {
            SUBR(sparseMat_create_fromRowFunc)(&A_,comm_,_Nglob,_Mglob,3,rowfunc,NULL,&iflag_);
            ASSERT_EQ(0,iflag_);
          }
          else
          {
            // check that the given context is compatible
            SUBR(sparseMat_create_fromRowFuncAndContext)(&A_,ctx,3,rowfunc,NULL,&iflag_);
            ASSERT_EQ(0,iflag_);
          }

          // "uninitialize"
          rowfunc(-2,NULL,NULL,NULL,NULL);
        }
        else if (MatNameEnumIsMatFunc(_MatName)==false)
        {
          SUBR(read_mat)(MatNameEnumToStr(_MatName),comm_,_Nglob,_Mglob,&A_,&iflag_);
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
        SUBR(sparseMat_get_range_map)(A_,&range_map,&iflag_);
        ASSERT_EQ(0,iflag_);

        // check the size of the maps
        phist_gidx map_nglob = 0;
        phist_map_get_global_length(range_map,&map_nglob,&iflag_);
        ASSERT_EQ(0,iflag_);
        ASSERT_EQ(_Nglob,map_nglob);

        phist_gidx map_mglob = 0;
        phist_map_get_global_length(domain_map,&map_mglob,&iflag_);
        ASSERT_EQ(0,iflag_);
        ASSERT_EQ(_Mglob,map_mglob);
        
        // check that the range- and domain map are the same. Our tests dassume
        // the matrix is square and symmetrically permuted.
        if (_Nglob==_Mglob)
        {
          phist_maps_compatible(range_map,domain_map,&iflag_);
          ASSERT_EQ(0,iflag_);
        }
        
        SUBR(sparseMat_get_context)(A_,&context_,&iflag_);
        ASSERT_EQ(0,iflag_);

        // now setup the map
        KernelTestWithMap<_Nglob>::SetUpTestCaseWithMap(range_map);
      }
    }

    virtual void SetUp()
    {
      KernelTest::SetUp();
      KernelTestWithMap<_Nglob>::SetUp();

      if( this->typeImplemented_ && !this->problemTooSmall_ )
      {
        phist_const_map_ptr map = NULL;
        SUBR(sparseMat_get_range_map)(A_,&map,&this->iflag_);
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


template<phist_gidx _Nglob,phist_gidx _Mglob,MATNAME_ENUM _MatName,int _multipleDefinitionCounter>
TYPE(sparseMat_ptr) KernelTestWithSparseMat<_ST_,_Nglob,_Mglob, _MatName,_multipleDefinitionCounter> :: A_;
template<phist_gidx _Nglob,phist_gidx _Mglob,MATNAME_ENUM _MatName,int _multipleDefinitionCounter>
phist_const_context_ptr KernelTestWithSparseMat<_ST_,_Nglob,_Mglob,_MatName,_multipleDefinitionCounter> :: context_;
