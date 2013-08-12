
/*! Test fixture. */
template<int _Nrows, int _Ncols>
class KernelTestWithSdMats<_ST_,_Nrows,_Ncols> : 
        public virtual KernelTestWithType< _ST_ >,
        public virtual KernelTest
  {

public:

  /*! Set up routine.
   */
virtual void SetUp()
  {
  KernelTest::SetUp();
  KernelTestWithType< _ST_ >::SetUp();
  if (this->typeImplemented_)
    {
    SUBR(sdMat_create)(&mat1_,this->nrows_,this->ncols_,this->comm_,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);
    SUBR(sdMat_extract_view)(mat1_,&mat1_vp_,&this->m_lda_,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);
    SUBR(sdMat_create)(&mat2_,this->nrows_,this->ncols_,this->comm_,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);
    SUBR(mvec_extract_view)(mat2_,&mat2_vp_,&this->m_lda_,&this->ierr_);
        ASSERT_EQ(0,this->ierr_);
    }
  }

  /*! Clean up.
   */
virtual void TearDown() 
  {
  if (this->typeImplemented_)
    {
    SUBR(sdMat_delete)(mat1_,&this->ierr_);
    SUBR(sdMat_delete)(mat2_,&this->ierr_);
    }
  KernelTestWithType< _ST_ >::TearDown();
  KernelTest::TearDown();
  }

static void PrintSdMat(std::ostream& os, std::string label, 
        ST* mat_vp, lidx_t lda, lidx_t stride,MPI_Comm mpi_comm)
  {
  KernelTestWithVectors<ST,_Nrows,_Ncols>::PrintVector(os, label, 
        mat_vp, (lidx_t)_Nrows, lda, stride,mpi_comm); 
  }
  
  TYPE(sdMat_ptr) mat1_, mat2_;
  _ST_ *mat1_vp_, *mat2_vp_;
  static const int nrows_=_Nrows;
  static const int ncols_=_Ncols;
  lidx_t m_lda_;
  };


template<int _Nrows, int _Ncols>
const int KernelTestWithSdMats<_ST_,_Nrows,_Ncols>::nrows_;

template<int _Nrows, int _Ncols>
const int KernelTestWithSdMats<_ST_,_Nrows,_Ncols>::ncols_;
