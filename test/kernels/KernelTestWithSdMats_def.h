
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
    SUBR(sdMat_extract_view)(mat2_,&mat2_vp_,&this->m_lda_,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);
    SUBR(sdMat_create)(&mat3_,this->nrows_,this->ncols_,this->comm_,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);
    SUBR(sdMat_extract_view)(mat3_,&mat3_vp_,&this->m_lda_,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);
    SUBR(sdMat_create)(&mat4_,this->nrows_,this->ncols_,this->comm_,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);
    SUBR(sdMat_extract_view)(mat4_,&mat4_vp_,&this->m_lda_,&this->ierr_);
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
    SUBR(sdMat_delete)(mat3_,&this->ierr_);
    SUBR(sdMat_delete)(mat4_,&this->ierr_);
    }
  KernelTestWithType< _ST_ >::TearDown();
  KernelTest::TearDown();
  }

static void PrintSdMat(std::ostream& os, std::string label, 
        ST* mat_vp, lidx_t lda, lidx_t stride,MPI_Comm mpi_comm)
  {
    int rank=0, np=1;
#ifdef PHIST_HAVE_MPI
    MPI_Comm_rank(mpi_comm,&rank);
    MPI_Comm_size(mpi_comm,&np);
#endif    
    if (rank==0)
      {
      os << std::endl<<label <<"="<<std::endl;
      os << "nproc  "<<np<<std::endl;
      os << "nrows  "<<nrows_<<std::endl;
      os << "ncols   "<<ncols_<<std::endl;
      os << "lda    "<<lda<<std::endl;
      os << "stride "<<stride<<std::endl;
#ifdef PHIST_SDMATS_ROW_MAJOR
      os << "row-major storage"<<std::endl;
#else
      os << "col-major storage"<<std::endl;
#endif      
      }
    for (int p=0;p<np;p++)
      {
      if (p==rank)
        {
        os << " @ rank "<<p<<" @"<<std::endl;
        for (int i=0;i<stride*nrows_;i+=stride)
          {
          for (int j=0;j<ncols_;j++)
            {
            os << mat_vp[MIDX(i,j,lda)]<<"  ";
            }//j
          os << std::endl;
          }//i
        os << std::flush;
        }//rank==p
#ifdef PHIST_HAVE_MPI
        MPI_Barrier(mpi_comm);
        MPI_Barrier(mpi_comm);
        MPI_Barrier(mpi_comm);
        MPI_Barrier(mpi_comm);
        MPI_Barrier(mpi_comm);
#endif
      }//p
    return;
  }
  
  TYPE(sdMat_ptr) mat1_, mat2_, mat3_, mat4_;
  _ST_ *mat1_vp_, *mat2_vp_, *mat3_vp_, *mat4_vp_;
  static const int nrows_=_Nrows;
  static const int ncols_=_Ncols;
  lidx_t m_lda_;
  };


template<int _Nrows, int _Ncols>
const int KernelTestWithSdMats<_ST_,_Nrows,_Ncols>::nrows_;

template<int _Nrows, int _Ncols>
const int KernelTestWithSdMats<_ST_,_Nrows,_Ncols>::ncols_;
