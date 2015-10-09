
/*! Test fixture. */
template<int _Nrows, int _Ncols, int _useViews>
class KernelTestWithSdMats<_ST_,_Nrows,_Ncols, _useViews> : 
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
  if (typeImplemented_)
  {
    bool align_p2=(useViews_==2);
    int pad_nr_pre=0, pad_nc_pre=0;
    int pad_nr_post=0, pad_nc_post=0;
    int nr_padded=nrows_, nc_padded=ncols_;
    if (useViews_)
    {
      pad_nr_pre=align_p2?8:3;
      pad_nr_post=1;
      pad_nc_pre=align_p2?8:5;
      pad_nc_post=1;
    // padding to power of 2
    int nr_pow_p2=(int)ceil(log((double)(nrows_+pad_nr_pre+pad_nr_post))/log(2.0))+1;
    int nc_pow_p2=(int)ceil(log((double)(ncols_+pad_nc_pre+pad_nc_post))/log(2.0))+1;
    int nr_padded_p2=(int)pow(2.0,(double)nr_pow_p2);
    int nc_padded_p2=(int)pow(2.0,(double)nc_pow_p2);
    nr_padded=align_p2?nr_padded_p2:nrows_+pad_nr_pre+pad_nr_post;
    nc_padded=align_p2?nc_padded_p2:ncols_+pad_nc_pre+pad_nc_post;
  }

    SUBR(sdMat_create)(&mem1_,nr_padded,nc_padded,this->comm_,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    SUBR(sdMat_create)(&mem2_,nr_padded,nc_padded,this->comm_,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    SUBR(sdMat_create)(&mem3_,nr_padded,nc_padded,this->comm_,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    SUBR(sdMat_create)(&mem4_,nr_padded,nc_padded,this->comm_,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    
    if (useViews_)
    {
      mat1_=NULL; mat2_=NULL; mat3_=NULL; mat4_=NULL;
      SUBR(sdMat_view_block)(mem1_,&mat1_,pad_nr_pre,pad_nr_pre+nrows_-1,
                             pad_nc_pre,pad_nc_pre+ncols_-1,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_view_block)(mem2_,&mat2_,pad_nr_pre,pad_nr_pre+nrows_-1,
                             pad_nc_pre,pad_nc_pre+ncols_-1,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_view_block)(mem3_,&mat3_,pad_nr_pre,pad_nr_pre+nrows_-1,
                             pad_nc_pre,pad_nc_pre+ncols_-1,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
      SUBR(sdMat_view_block)(mem4_,&mat4_,pad_nr_pre,pad_nr_pre+nrows_-1,
                             pad_nc_pre,pad_nc_pre+ncols_-1,&this->iflag_);
      ASSERT_EQ(0,this->iflag_);
    }
    else
    {
      mat1_=mem1_;
      mat2_=mem2_;
      mat3_=mem3_;
      mat4_=mem4_;
    }
    
    SUBR(sdMat_extract_view)(mat1_,&mat1_vp_,&this->m_lda_,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    SUBR(sdMat_extract_view)(mat2_,&mat2_vp_,&this->m_lda_,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    SUBR(sdMat_extract_view)(mat3_,&mat3_vp_,&this->m_lda_,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    SUBR(sdMat_extract_view)(mat4_,&mat4_vp_,&this->m_lda_,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
  }
}

  /*! Clean up.
   */
virtual void TearDown() 
{
  if (typeImplemented_)
  {
    if (useViews_)
    {
      SUBR(sdMat_delete)(mat1_,&this->iflag_);
      SUBR(sdMat_delete)(mat2_,&this->iflag_);
      SUBR(sdMat_delete)(mat3_,&this->iflag_);
      SUBR(sdMat_delete)(mat4_,&this->iflag_);    
    }
    SUBR(sdMat_delete)(mem4_,&this->iflag_);
    SUBR(sdMat_delete)(mem3_,&this->iflag_);
    SUBR(sdMat_delete)(mem2_,&this->iflag_);
    SUBR(sdMat_delete)(mem1_,&this->iflag_);
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

  static bool pointerUnchanged(TYPE(sdMat_ptr) V, ST* expected_location, int expected_lda)
  {
    int iflag;
    lidx_t lda;
    ST* ptr;
    SUBR(sdMat_extract_view)(V,&ptr,&lda,&iflag);
    return ( (iflag==0)&&(lda==expected_lda)&&(ptr==expected_location) );
  }

  
  TYPE(sdMat_ptr) mem1_, mem2_, mem3_, mem4_;
  TYPE(sdMat_ptr) mat1_, mat2_, mat3_, mat4_;
  _ST_ *mat1_vp_, *mat2_vp_, *mat3_vp_, *mat4_vp_;
  static const int nrows_=_Nrows;
  static const int ncols_=_Ncols;
  static const int useViews_=_useViews;
  lidx_t m_lda_;
};


template<int _Nrows, int _Ncols, int _useViews>
const int KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews>::nrows_;

template<int _Nrows, int _Ncols, int _useViews>
const int KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews>::ncols_;

template<int _Nrows, int _Ncols, int _useViews>
const int KernelTestWithSdMats<_ST_,_Nrows,_Ncols,_useViews>::useViews_;

