/*! Test fixture. */
template<gidx_t _Nglob, int _Nvec, bool _useViews>
class KernelTestWithVectors<_ST_,_Nglob,_Nvec, _useViews> : 
        public virtual KernelTestWithType< _ST_ >,
        public virtual KernelTestWithMap<_Nglob> 
  {

private:

void createVecs()
{
  int pad_pre=0;
  int pad_post=0;
  // vectors created with the same function should get the same stride (lda)
  lidx_t lda;
  if (useViews_)
  {
    pad_pre=2;
    pad_post=3;
  }
  int pad=pad_pre+pad_post;
  SUBR(mvec_create)(&vec1_,this->map_,nvec_+pad,&this->iflag_);
  ASSERT_EQ(0,this->iflag_);

  SUBR(mvec_create)(&vec2_,this->map_,nvec_+pad,&this->iflag_);
  ASSERT_EQ(0,this->iflag_);

  SUBR(mvec_create)(&vec3_,this->map_,nvec_+pad,&this->iflag_);
  ASSERT_EQ(0,this->iflag_);
  
  // if requested, set vecX to views of memX
  if (useViews_)
  {
    mem1_=vec1_; vec1_=NULL;
    mem2_=vec2_; vec2_=NULL;
    mem3_=vec3_; vec3_=NULL;
    SUBR(mvec_view_block)(mem1_,&vec1_,pad_pre,pad_pre+nvec_-1,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    SUBR(mvec_view_block)(mem2_,&vec2_,pad_pre,pad_pre+nvec_-1,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    SUBR(mvec_view_block)(mem3_,&vec3_,pad_pre,pad_pre+nvec_-1,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
  }

  // extract raw data views and check that the stride (lda) is the same for all mvecs
  // with the same number of cols

  SUBR(mvec_extract_view)(vec1_,&vec1_vp_,&lda_,&this->iflag_);
  ASSERT_EQ(0,this->iflag_);
  lda=lda_;
  
  SUBR(mvec_extract_view)(vec2_,&vec2_vp_,&lda_,&this->iflag_);
  ASSERT_EQ(0,this->iflag_);
  ASSERT_EQ(lda,lda_);
  
  SUBR(mvec_extract_view)(vec3_,&vec3_vp_,&lda_,&this->iflag_);
  ASSERT_EQ(0,this->iflag_);
  ASSERT_EQ(lda,lda_);
  stride_=1;
}

void deleteVecs()
{
  SUBR(mvec_delete)(this->vec1_,&this->iflag_);
  ASSERT_EQ(0,this->iflag_);
  SUBR(mvec_delete)(this->vec2_,&this->iflag_);
  ASSERT_EQ(0,this->iflag_);
  SUBR(mvec_delete)(this->vec3_,&this->iflag_);
  ASSERT_EQ(0,this->iflag_);

  // delete memory blocks if vecX are views
  if (useViews_)
  {
    SUBR(mvec_delete)(this->mem1_,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    SUBR(mvec_delete)(this->mem2_,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
    SUBR(mvec_delete)(this->mem3_,&this->iflag_);
    ASSERT_EQ(0,this->iflag_);
  }
}

public:

  /*! Set up routine.
   */
virtual void SetUp()
  {
  KernelTestWithType< ST >::SetUp();
  mem1_=NULL; mem2_=NULL; mem3_=NULL;
  if (typeImplemented_)
    {
    lidx_t lda;
    // vectors created with the same function should get the same stride (lda)
    KernelTestWithMap<_Nglob>::SetUp();
    createVecs();
    }
  }

  /*! Clean up.
   */
virtual void TearDown() 
  {
  if (this->typeImplemented_)
    {
    deleteVecs();
    KernelTestWithMap<_Nglob>::TearDown();
    }
  KernelTestWithType< ST >::TearDown();
  }

  /*! Replace the map and rebuild vectors
   */
virtual void replaceMap(const_map_ptr_t map)
  {
    if (this->typeImplemented_)
      {
      deleteVecs();

      KernelTestWithMap<_Nglob>::replaceMap(map);

      createVecs();
      }
  }

// tolerance for tests depending on the vector length
inline static MT releps(TYPE(const_mvec_ptr) V=NULL)
  {
  if (V==NULL) return std::sqrt((MT)_Nglob*mt::eps());
  int nvec,iflag;
  SUBR(mvec_num_vectors)(V,&nvec,&iflag);
  MT *nrms = new MT[nvec];
  MT max_nrm=0;
  SUBR(mvec_norm2)(V,nrms,&iflag);
  for (int i=0;i<nvec;i++)
    {
    max_nrm = std::max(max_nrm,nrms[i]);
    }
  if (max_nrm<4*mt::eps()) max_nrm=mt::sqrt((MT)_Nglob);
  delete[] nrms;
  return max_nrm*mt::eps();
  }

//! in-place reduction operation on scalar data type (for testing with MPI)
static int global_sum(ST* value, int count, MPI_Comm mpi_comm)
  {
  int iflag=0;
#ifdef PHIST_HAVE_MPI
        ST* gvalue = new ST[count];
        for (int i=0;i<count;i++) gvalue[i]=st::zero();
        iflag=MPI_Allreduce(value,gvalue,count,
                st::mpi_type(), MPI_SUM, mpi_comm);
        for (int i=0;i<count;i++) value[i]=gvalue[i];
        delete [] gvalue;
#endif
  return iflag;
  }

//! in-place reduction operation on scalar data type (for testing with MPI)
static int global_msum(MT* value, int count, MPI_Comm mpi_comm)
  {
  int iflag=0;
#ifdef PHIST_HAVE_MPI
        ST* gvalue = new MT[count];        
        iflag=MPI_Allreduce(value,gvalue,count,
                mt::mpi_type(), MPI_SUM, mpi_comm);
        for (int i=0;i<count;i++) value[i]=gvalue[i];
        delete [] gvalue;
#endif
  return iflag;
  }

  //! tests if each column of an mv is normalized in the 2-norm
  static MT ColsAreNormalized(const ST* vec_vp, lidx_t nloc, lidx_t lda, lidx_t stride,
        MPI_Comm mpi_comm)
    {
    MT res=1.0;
    // see if all columns in vec2 have 2-norm 1
    ST *norms = new ST[nvec_];
    for (int j=0;j<nvec_;j++)
      {
      ST sum=st::zero();
      for (int i=0;i<stride*nloc;i+=stride)
        {
        ST val=vec_vp[VIDX(i,j,lda)];
        sum+=st::conj(val)*val; 
        }
      norms[j]=sum;
      }
    global_sum(norms,nvec_,mpi_comm);
    for (int j=0;j<nvec_;j++) norms[j]=st::sqrt(norms[j]);
    res=ArrayEqual(norms,nvec_,1,nvec_,1,st::one());
    delete [] norms;
    return res;
    }

  // check if vectors are mutually orthogonal after QR factorization
  static MT ColsAreOrthogonal(ST* vec_vp, lidx_t nloc, lidx_t lda, lidx_t stride,MPI_Comm mpi_comm) 
  {
    MT res=mt::one();
    int nsums=(nvec_*nvec_-nvec_)/2;
    if(nsums == 0)
      return mt::one();
    else
    {
      ST sums[nsums+1];
      int k=0;
      for (int j1=0;j1<nvec_;j1++)
        for (int j2=j1+1;j2<nvec_;j2++)
        {
          ST sum=st::zero();
          for (int i=0;i<stride*nloc;i+=stride)
          {
            ST val1=vec_vp[VIDX(i,j1,lda)];
            ST val2=vec_vp[VIDX(i,j2,lda)];
            sum+=val1*st::conj(val2);
          }
          sums[k++]=sum;
        }
      global_sum(sums,nsums,mpi_comm);
      res=ArrayEqual(sums,nsums,1,nsums,1,st::zero());
      return res;
    }
  }

  // check if vectors are mutually orthogonal after QR factorization
  static void PrintVector(std::ostream& os, std::string label, 
        ST* vec_vp, lidx_t nloc, lidx_t lda, lidx_t stride,MPI_Comm mpi_comm) 
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
      os << "nglob  "<<_Nglob<<std::endl;
      os << "nloc   "<<nloc<<std::endl;
      os << "lda    "<<lda<<std::endl;
      os << "stride "<<stride<<std::endl;
#ifdef PHIST_MVECS_ROW_MAJOR
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
        for (int i=0;i<stride*nloc;i+=stride)
          {
          for (int j=0;j<nvec_;j++)
            {
            os << vec_vp[VIDX(i,j,lda)]<<"  ";
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
    
  static bool pointerUnchanged(TYPE(mvec_ptr) V, ST* expected_location, int expected_lda)
  { 
    int iflag;
    lidx_t lda;
    ST* ptr;
    SUBR(mvec_extract_view)(V,&ptr,&lda,&iflag);
    return ( (iflag==0)&&(lda==expected_lda)&&(ptr==expected_location) );
  }

  //! if _useViews=true, these are larger memory blocks
  //! holding vecX_ as an inner view. Otherwise they are NULL.
  TYPE(mvec_ptr) mem1_, mem2_, mem3_;

  TYPE(mvec_ptr) vec1_, vec2_, vec3_;
  ST *vec1_vp_, *vec2_vp_, *vec3_vp_;
  static const int nvec_=_Nvec;
  static const int useViews_=_useViews;
  lidx_t lda_, stride_;
  };

template<gidx_t n, int nvec, bool useViews>
const int KernelTestWithVectors<_ST_,n,nvec,useViews>::nvec_;

template<gidx_t n, int nvec, bool useViews>
const int KernelTestWithVectors<_ST_,n,nvec,useViews>::useViews_;

