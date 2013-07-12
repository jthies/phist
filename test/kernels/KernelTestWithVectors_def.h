/*! Test fixture. */
template<int _Nglob, int _Nvec>
class KernelTestWithVectors<_ST_,_Nglob,_Nvec> : 
        public virtual KernelTestWithType< _ST_ >,
        public virtual KernelTestWithMap<_Nglob> 
  {

public:

  /*! Set up routine.
   */
virtual void SetUp()
  {
  KernelTestWithType< ST >::SetUp();
  if (this->typeImplemented_)
    {
    KernelTestWithMap<_Nglob>::SetUp();
    _SUBR_(mvec_create)(&vec1_,this->map_,this->nvec_,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);
    _SUBR_(mvec_extract_view)(vec1_,&vec1_vp_,&this->lda_,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);
    _SUBR_(mvec_create)(&vec2_,this->map_,this->nvec_,&this->ierr_);
    ASSERT_EQ(0,this->ierr_);
    _SUBR_(mvec_extract_view)(vec2_,&vec2_vp_,&this->lda_,&this->ierr_);
        ASSERT_EQ(0,this->ierr_);
    stride_=1;
    }
  }

  /*! Clean up.
   */
virtual void TearDown() 
  {
  if (this->typeImplemented_)
    {
    _SUBR_(mvec_delete)(vec1_,&this->ierr_);
    _SUBR_(mvec_delete)(vec2_,&this->ierr_);
    KernelTestWithMap<_Nglob>::TearDown();
    }
  KernelTestWithType< ST >::TearDown();
  }

//! in-place reduction operation on scalar data type (for testing with MPI)
static int global_sum(ST* value, int count, MPI_Comm mpi_comm)
  {
  int ierr=0;
#ifdef PHIST_HAVE_MPI
        ST* gvalue = new ST[count];
        for (int i=0;i<count;i++) gvalue[i]=st::zero();
        ierr=MPI_Allreduce(value,gvalue,count,
                st::mpi_type(), MPI_SUM, mpi_comm);
        for (int i=0;i<count;i++) value[i]=gvalue[i];
        delete [] gvalue;
#endif
  return ierr;
  }

//! in-place reduction operation on scalar data type (for testing with MPI)
static int global_msum(MT* value, int count, MPI_Comm mpi_comm)
  {
  int ierr=0;
#ifdef PHIST_HAVE_MPI
        ST* gvalue = new MT[count];        
        ierr=MPI_Allreduce(value,gvalue,count,
                mt::mpi_type(), MPI_SUM, mpi_comm);
        for (int i=0;i<count;i++) value[i]=gvalue[i];
        delete [] gvalue;
#endif
  return ierr;
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
        ST val=vec_vp[j*lda+i];
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
      ST sums[nsums];
      int k=0;
      for (int j1=0;j1<nvec_;j1++)
      for (int j2=j1+1;j2<nvec_;j2++)
        {
        ST sum=st::zero();
        for (int i=0;i<stride*nloc;i+=stride)
          {
          ST val1=vec_vp[j1*lda+i];
          ST val2=vec_vp[j2*lda+i];
          sum+=val1*st::conj(val2);
          }
        sums[k++]=sum;
        }
      global_sum(sums,nsums,mpi_comm);
      res=ArrayEqual(sums,nsums,1,nsums,1,st::zero());
      return res;
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
            os << vec_vp[j*lda+i]<<"  ";
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


  _TYPE_(mvec_ptr) vec1_, vec2_;
  ST *vec1_vp_, *vec2_vp_;
  static const int nvec_=_Nvec;
  lidx_t lda_, stride_;
  };

template<int n, int nvec>
const int KernelTestWithVectors<_ST_,n,nvec>::nvec_;

