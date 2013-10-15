#include "phist_macros.h"

#include "phist_typedefs.h"
#include "phist_kernels.h"

#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_RCP.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_MatrixIO.hpp"

#include "BelosTpetraAdapter.hpp"
#include "BelosTsqrOrthoManager.hpp"

extern "C" {

// we implement all the four types
void SUBR(type_avail)(int* ierr)
  {
  *ierr=0;
  }

// \name Matrix input from a file

//@{

//! read a matrix from a MatrixMarket (ASCII) file
void SUBR(crsMat_read_mm)(TYPE(crsMat_ptr)* vA, const char* filename,int* ierr)
  {
  Tpetra::MatrixMarket::Reader<Traits<_ST_>::crsMat_t> reader;
  std::string fstring(filename);

  
  Teuchos::RCP<Traits<_ST_>::crsMat_t> A;

  Teuchos::ParameterList nodeParams;
  Teuchos::RCP<node_t> node = Teuchos::rcp(new node_t(nodeParams));
  Teuchos::RCP<const comm_t> comm = Teuchos::DefaultComm<int>::getComm();
  
  _TRY_CATCH_(A=reader.readSparseFile(fstring,comm,node),*ierr);
  Teuchos::Ptr<Traits<_ST_>::crsMat_t> Aptr = A.release();
  *vA = (TYPE(crsMat_ptr))(Aptr.get());
  }

//! read a matrix from a Ghost CRS (binary) file.
void SUBR(crsMat_read_bin)(TYPE(crsMat_ptr)* vA, const char* filename,int* ierr)
  {
  // TODO - not implemented (should read the binary file format defined by ghost)
  *ierr=-99;
  }

//! read a matrix from a Harwell-Boeing (HB) file
void SUBR(crsMat_read_hb)(TYPE(crsMat_ptr)* vA, const char* filename,int* ierr)
  {
  *ierr=0;
  std::string fname(filename);
  Teuchos::RCP<const comm_t> comm = Teuchos::DefaultComm<int>::getComm();
  Teuchos::ParameterList nodeParams;
  
  Teuchos::RCP<Traits<_ST_>::crsMat_t> A;
  Teuchos::RCP<node_t> node = Teuchos::rcp(new node_t(nodeParams));
  Teuchos::RCP<map_t> rowMap=Teuchos::null; // assume linear distribution for now
  Teuchos::RCP<Teuchos::ParameterList> params=Teuchos::null;
  //_TRY_CATCH_(Tpetra::Utils::readHBMatrix(fname,comm,node,A,rowMap, params),*ierr);
  _TRY_CATCH_(Tpetra::Utils::readHBMatrix(fname,comm,node,A),*ierr);
  Teuchos::Ptr<Traits<_ST_>::crsMat_t> Aptr = A.release();
  *vA = (TYPE(crsMat_ptr))(Aptr.get());
  }
//!@}

//! \name get information about the data distribution in a matrix (maps)

//!@{
//! get the row distribution of the matrix
void SUBR(crsMat_get_row_map)(TYPE(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::crsMat_t, A, vA, *ierr);
  *vmap = (const_map_ptr_t)(A->getRowMap().get());
  }

//! get column distribution of a matrix
void SUBR(crsMat_get_col_map)(TYPE(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::crsMat_t, A, vA, *ierr);
  *vmap = (const_map_ptr_t)(A->getColMap().get());
  }

//! get the map for vectors x in y=A*x
void SUBR(crsMat_get_domain_map)(TYPE(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::crsMat_t, A, vA, *ierr);
  *vmap = (const_map_ptr_t)(A->getDomainMap().get());
  }

//! get the map for vectors y in y=A*x
void SUBR(crsMat_get_range_map)(TYPE(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::crsMat_t, A, vA, *ierr);
  *vmap = (const_map_ptr_t)(A->getRangeMap().get());
  }
//@}

//! \name constructors

//@{
//! create a block-vector. The entries are stored contiguously
//! at val in column major ordering.
void SUBR(mvec_create)(TYPE(mvec_ptr)* vV, const_map_ptr_t vmap, int nvec, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const map_t, map, vmap, *ierr);
  Teuchos::RCP<const map_t> map_ptr = Teuchos::rcp(map,false);
  Traits<_ST_>::mvec_t* V = new Traits<_ST_>::mvec_t(map_ptr,nvec);
  *vV=(TYPE(mvec_ptr))(V);
  }

//! create a block-vector as view of raw data. The map tells the object
//! how many rows it should 'see' in the data (at most lda, the leading
//! dimension of the 2D array values).
void SUBR(mvec_create_view)(TYPE(mvec_ptr)* vV, const_map_ptr_t vmap, 
        _ST_* values, lidx_t lda, int nvec,
        int* ierr)
  {
  _CAST_PTR_FROM_VOID_(const map_t, map, vmap, *ierr);
  Teuchos::RCP<const map_t> map_ptr = Teuchos::rcp(map,false);
  Teuchos::ArrayView<_ST_> val_ptr(values,lda*nvec);
  Traits<_ST_>::mvec_t* V = new Traits<_ST_>::mvec_t(map_ptr,val_ptr,lda,nvec);
  *vV=(TYPE(mvec_ptr))(V);  
  }

//! create a serial dense n x m matrix on all procs, with column major
//! ordering.
void SUBR(sdMat_create)(TYPE(sdMat_ptr)* vM, int nrows, int ncols, 
        const_comm_ptr_t vcomm, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const comm_t, comm, vcomm, *ierr);

  //TODO - add node arg

  // create local map
  Teuchos::RCP<map_t> localMap =
        Teuchos::rcp(new map_t(nrows, 0, Teuchos::rcp(comm,false), Tpetra::LocallyReplicated));
  Traits<_ST_>::sdMat_t* M = new Traits<_ST_>::mvec_t(localMap,ncols);
  *vM=(TYPE(sdMat_ptr))(M);
  }

//@}

//! retrieve local length of the vectors in V
void SUBR(mvec_my_length)(TYPE(const_mvec_ptr) vV, lidx_t* len, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,V,vV,*ierr);
  *len = V->getLocalLength();
  }

//! retrieve the map of the vectors in V
void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) vV, const_map_ptr_t* vmap, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,V,vV,*ierr);
  *vmap=(const_map_ptr_t)(V->getMap().get());
  }

//! retrieve the comm used for MPI communication in V
void SUBR(mvec_get_comm)(TYPE(const_mvec_ptr) vV, const_comm_ptr_t* vcomm, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,V,vV,*ierr);
  *vcomm=(const_comm_ptr_t)(V->getMap()->getComm().get());
  }


//! retrieve number of vectors/columns in V
void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) vV, int* nvec, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,V,vV,*ierr);
  *nvec = V->getNumVectors();
  }

//! get number of cols in local dense matrix
void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) vM, int* nrows, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::sdMat_t,M,vM,*ierr);
  *nrows = M->getLocalLength();
  }

//! get number of cols in local dense matrix
void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) vM, int* ncols, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::sdMat_t,M,vM,*ierr);
  *ncols = M->getNumVectors();
  }


//! extract view from multi-vector
void SUBR(mvec_extract_view)(TYPE(mvec_ptr) vV, _ST_** val, lidx_t* lda, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,V,vV,*ierr);
  Teuchos::ArrayRCP<_ST_> val_ptr = V->get1dViewNonConst();
  *val = val_ptr.getRawPtr();
  *lda = V->getStride();
  }

//! extract view from serial dense matrix
void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) vM, _ST_** val, lidx_t* lda, int* ierr)
  {
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::sdMat_t,M,vM,*ierr);
  Teuchos::ArrayRCP<_ST_> valptr = M->get1dViewNonConst();
  *val = valptr.getRawPtr();
  *lda=M->getStride();
  *ierr=0; //TODO - how do we get LDA?
  }

//! get a new vector that is a view of some columns of the original one,
//! Vblock = V(:,jmin:jmax). The new object Vblock is created but does not
//! allocate memory for the vector entries, instead using the entries from V
//! directly. When mvec_delete(Vblock) is called, the library has to take care
//! that the value array is not deleted 
void SUBR(mvec_view_block)(TYPE(mvec_ptr) vV,
                             TYPE(mvec_ptr)* vVblock,
                             int jmin, int jmax, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,V,vV,*ierr);
  Teuchos::RCP<Traits<_ST_>::mvec_t> Vblock;
  _TRY_CATCH_(Vblock = V->subViewNonConst(Teuchos::Range1D(jmin,jmax)),*ierr);
  if (*vVblock!=NULL)
    {
    _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,tmp,*vVblock,*ierr);
    delete [] tmp;
    }
  *vVblock = (TYPE(mvec_ptr))(Vblock.release().get());                        
  }

//! get a new vector that is a copy of some columns of the original one,  
//! Vblock = V(:,jmin:jmax). The object Vblock must be created beforehand 
//! and the corresponding columns of V are copied into the value array    
//! of Vblock. V is not modified.
void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) vV,
                             TYPE(mvec_ptr) vVblock,
                             int jmin, int jmax, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,Vblock,vVblock,*ierr);
  // get a view of the columns of V first
  Teuchos::RCP<const Traits<_ST_>::mvec_t> Vcols;
  _TRY_CATCH_(Vcols = V->subView(Teuchos::Range1D(jmin,jmax)),*ierr);
  *Vblock = *Vcols; // copy operation
  }

//! given a multi-vector Vblock, set V(:,jmin:jmax)=Vblock by copying the corresponding
//! vectors. Vblock is not modified.
void SUBR(mvec_set_block)(TYPE(mvec_ptr) vV,
                             TYPE(const_mvec_ptr) vVblock,
                             int jmin, int jmax, int* ierr)
  {
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,Vblock,vVblock,*ierr);
  // get a view of the columns of V first
  Teuchos::RCP<Traits<_ST_>::mvec_t> Vcols;
  _TRY_CATCH_(Vcols = V->subViewNonConst(Teuchos::Range1D(jmin,jmax)),*ierr);
  // copy operation
  _TRY_CATCH_(*Vcols = *Vblock, *ierr);
  }

//! get a new matrix that is a view of some rows and columns of the original one,
//! Mblock = M(imin:imax,jmin:jmax). 
void SUBR(sdMat_view_block)(TYPE(sdMat_ptr) vM,
                             TYPE(sdMat_ptr)* vMblock,
                             int imin, int imax, int jmin, int jmax, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::sdMat_t,M,vM,*ierr);
  Traits<_ST_>::sdMat_t* Mblock;
  int lda=M->getStride();
  int nrows=imax-imin+1;
  int ncols=jmax-jmin+1;
  Teuchos::ArrayRCP<_ST_> valptr0 = M->get1dViewNonConst();

  _ST_* val = valptr0.getRawPtr();
  int len=lda*nrows - imin;
  Teuchos::ArrayView<_ST_> valptr(val+imin,len);

  Teuchos::SerialComm<int> scomm;
  Teuchos::RCP<map_t> smap =
        Teuchos::rcp(new map_t(nrows, 0, Teuchos::rcpFromRef< const Teuchos::Comm<int> >(scomm), 
                Tpetra::LocallyReplicated));

  _TRY_CATCH_(Mblock = new Traits<_ST_>::sdMat_t(smap,valptr,lda,ncols),*ierr);

  if (*vMblock!=NULL)
    {
    _CAST_PTR_FROM_VOID_(Traits<_ST_>::sdMat_t,tmp,*vMblock,*ierr);
    delete [] tmp;
    }


  *vMblock = (TYPE(sdMat_ptr))Mblock;                        
  }

//! get a new matrix that is a copy of some rows and columns of the original one,  
//! Mblock = M(imin:imax,jmin:jmax). The object Mblock must be created beforehand 
//! and the corresponding columns of M are copied into the value array    
//! of Mblock. M is not modified.
void SUBR(sdMat_get_block)(TYPE(const_sdMat_ptr) M, 
                             TYPE(sdMat_ptr) Mblock,
                             int imin, int imax, int jmin, int jmax, int* ierr)
  {
  *ierr=-99;
  }

//! given a serial dense matrix Mblock, set M(imin:imax,jmin:jmax)=Mblock by 
//! copying the corresponding elements. Mblock is not modified.
void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) M, 
                             TYPE(const_sdMat_ptr) Mblock,
                             int imin, int imax, int jmin, int jmax, int* ierr)
  {
  *ierr=-99;
  }


//! \name destructors

//@{

//!
void SUBR(crsMat_delete)(TYPE(crsMat_ptr) vA, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::crsMat_t,A,vA,*ierr);
  delete A;
  }

//!
void SUBR(mvec_delete)(TYPE(mvec_ptr) vV, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,V,vV,*ierr);
  delete V;
  }

//!
void SUBR(sdMat_delete)(TYPE(sdMat_ptr) vM, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,M,vM,*ierr);
  delete M;
  }

//@}

//! \name Numerical functions
//!@{

//! put scalar value into all elements of a multi-vector
void SUBR(mvec_put_value)(TYPE(mvec_ptr) vV, _ST_ value, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,V,vV,*ierr);
  _TRY_CATCH_(V->putScalar(value),*ierr);
  }

//! put scalar value into all elements of a multi-vector
void SUBR(sdMat_put_value)(TYPE(sdMat_ptr) vM, _ST_ value, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::sdMat_t,M,vM,*ierr);
  _TRY_CATCH_(M->putScalar(value),*ierr);
  }

//! put random numbers into all elements of a multi-vector
void SUBR(mvec_random)(TYPE(mvec_ptr) vV, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,V,vV,*ierr);
  _TRY_CATCH_(V->randomize(),*ierr);
  }

void SUBR(mvec_print)(TYPE(const_mvec_ptr) vV, int* ierr)
  {
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,V,vV,*ierr);
  std::cout << *V;
  }

void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) vM, int* ierr)
  {
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::sdMat_t,M,vM,*ierr);
  std::cout << *M;
  }


//! put random numbers into all elements of a serial dense matrix
void SUBR(sdMat_random)(TYPE(sdMat_ptr) vM, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,M,vM,*ierr);
  _TRY_CATCH_(M->randomize(),*ierr);
  }

  //! compute the 2-norm) of each column of v                   
  //! (vnrm[i] must be pre-allocated by caller)
  void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) vV,
                            _MT_* vnrm, int* ierr) 
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,V,vV,*ierr);
int nvec = V->getNumVectors();
  Teuchos::ArrayView<_MT_> norms(vnrm,nvec);
  _TRY_CATCH_(V->norm2(norms),*ierr);
  return;
  }


  //! normalize (in the 2-norm) each column of v and return ||v||_2
  //! for each vector i in vnrm[i] (must be pre-allocated by caller)
  void SUBR(mvec_normalize)(TYPE(mvec_ptr) vV,
                            _MT_* vnrm, int* ierr) 
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,V,vV,*ierr);
int nvec = V->getNumVectors();
  Teuchos::ArrayView<_MT_> norms(vnrm,nvec);
  Teuchos::Array<_ST_> scaling(nvec);
  _TRY_CATCH_(V->norm2(norms),*ierr);
  for (int i=0;i<nvec;i++)
    {
    scaling[i]=1.0/norms[i];
    }
  _TRY_CATCH_(V->scale(scaling),*ierr);
  return;
  }

//! scale each column i of v and by scalar
void SUBR(mvec_scale)(TYPE(mvec_ptr) vV, 
                            _ST_ scalar, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,V,vV,*ierr);
  _TRY_CATCH_(V->scale(scalar),*ierr);
  return;
  }

//! scale each column i of v and by scalar[i]
void SUBR(mvec_vscale)(TYPE(mvec_ptr) vV, 
                            _ST_* scalar, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,V,vV,*ierr);
int nvec = V->getNumVectors();
  Teuchos::ArrayView<_ST_> scal(scalar,nvec);
  _TRY_CATCH_(V->scale(scal),*ierr);
  return;
  }

//! y=alpha*x+beta*y
void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) vX,
                            _ST_ beta,  TYPE(mvec_ptr)       vY, 
                            int* ierr)
  {
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,X,vX,*ierr);
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,Y,vY,*ierr);
  _TRY_CATCH_(Y->update(alpha,*X,beta),*ierr);
  }

//! y[i]=alpha[i]*x[i]+beta*y[i], i=1..nvec
void SUBR(mvec_vadd_mvec)(const _ST_ alpha[], TYPE(const_mvec_ptr) vX,
                            _ST_ beta,  TYPE(mvec_ptr)       vY, 
                            int* ierr)
  {
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,X,vX,*ierr);
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,Y,vY,*ierr);
  
  for (int i=0;i<X->getNumVectors(); i++)
    {
    _TRY_CATCH_(Y->getVectorNonConst(i)->update(alpha[i],*X->getVector(i), beta),*ierr);
    }
  }

//! B=alpha*A+beta*B
void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA,
                            _ST_ beta,  TYPE(sdMat_ptr)       vB, 
                            int* ierr)
  {
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::sdMat_t,A,vA,*ierr);
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::sdMat_t,B,vB,*ierr);
  _TRY_CATCH_(B->update(alpha,*A,beta),*ierr);
  }


//! y=alpha*A*x+beta*y.
void SUBR(crsMat_times_mvec)(_ST_ alpha, TYPE(const_crsMat_ptr) vA, 
                                        TYPE(const_mvec_ptr) vx, 
                                        _ST_ beta, TYPE(mvec_ptr) vy, 
                                        int* ierr)
  {
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::crsMat_t,A,vA,*ierr);
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,x,vx,*ierr);
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,y,vy,*ierr);
  PHIST_OUT(9,"alpha=%g+i%g\n",st::real(alpha),st::imag(alpha));
  PHIST_OUT(9,"beta=%g+i%g\n",st::real(beta),st::imag(beta));
  Traits<_ST_>::crsMVM_t spMVM(Teuchos::rcp(A,false));
  _TRY_CATCH_(spMVM.apply(*x,*y,Teuchos::NO_TRANS,alpha,beta),*ierr);

  }

//! dot product of vectors v_i and w_i, i=1..numvecs
void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) vv, TYPE(const_mvec_ptr) vw, _ST_* s, 
int* 
ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,v,vv,*ierr);
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,w,vw,*ierr);
  Teuchos::ArrayView<_ST_> dots(s,v->getNumVectors());
  _TRY_CATCH_(v->dot(*w,dots),*ierr);
  }

//! dense tall skinny matrix-matrix product yielding a serial dense matrix
//! C=alpha*V'*W+beta*C. C is replicated on all MPI processes sharing V and W.
void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) vV, 
                           TYPE(const_mvec_ptr) vW, _ST_ beta, 
                           TYPE(sdMat_ptr) vC, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,W,vW,*ierr);
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::sdMat_t,C,vC,*ierr);
#ifdef _IS_COMPLEX_
  _TRY_CATCH_(C->multiply(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS,alpha,*V,*W,beta),*ierr);
#else
  _TRY_CATCH_(C->multiply(Teuchos::TRANS, Teuchos::NO_TRANS,alpha,*V,*W,beta),*ierr);
#endif
  }

//! n x m multi-vector times m x m dense matrix gives n x m multi-vector,
//! W=alpha*V*C + beta*W
void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) vV,
                                       TYPE(const_sdMat_ptr) vC,
                           _ST_ beta,  TYPE(mvec_ptr) vW,
                                       int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::mvec_t,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::sdMat_t,C,vC,*ierr);
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,W,vW,*ierr);
  _TRY_CATCH_(W->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,
        alpha, *V, *C, beta), *ierr);
  }

//! n x m serial dense matrix times m x k serial dense matrix gives n x k multi-vector,
//! C=alpha*V*W + beta*C
void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vV,
                                           TYPE(const_sdMat_ptr) vW,
                               _ST_ beta, TYPE(sdMat_ptr) vC,
                                       int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::sdMat_t,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(const Traits<_ST_>::sdMat_t,W,vW,*ierr);
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::sdMat_t,C,vC,*ierr);
  _TRY_CATCH_(C->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,
        alpha, *V, *W, beta), *ierr);
  }


//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.   
//! Q is computed in place of V. If V does not have full rank, ierr>0   
//! indicates the dimension of the null-space of V. The first m-ierr    
//! columns of Q are an orthogonal basis of the column space of V, the  
//! remaining columns form a basis for the null space.
void SUBR(mvec_QR)(TYPE(mvec_ptr) vV, TYPE(sdMat_ptr) vR, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::mvec_t,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(Traits<_ST_>::sdMat_t,R,vR,*ierr);
  
  if (R->isConstantStride()==false)
    {
    *ierr = -1;
    return;
    }

  int stride = R->getStride();
  int nrows = R->getLocalLength();
  int ncols = R->getNumVectors();
    
#ifdef TESTING
  _CHECK_ZERO_(nrows-ncols,*ierr);
  _CHECK_ZERO_(nrows-V->getNumVectors(),*ierr);
#endif  

  Teuchos::RCP<Traits< _ST_ >::Teuchos_sdMat_t> R_view
        = Traits< _ST_ >::CreateTeuchosViewNonConst(Teuchos::rcp(R,false),ierr);
  if (*ierr) return;      
  
  Belos::TsqrOrthoManager< _ST_ , Traits< _ST_ >::mvec_t> tsqr("phist");
  Teuchos::RCP<const Teuchos::ParameterList> valid_params = 
        tsqr.getValidParameters();
  // faster but numerically less robust settings:
  Teuchos::RCP<const Teuchos::ParameterList> fast_params = 
        tsqr.getFastParameters();
  Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp
        (new Teuchos::ParameterList(*valid_params));
  params->set("randomizeNullSpace",true);
  tsqr.setParameterList(params);

  int rank;
  _TRY_CATCH_(rank = tsqr.normalize(*V,R_view),*ierr);  
  *ierr = ncols-rank;// return positive number if rank not full.
  return;
  }

//!@}

} // extern "C"

