extern "C" {

// we implement only the double precision real type D
void SUBR(type_avail)(int* ierr)
  {
  *ierr=0;
  }

// \name Matrix input from a file
//@{


//! read a matrix from a MatrixMarket (ASCII) file
void SUBR(crsMat_read_mm)(TYPE(crsMat_ptr)* vA, const char* filename,int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr = -99; // not implemented in ghost, use converter script to bin crs
  }

//! read a matrix from a Ghost CRS (binary) file.
void SUBR(crsMat_read_bin)(TYPE(crsMat_ptr)* vA, const char* filename,int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"

  ghost_mat_t* mat;
  ghost_context_t *ctx;

  ghost_mtraits_t *mtraits=new ghost_mtraits_t;
        mtraits->format = GHOST_SPM_FORMAT_CRS;
        mtraits->datatype = st::ghost_dt;

  ctx = ghost_createContext(GHOST_GET_DIM_FROM_MATRIX,GHOST_GET_DIM_FROM_MATRIX,
        GHOST_CONTEXT_DEFAULT,(char*)filename,MPI_COMM_WORLD,1.0);
  mat = ghost_createMatrix(ctx,mtraits,1);                               
  mat->fromFile(mat,const_cast<char*>(filename));
#if DEBUG > 0  
  ghost_printContextInfo(ctx);
  ghost_printMatrixInfo(mat);
#endif
  *vA = (TYPE(crsMat_ptr))mat;  
  }

//! read a matrix from a Harwell-Boeing (HB) file
void SUBR(crsMat_read_hb)(TYPE(crsMat_ptr)* vA, const char* filename,int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr = -99; // not implemented in ghost, use converter script to bin crs
  }

//!@}

//! \name get information about the data distribution in a matrix (maps)

//!@{
//! get the row distribution of the matrix
void SUBR(crsMat_get_row_map)(TYPE(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const ghost_mat_t,A,vA,*ierr);
  ghost_map_t* map = new ghost_map_t;
  map->ctx = A->context;
  map->vtraits_template=phist_default_vtraits();
  map->vtraits_template->flags=GHOST_VEC_LHS;
  *vmap = (const_map_ptr_t)map;
  }

//! get column distribution of a matrix
//! we currently treat all maps as the same as we don't allow any fancy
//! operations using them anyway and ghost can handle both halo'd (colmap)
//! and standard (rowmap) vectors in the mvm.
void SUBR(crsMat_get_col_map)(TYPE(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const ghost_mat_t,A,vA,*ierr);
  ghost_map_t* map = new ghost_map_t;
  map->ctx = A->context;
  map->vtraits_template=phist_default_vtraits();
  map->vtraits_template->flags=GHOST_VEC_RHS;
  *vmap = (const_map_ptr_t)map;
  }

//! get the map for vectors x in y=A*x
//! we currently treat all maps as the same as we don't allow any fancy
//! operations using them anyway and ghost can handle both halo'd (colmap)
//! and standard (rowmap) vectors in the mvm.
void SUBR(crsMat_get_domain_map)(TYPE(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  SUBR(crsMat_get_col_map)(vA,vmap,ierr);
  }

//! get the map for vectors y in y=A*x
//! we currently treat all maps as the same as we don't allow any fancy
//! operations using them anyway and ghost can handle both halo'd (colmap)
//! and standard (rowmap) vectors in the mvm.
void SUBR(crsMat_get_range_map)(TYPE(const_crsMat_ptr) vA, const_map_ptr_t* vmap, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  SUBR(crsMat_get_row_map)(vA,vmap,ierr);
  }
//@}

//! \name constructors

//@{
//! create a block-vector. The entries are stored contiguously
//! at val in column major ordering.
void SUBR(mvec_create)(TYPE(mvec_ptr)* vV, 
        const_map_ptr_t vmap, int nvec, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const ghost_map_t, map,vmap,*ierr);
  ghost_vec_t* result;
  ghost_vtraits_t *vtraits = ghost_cloneVtraits(map->vtraits_template);
        vtraits->nvecs=nvec;
        vtraits->datatype = st::ghost_dt;
  result=ghost_createVector(map->ctx,vtraits);
  ST zero = st::zero();
  // this allocates the vector and fills it with zeros
  result->fromScalar(result,&zero);
  PHIST_OUT(9,"mvec nrows: %ld\n",result->traits->nrows);
  *vV=(TYPE(mvec_ptr))(result);
  }

//! create a block-vector as view of raw data. The map tells the object
//! how many rows it should 'see' in the data (at most lda, the leading
//! dimension of the 2D array values).
void SUBR(mvec_create_view)(TYPE(mvec_ptr)* vV, const_map_ptr_t vmap, 
        _ST_* values, lidx_t lda, int nvec,
        int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=-99;
/*
  _CAST_PTR_FROM_VOID_(const map_t, map, vmap, *ierr);
  Teuchos::RCP<const map_t> map_ptr = Teuchos::rcp(map,false);
  Teuchos::ArrayView<_ST_> val_ptr(values,lda*nvec);
  Traits<_ST_>::mvec_t* V = new Traits<_ST_>::mvec_t(map_ptr,val_ptr,lda,nvec);
  *vV=(TYPE(mvec_ptr))(V);  
*/  
}


//! create a serial dense n x m matrix on all procs, with column major
//! ordering.
void SUBR(sdMat_create)(TYPE(sdMat_ptr)* vM, int nrows, int ncols, 
        const_comm_ptr_t vcomm, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const MPI_Comm, comm,vcomm,*ierr);
  ghost_vec_t* result;
  ghost_vtraits_t *dmtraits=new ghost_vtraits_t;
        dmtraits->flags = GHOST_VEC_DEFAULT; 
        dmtraits->nrows=nrows;
        dmtraits->nrowshalo=nrows;
        dmtraits->nrowspadded=ghost_pad(nrows,GHOST_PAD_MAX);
        dmtraits->nvecs=ncols;
        dmtraits->datatype=st::ghost_dt;

#ifdef TESTING
PHIST_CHK_IERR(comm==NULL,*ierr);
if (*comm==MPI_COMM_WORLD)
  {
  PHIST_OUT(2, "create sdMat with MPI_COMM_WORLD");
  }
else if (*comm==MPI_COMM_SELF)
  {
  PHIST_OUT(2, "create sdMat with MPI_COMM_SELF");
  }
else
  {
  PHIST_OUT(2, "create sdMat with non-standard comm");
  }
#endif

  // I think the sdMat should not have a context
  ghost_context_t* ctx=NULL;
/*
  ghost_context_t* ctx=ghost_createContext(nrows, ncols, GHOST_CONTEXT_DEFAULT, 
        NULL, *comm, 1.0);
*/
  result=ghost_createVector(ctx,dmtraits);
  ST zero = st::zero();
  result->fromScalar(result,&zero);
  *vM=(TYPE(sdMat_ptr))result;
  }

//@}

//! retrieve local length of the vectors in V
void SUBR(mvec_my_length)(TYPE(const_mvec_ptr) vV, lidx_t* len, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;
  _CAST_PTR_FROM_VOID_(const ghost_vec_t,V,vV,*ierr);
  PHIST_OUT(9,"vV @ %p",vV);
  PHIST_OUT(9,"V @ %p",V);
  *len = V->traits->nrows;
  }

//! retrieve the map of the vectors in V
void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) vV, const_map_ptr_t* vmap, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const ghost_vec_t,V,vV,*ierr);
  ghost_map_t* map = new ghost_map_t;
  map->ctx=V->context; 
  *vmap=(const_map_ptr_t)map;
  }

//! retrieve the comm used for MPI communication in V
void SUBR(mvec_get_comm)(TYPE(const_mvec_ptr) vV, const_comm_ptr_t* vcomm, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const ghost_vec_t,V,vV,*ierr);

  if (V->context!=NULL)
    {
    *vcomm=(const_comm_ptr_t)V->context->mpicomm;
    }
  else
    {
    MPI_Comm* comm = new MPI_Comm;
    *comm=MPI_COMM_SELF;
    *vcomm=(const_comm_ptr_t)(comm);
    }
  return;
  }

//! retrieve number of vectors/columns in V
void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) vV, int* nvec, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr = 0;
  _CAST_PTR_FROM_VOID_(const ghost_vec_t,V,vV,*ierr);
  *nvec = V->traits->nvecs;
  }

//! get number of cols in local dense matrix
void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) vM, int* nrows, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const ghost_vec_t,M,vM,*ierr);
  *nrows = M->traits->nrows;
  }
  
//! get number of cols in local dense matrix
void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) vM, int* ncols, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(const ghost_vec_t,M,vM,*ierr);
  *ncols = M->traits->nvecs;
  }


void SUBR(mvec_extract_view)(TYPE(mvec_ptr) vV, _ST_** val, lidx_t* lda, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  _CAST_PTR_FROM_VOID_(ghost_vec_t,V, vV, *ierr);
  *val = (_ST_*)V->val;
  *lda = V->traits->nrowspadded;
  }

void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) vM, _ST_** val, lidx_t* lda, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  _CAST_PTR_FROM_VOID_(ghost_vec_t,M, vM, *ierr);
  *val = (_ST_*)M->val;
  *lda = M->traits->nrowspadded;
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
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,V,vV,*ierr);
  // TODO - should we pass in nrows, nrowshalo for nr
  ghost_vec_t *Vblock = V->viewVec(V, jmax-jmin+1, 0);
  if (vVblock!=NULL)
    {
    _CAST_PTR_FROM_VOID_(ghost_vec_t,tmp,vVblock,*ierr);
    tmp->destroy(tmp);
    }
  *vVblock = (TYPE(mvec_ptr))Vblock;
  }

//! get a new vector that is a copy of some columns of the original one,  
//! Vblock = V(:,jmin:jmax). The object Vblock must be created beforehand 
//! and the corresponding columns of V are copied into the value array    
//! of Vblock. V is not modified.
void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) vV,
                             TYPE(mvec_ptr) vVblock,
                             int jmin, int jmax, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(ghost_vec_t,Vblock,vVblock,*ierr);
  *ierr=0;
  //TODO check bounds of Vblock
  Vblock->fromVec(Vblock,V,jmin);
  }

//! given a multi-vector Vblock, set V(:,jmin:jmax)=Vblock by copying the corresponding
//! vectors. Vblock is not modified.
void SUBR(mvec_set_block)(TYPE(mvec_ptr) vV,
                             TYPE(const_mvec_ptr) vVblock,
                             int jmin, int jmax, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(ghost_vec_t,Vblock,vVblock,*ierr);
  // TODO - bounds checking
  // create a view of the requested columns of V
  ghost_vec_t *Vcols = Vcols->viewVec(V,jmax-jmin+1,jmin);
  // copy the data
  Vcols->fromVec(Vcols,Vblock,0);
  // delete the view
  Vcols->destroy(Vcols);
  }

//! get a new sdMat that is a view of some rows and columns of the original one,
//! Mblock = M(imin:imax,jmin:jmax). The new object Vblock is created but does not
//! allocate memory for the vector entries, instead using the entries from V
//! directly.
void SUBR(sdMat_view_block)(TYPE(mvec_ptr) vM, TYPE(mvec_ptr)* vMblock,
                             int imin, int imax, int jmin, int jmax, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,M,vM,*ierr);
  // first just create a view of the corresponding columns
  ghost_vec_t *Mblock = M->viewVec(M, jmax-jmin+1, 0);
  // adjust the offset and the number of rows seen by the object
  Mblock->traits->nrows=imax-imin+1;
  std::ptrdiff_t offset=(std::ptrdiff_t)imin*(std::ptrdiff_t)ghost_sizeofDataType(Mblock->traits->datatype);
  Mblock->val=(void*)((char*)(Mblock->val)+offset);
  if (vMblock!=NULL)
    {
    _CAST_PTR_FROM_VOID_(ghost_vec_t,tmp,vMblock,*ierr);
    tmp->destroy(tmp);
    }
  *vMblock = (TYPE(mvec_ptr))Mblock;
  }

//! get a new matrix that is a copy of some rows and columns of the original one,  
//! Mblock = M(imin:imax,jmin:jmax). The object Mblock must be created beforehand 
//! and the corresponding columns of M are copied into the value array    
//! of Mblock. M is not modified.
void SUBR(sdMat_get_block)(TYPE(const_sdMat_ptr) M, 
                             TYPE(sdMat_ptr) Mblock,
                             int imin, int imax, int jmin, int jmax, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  //TODO - do we need this in ghost?
  *ierr=-99;
  }

//! given a serial dense matrix Mblock, set M(imin:imax,jmin:jmax)=Mblock by 
//! copying the corresponding elements. Mblock is not modified.
void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) M, 
                             TYPE(const_sdMat_ptr) Mblock,
                             int imin, int imax, int jmin, int jmax, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=-99;
  }

//! \name destructors

//@{

//!
void SUBR(crsMat_delete)(TYPE(crsMat_ptr) vA, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_mat_t,A,vA,*ierr);
  A->destroy(A);
  }

//!
void SUBR(mvec_delete)(TYPE(mvec_ptr) vV, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,V,vV,*ierr);
  V->destroy(V);
  }

//!
void SUBR(sdMat_delete)(TYPE(sdMat_ptr) vM, int* ierr)
  {
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,M,vM,*ierr);
  M->destroy(M);
  }

//@}

//! \name Numerical functions
//!@{

//! put scalar value into all elements of a multi-vector
void SUBR(mvec_put_value)(TYPE(mvec_ptr) vV, _ST_ value, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,V,vV,*ierr);
  PHIST_OUT(9,"put value, V @ %p. V->traits->nrows=%ld\n",V,V->traits->nrows);
  V->fromScalar(V,(void*)&value);
  }

//! put scalar value into all elements of a multi-vector
void SUBR(sdMat_put_value)(TYPE(sdMat_ptr) vV, _ST_ value, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,V,vV,*ierr);
  PHIST_OUT(9,"sdMat put value, M @ %p. V->traits->nrows=%ld\n",V,V->traits->nrows);
  V->fromScalar(V,(void*)&value);
  }

//! put random numbers into all elements of a multi-vector
void SUBR(mvec_random)(TYPE(mvec_ptr) vV, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,V,vV,*ierr);
  V->fromRand(V);
  }

void SUBR(mvec_print)(TYPE(const_mvec_ptr) vV, int* ierr)
  {
  _CAST_PTR_FROM_VOID_(ghost_vec_t,V,vV,*ierr);
  V->print(V);
  }

void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) vM, int* ierr)
  {
  _CAST_PTR_FROM_VOID_(ghost_vec_t,M,vM,*ierr);
  M->print(M);
  }


//! put random numbers into all elements of a serial dense matrix
void SUBR(sdMat_random)(TYPE(sdMat_ptr) vM, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,M,vM,*ierr);
  M->fromRand(M);
  }

//! \name Numerical functions

//! compute the 2-norm) of each column of v                   
//! (vnrm[i] must be pre-allocated by caller)
  void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) vV,
                            _MT_* vnrm, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"  
  int i,nv;
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,V,vV,*ierr);  
  nv=V->traits->nvecs;
  _ST_ tmp[V->traits->nvecs];
  V->dotProduct(V,V,tmp);
  for (i=0;i<nv;i++) vnrm[i]=mt::sqrt(st::real(tmp[i]));
  return;
  }

  //! normalize (in the 2-norm) each column of v and return ||v||_2
  //! for each vector i in vnrm[i] (must be pre-allocated by caller)
  void SUBR(mvec_normalize)(TYPE(mvec_ptr) vV,
                            _MT_* vnrm, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"  
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,V,vV,*ierr);  
  // TODO - this call doesn't return the norm as we wish
  //V->normalize(V);    
  PHIST_CHK_IERR(SUBR(mvec_norm2)(vV,vnrm,ierr),*ierr);
  _ST_ inrm[V->traits->nvecs];
  for (int i=0;i<V->traits->nvecs;i++) inrm[i]=st::one()/vnrm[i];
  PHIST_CHK_IERR(SUBR(mvec_vscale)(vV,inrm,ierr),*ierr);
  return;
  }

//! scale each column i of v and by scalar[i]
void SUBR(mvec_scale)(TYPE(mvec_ptr) vV, 
                            _ST_ scalar, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,V,vV,*ierr);  
  V->scale(V,(void*)&scalar);
  return;
  }

//! scale each column i of v and by scalar[i]
void SUBR(mvec_vscale)(TYPE(mvec_ptr) vV, 
                            _ST_* scalar, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"  
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,V,vV,*ierr);  
  V->vscale(V,(void*)scalar);
  return;
  }

//! y=alpha*x+beta*y
void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) vX,
                            _ST_ beta,  TYPE(mvec_ptr)       vY, 
                            int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,X,vX,*ierr);
  _CAST_PTR_FROM_VOID_(ghost_vec_t,Y,vY,*ierr);
  Y->axpby(Y,X,(void*)&alpha,(void*)&beta);
  }

//! y[j]=alpha[j]*x[j]+beta[j]*y[j] for all columns j
void SUBR(mvec_vadd_mvec)(const _ST_ *alpha, TYPE(const_mvec_ptr) vX,
                            _ST_ beta,  TYPE(mvec_ptr)       vY, 
                            int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,X,vX,*ierr);
  _CAST_PTR_FROM_VOID_(ghost_vec_t,Y,vY,*ierr);
  Y->vaxpby(Y,X,(void*)&alpha,(void*)&beta);
  }

//! B=alpha*A+beta*B
void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA,
                            _ST_ beta,  TYPE(sdMat_ptr)       vB,
                            int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  SUBR(mvec_add_mvec)(alpha,vA,beta,vB,ierr);
  }

//! y=alpha*A*x+beta*y.
void SUBR(crsMat_times_mvec)(_ST_ alpha, TYPE(const_crsMat_ptr) vA, TYPE(const_mvec_ptr) vx, 
_ST_ beta, TYPE(mvec_ptr) vy, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_mat_t,A,vA,*ierr);
  _CAST_PTR_FROM_VOID_(ghost_vec_t,x,vx,*ierr);
  _CAST_PTR_FROM_VOID_(ghost_vec_t,y,vy,*ierr);
  if (alpha==st::zero())
    {
    // no MVM needed
    if (beta==st::zero())
      {
      PHIST_CHK_IERR(SUBR(mvec_put_value)(vy,beta,ierr),*ierr);
      }
    else if (beta!=st::one())
      {
      y->scale(y,(void*)&beta);
      }
    }
  else
    {
    int spMVM_opts=GHOST_SPMVM_DEFAULT;
    void* old_scale = A->traits->scale;
    if (alpha!=st::one())
      {
      // TODO: this fails for some reason!
      A->traits->scale = (void*)&alpha;
      spMVM_opts|=GHOST_SPMVM_APPLY_SCALE;
      // copy input vector and scale 
      /*
      x=x->clone(x,x->traits->nvecs,0);
      x->scale(x,(void*)&alpha);
      */
      }
    if (beta!=st::zero())
      {
      spMVM_opts|=GHOST_SPMVM_AXPY;
      if (beta!=st::one())
        {
        y->scale(y,&beta);
        }
      }
    *ierr=ghost_spmvm(A->context,y,A,x,&spMVM_opts);
    if (alpha!=st::one())
      {
      A->traits->scale = old_scale;
      //x->destroy(x); // x has been cloned
      }
    }
  }

//! dot product of vectors v_i and w_i, i=1..numvecs
void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) vV, TYPE(const_mvec_ptr) vW, _ST_* s, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(ghost_vec_t,W,vW,*ierr);
  V->dotProduct(V,W,s);
  }

//! dense tall skinny matrix-matrix product yielding a serial dense matrix
//! C=alpha*V'*W+beta*C. C is replicated on all MPI processes sharing V and W.
void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) vV, TYPE(const_mvec_ptr) vW, _ST_ beta, TYPE(sdMat_ptr) vC, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(ghost_vec_t,W,vW,*ierr);
  _CAST_PTR_FROM_VOID_(ghost_vec_t,C,vC,*ierr);
  char trans;
#ifdef _IS_COMPLEX_
  trans='C';
#else
  trans='T';
#endif  
  *ierr=ghost_gemm(&trans,V,W,C,(void*)&alpha,(void*)&beta,GHOST_GEMM_ALL_REDUCE);
  }


//! n x m multi-vector times m x m dense matrix gives n x m multi-vector,
//! W=alpha*V*C + beta*W
void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) vV,
                                       TYPE(const_sdMat_ptr) vC,
                           _ST_ beta,  TYPE(mvec_ptr) vW,
                                       int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(ghost_vec_t,C,vC,*ierr);
  _CAST_PTR_FROM_VOID_(ghost_vec_t,W,vW,*ierr);
  char trans='N';
  // note: C is replicated, so this operation is a purely local one.
  *ierr=ghost_gemm(&trans,V,C,W,(void*)&alpha,(void*)&beta,GHOST_GEMM_NO_REDUCE);
  }

//! n x m serial dense matrix times m x k serial dense matrix gives n x k sdMat,
//! C=alpha*V*W + beta*C (serial XGEMM wrapper)
void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vV,
                                         TYPE(const_sdMat_ptr) vW,
                              _ST_ beta, TYPE(sdMat_ptr) vC,
                                         int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(ghost_vec_t,W,vW,*ierr);
  _CAST_PTR_FROM_VOID_(ghost_vec_t,C,vC,*ierr);
  char trans='N';
  *ierr=ghost_gemm(&trans, V, W, C, (void*)&alpha, (void*)&beta, GHOST_GEMM_NO_REDUCE);
  }

//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.   
//! Q is computed in place of V. If V does not have full rank, ierr>0   
//! indicates the dimension of the null-space of V. The first m-ierr    
//! columns of Q are an orthogonal basis of the column space of V, the  
//! remaining columns form a basis for the null space.  
void SUBR(mvec_QR)(TYPE(mvec_ptr) vV, TYPE(sdMat_ptr) vR, int* ierr)
  {
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  _CAST_PTR_FROM_VOID_(ghost_vec_t,V,vV,*ierr);
  _CAST_PTR_FROM_VOID_(ghost_vec_t,R,vR,*ierr);

  int stride = R->traits->nrowspadded;
  int nrows = R->traits->nrows;
  int ncols = R->traits->nvecs;
    
  // wrapper class for ghost_vec_t for calling Belos.
  // The wrapper does not own the vector so it doesn't destroy it.
  phist::GhostMV mv_V(V,false);
    
#ifdef TESTING
  _CHECK_ZERO_(nrows-ncols,*ierr);
  _CHECK_ZERO_(nrows-(V->traits->nvecs),*ierr);
#endif  

  PHIST_OUT(9,"create Teuchos view of R");
  Teuchos::RCP<Traits<_ST_ >::Teuchos_sdMat_t> R_view;
  PHIST_CHK_IERR(R_view = Traits<_ST_ >::CreateTeuchosViewNonConst
        (Teuchos::rcp(R,false),ierr),*ierr);

  PHIST_OUT(9,"create TSQR ortho manager");  
  Belos::TsqrOrthoManager<_ST_, phist::GhostMV> tsqr("phist/ghost");
  Teuchos::RCP<const Teuchos::ParameterList> valid_params = 
        tsqr.getValidParameters();
  // faster but numerically less robust settings:
  Teuchos::RCP<const Teuchos::ParameterList> fast_params = 
        tsqr.getFastParameters();
  Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp
        (new Teuchos::ParameterList(*valid_params));
  params->set("randomizeNullSpace",true);
  PHIST_OUT(9,"set TSQR parameters");
  tsqr.setParameterList(params);

  int rank;
  _TRY_CATCH_(rank = tsqr.normalize(mv_V,R_view),*ierr);
  PHIST_OUT(9,"V has %d columns and rank %d",ncols,rank);
  *ierr = ncols-rank;// return positive number if rank not full.
  return;
  }


//!@}

}// extern "C"
