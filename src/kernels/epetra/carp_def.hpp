//! kernels to implement CARP for the matrix sigma[j]I-A.

typedef struct {
Teuchos::RCP<Epetra_Vector> invProcWeights_;
Teuchos::RCP<Epetra_MultiVector> xLoc_;
Teuchos::RCP<Epetra_MultiVector> xLoc_i_;
} carpWork_t;

// local KACZ sweep on (A-sigma*I)x=b, b may be NULL
void private_kacz_sweep(const Epetra_CrsMatrix* A, 
                        const double sigma_r[],
                              Epetra_MultiVector* xLoc, 
                        const Epetra_MultiVector* rhs,
                        double const omega[],
                        int istart, int iend, int istep,
                        int *iflag);

// local KACZ sweep on [A-sigma*I, Q; Q' 0][xLoc; q]=[b; 0], b may be NULL
void private_kacz_sweep_aug(const Epetra_CrsMatrix* A, 
                        const double sigma_r[],
                        const Epetra_MultiVector* Q,
                              Epetra_MultiVector* xLoc, 
                              Epetra_MultiVector* q,
                        const Epetra_MultiVector* rhs,
                        double const omega[],
                        int istart, int iend, int istep,
                        int *iflag);

// local KACZ sweep on (A-(sigma_r+i*sigma_i)I) (x_r+i*x_i)=b, where b is real and may be NULL
void private_kacz_sweep_rc(const Epetra_CrsMatrix* A, 
                        const double sigma_r[],
                        const double sigma_i[],
                              Epetra_MultiVector* xLoc_r, 
                              Epetra_MultiVector* xLoc_i, 
                        const Epetra_MultiVector* rhs,
                        double const omega[],
                        int istart, int iend, int istep,
                        int *iflag);

// local KACZ sweep on [(A-(sigma_r+i*sigma_i)I), Q; Q' 0] [(x_r+i*x_i); q_r+i*q_i]=[b; 0], where b is real and may be NULL
void private_kacz_sweep_aug_rc(const Epetra_CrsMatrix* A, 
                        const double sigma_r[],
                        const double sigma_i[],
                        const Epetra_MultiVector* Q,
                              Epetra_MultiVector* xLoc_r, 
                              Epetra_MultiVector* xLoc_i, 
                              Epetra_MultiVector* q_r, 
                              Epetra_MultiVector* q_i, 
                        const Epetra_MultiVector* rhs,
                        double const omega[],
                        int istart, int iend, int istep,
                        int *iflag);

//! create data structures needed for subsequent calls to carp_sweep: 
//! 
//! input: 
//!
//! numShifts: number of shifts sigma 
//! sigma_r, sigma_i: possibly complex shifts sigma=sigma_r+i*sigma_i
//!
//! output:
//! 
//! *nrms_ai2i and *work should be NULL on input and not touched while
//! carp_sweep is being called. If no longer needed, they should be cleaned
//! up using carp_destroy.
//!
//! If the shifts sigma change, carp_destroy and carp_setup should be
//! used to rebuild the working objects.

extern "C" void SUBR(carp_setup)(TYPE(const_sparseMat_ptr) vA, int numShifts, 
        _MT_ const sigma[],
        void** work, int* iflag)
{
  PHIST_CHK_IERR(SUBR(carp_setup_rc)(vA,numShifts,sigma,NULL,work,iflag),*iflag);
}

extern "C" void SUBR(carp_setup_rc)(TYPE(const_sparseMat_ptr) vA, int numShifts, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        void** work, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;

  PHIST_CAST_PTR_FROM_VOID(const Epetra_CrsMatrix, A, vA, *iflag);
  
  const Epetra_Map& rowMap = A->RowMap();
  const Epetra_Map& colMap = A->ColMap();
  
  if (!rowMap.SameAs(A->RangeMap()) || 
      !rowMap.SameAs(A->DomainMap()))
      {
        PHIST_SOUT(PHIST_ERROR,"Epetra CARP implementation currently needs a simple square matrix.\n");
        *iflag=PHIST_NOT_IMPLEMENTED;
        return;
      }
  
/* commented out because we prefer to compute the norms on-the-fly
  Epetra_Vector* diagA = new Epetra_Vector(colMap);
  PHIST_CHK_IERR(*iflag=A->ExtractDiagonalCopy(*diagA),*iflag);

  // compute squared two-norm of each row of A
  int* col;
  int len;
  double* val;

  for (int i=0;i<A->NumMyRows();i++)
  {
    *iflag=A->ExtractMyRowView(i,len,val,col);
    (*rowScaling)[i]=0.0;
    for (int j=0;j<len;j++)
    {
      rowScaling[0][i] += val[j]*val[j];
    }
    // we do not include the diagonal in this scaling, it must
    // explicitly be added in carp_fb
    rowScaling[0][i] -= (*diagA)[i]*(*diagA)[i];
  }
  */

  carpWork_t *dat=new carpWork_t;
  dat->xLoc_=Teuchos::null; // (re)allocated in sweep kernel depending on #rhs
  dat->xLoc_i_=Teuchos::null; 
  
  dat->invProcWeights_=Teuchos::rcp(new Epetra_Vector(A->ColMap()));
  
  PHIST_CHK_IERR(SUBR(mvec_put_value)((void*)dat->invProcWeights_.get(),st::one(),iflag),*iflag);
  if (A->Importer()!=NULL)
  {
    Teuchos::RCP<Epetra_Vector> locWeights = Teuchos::rcp(new Epetra_Vector(A->RowMap()));
    PHIST_CHK_IERR(*iflag=locWeights->Export(*dat->invProcWeights_,*(A->Importer()),Add),*iflag);
    PHIST_CHK_IERR(*iflag=dat->invProcWeights_->Import(*locWeights,*(A->Importer()),Insert),*iflag);
    for (int i=0; i<dat->invProcWeights_->MyLength(); i++)
    {
      double& ipw_i = (*dat->invProcWeights_)[i];
      ipw_i=1.0/ipw_i;
    }
  }
  *work=(void*)dat;
}

//! forward/backward sweep of Kaczmarz/CARP algorithm (SSOR sweep on the normal equations),
//! with matrix sigma[j]*I-A applied to the columns of mvec X[j] with a single rhs B. The
//! input arguments nrms_ai2i and work must be unchanged from the _setup routine. For each
//! shift sigma[j], a separate relaxation parameter omega[j] must be provided.
extern "C" void SUBR(carp_sweep)(TYPE(const_sparseMat_ptr) vA, 
        _ST_ const sigma[],
        TYPE(const_mvec_ptr) vrhs,
        TYPE(mvec_ptr) X,
        void* const work,
        _MT_ const * omega, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  // implementation is in one place in the aug_rc variant, which will call
  // the appropriate KACZ kernel and take care of the averaging (CARP part)
  PHIST_CHK_IERR(SUBR(carp_sweep_aug_rc)(vA, sigma, NULL, NULL, vrhs, X, NULL,
                NULL, NULL, work, omega, iflag), *iflag);
}

extern "C" void SUBR(carp_sweep_rc)(TYPE(const_sparseMat_ptr) vA, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) vrhs,
        TYPE(mvec_ptr) X_r, TYPE(mvec_ptr) X_i,
        void* const work,
        _MT_ const * omega, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(carp_sweep_aug_rc)(vA, sigma_r, sigma_i, NULL, vrhs, X_r, X_i,
                NULL, NULL, work, omega, iflag), *iflag);
}

void SUBR(carp_sweep_aug)(TYPE(const_sparseMat_ptr) vA,
        _ST_ const sigma[],
        TYPE(const_mvec_ptr) vQ,
        TYPE(const_mvec_ptr) vrhs,
        TYPE(mvec_ptr) X,
        TYPE(sdMat_ptr) vq,
        void* const work,
        _MT_ const * omega, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CHK_IERR(SUBR(carp_sweep_aug_rc)(vA, sigma,NULL, NULL, vrhs, X, NULL,
                vq, NULL, work, omega, iflag), *iflag);
}

void SUBR(carp_sweep_aug_rc)(TYPE(const_sparseMat_ptr) vA,
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) vQ,
        TYPE(const_mvec_ptr) vrhs,
        TYPE(mvec_ptr) X_r, TYPE(mvec_ptr) X_i,
        TYPE(sdMat_ptr) vq_r, TYPE(sdMat_ptr) vq_i,
        void* const work,
        _MT_ const * omega, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Epetra_CrsMatrix, A, vA, *iflag);
  const Epetra_MultiVector* rhs=(const Epetra_MultiVector*)vrhs;
  const Epetra_MultiVector* Q=(const Epetra_MultiVector*)vQ;
        Epetra_MultiVector* q_r=(Epetra_MultiVector*)vq_r;
        Epetra_MultiVector* q_i=(Epetra_MultiVector*)vq_i;
  
    carpWork_t* dat=(carpWork_t*)work;

  int nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(X_r,&nvec,iflag),*iflag);

  bool rc_variant=(X_i!=NULL); // variant with real A and rhs but complex diagonal shift sigma
  bool aug_variant=(Q!=NULL);
  
  if (aug_variant)
  {
    if (q_r==NULL || (q_i==NULL&&rc_variant))
    {
      PHIST_SOUT(PHIST_ERROR,"augmented CARP variant requires q_r (and q_i if X_i is given)");
      *iflag=PHIST_INVALID_INPUT;
      return;
    }
#ifdef TESTING
    double norms[Q->NumVectors()];
    PHIST_CHK_IERR(*iflag=Q->Norm2(norms),*iflag);
    double norm_err=0.0;
    for (int i=0; i<Q->NumVectors(); i++) norm_err=std::max(std::abs(1.0-norms[i]),norm_err);
    if (norm_err>1.0e-14) PHIST_SOUT(PHIST_WARNING,"augmented CARP kernel assumes Q is normalized, found max_j(abs(||Q(:,j)||_2-1.0))=%e\n",norm_err);
#endif
  }
  
  PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector, sol, X_r, *iflag);
  Epetra_MultiVector* sol_i=NULL;
  if (rc_variant)
  {
    sol_i=(Epetra_MultiVector*)X_i;
    if (sigma_i==NULL)
    {
      *iflag=PHIST_BAD_CAST;
      PHIST_SOUT(PHIST_ERROR,"sigma_i was NULL on input, but X_i was given\n");
    }
  }

    if (rhs!=NULL)
    {
      if (sol->Map().SameAs(rhs->Map())==false)
      {
        PHIST_SOUT(PHIST_ERROR,"CARP: sol and rhs do not match \n"
           "(file %s, line %d)\n",__FILE__,__LINE__);
      }
    }
    if (A->RowMap().SameAs(sol->Map())==false)
    {
      PHIST_SOUT(PHIST_ERROR,"CARP: vectors and matrix do not match or case not "
       "implemented\n(file %s, line %d)\n",__FILE__,__LINE__);
    }

    if (A->Importer()==NULL)
    {
      dat->xLoc_=Teuchos::rcp(sol,false);
      if (rc_variant) dat->xLoc_i_=Teuchos::rcp(sol_i,false);
      else dat->xLoc_i_ = Teuchos::null;
    }
    else if (dat->xLoc_!=Teuchos::null)
    {
      if (sol->NumVectors()!=dat->xLoc_->NumVectors())
      {
        dat->xLoc_=Teuchos::null;
        dat->xLoc_i_=Teuchos::null;
      }
    }
      
  
    if (dat->xLoc_==Teuchos::null)
    {
      PHIST_SOUT(PHIST_DEBUG,"(re-)allocate temporary import vector(s)");
      dat->xLoc_ = Teuchos::rcp(new Epetra_MultiVector(A->ColMap(),sol->NumVectors()));
    }
    if (sol_i!=NULL && dat->xLoc_i_==Teuchos::null)
    {
      dat->xLoc_i_ = Teuchos::rcp(new Epetra_MultiVector(A->ColMap(),sol_i->NumVectors()));
    }
  
    // import the sol vector into the column map of A
    if (dat->xLoc_.get()!=sol)
    {
      PHIST_CHK_IERR(*iflag=dat->xLoc_->Import(*sol, *(A->Importer()),Insert),*iflag);
    }
    if (rc_variant && dat->xLoc_i_.get()!=sol_i)
    {
      PHIST_CHK_IERR(*iflag=dat->xLoc_i_->Import(*sol_i, *(A->Importer()),Insert),*iflag);
    }

    ///////////////////////////////////////////////////////////
    // local forward Kaczmarz sweep                         //
    ///////////////////////////////////////////////////////////

    if (rc_variant)
    {
      if (aug_variant)
      {
        PHIST_CHK_IERR(private_kacz_sweep_aug_rc(A, sigma_r,sigma_i,Q,
                  dat->xLoc_.get(), dat->xLoc_i_.get(), q_r, q_i,
                  rhs,omega,0,A->NumMyRows(),+1,iflag),*iflag);
      }
      else
      {
        PHIST_CHK_IERR(private_kacz_sweep_rc(A, sigma_r,sigma_i,
                  dat->xLoc_.get(), dat->xLoc_i_.get(),
                  rhs,omega,0,A->NumMyRows(),+1,iflag),*iflag);
      }
    }
    else
    {
      if (aug_variant)
      {
        PHIST_CHK_IERR(private_kacz_sweep_aug(A, sigma_r, Q,
                  dat->xLoc_.get(), q_r, rhs,omega,
                  0,A->NumMyRows(),+1,iflag),*iflag);
      }
      else
      {
        PHIST_CHK_IERR(private_kacz_sweep(A, sigma_r,
                  dat->xLoc_.get(),rhs,omega,
                  0,A->NumMyRows(),+1,iflag),*iflag);
      }
    }
    ///////////////////////////////////////////////////////////
    // component averaging                                   //
    ///////////////////////////////////////////////////////////

    // average overlapping nodes into X. As explained by Gordon & Gordon (SISC 27, 2005), this
    // parallel algorithm (which they call CARP) is equivalent to Kaczmarz in a superspace of 
    // R^n in which the overlapping elements appear multiple times.
    if (dat->xLoc_.get()!=sol)
    {
      for (int i=0; i<dat->xLoc_->MyLength();i++)
      {
        for (int j=0; j< dat->xLoc_->NumVectors(); j++)
        {
          (*dat->xLoc_)[j][i]*=(*dat->invProcWeights_)[i];
          if (rc_variant)
          {
            (*dat->xLoc_i_)[j][i]*=(*dat->invProcWeights_)[i];
          }
        }
      }
      PHIST_CHK_IERR(*iflag=sol->Export(*dat->xLoc_, *(A->Importer()), Add),*iflag);
      PHIST_CHK_IERR(*iflag=dat->xLoc_->Import(*sol, *A->Importer(),Insert),*iflag);
    }
    if (rc_variant && dat->xLoc_i_.get()!=sol_i)
    {
      PHIST_CHK_IERR(*iflag=sol_i->Export(*dat->xLoc_i_, *(A->Importer()), Average),*iflag);
      PHIST_CHK_IERR(*iflag=dat->xLoc_i_->Import(*sol_i, *A->Importer(),Insert),*iflag);
    }

    ///////////////////////////////////////////////////////////
    // local backward Kaczmarz sweep                         //
    ///////////////////////////////////////////////////////////
    if (rc_variant)
    {
      if (aug_variant)
      {
        PHIST_CHK_IERR(private_kacz_sweep_aug_rc(A,sigma_r,sigma_i,Q,
                  dat->xLoc_.get(),dat->xLoc_i_.get(),
                  q_r, q_i,
                  rhs,omega,A->NumMyRows()-1,-1,-1,iflag),*iflag);
      }
      else
      {
        PHIST_CHK_IERR(private_kacz_sweep_rc(A,sigma_r,sigma_i,
                  dat->xLoc_.get(),dat->xLoc_i_.get(),
                  rhs,omega,A->NumMyRows()-1,-1,-1,iflag),*iflag);
      }
    }
    else
    {
      if (aug_variant)
      {
        PHIST_CHK_IERR(private_kacz_sweep_aug(A,sigma_r,Q,
                  dat->xLoc_.get(), q_r,
                  rhs,omega,A->NumMyRows()-1,-1,-1,iflag),*iflag);
      }
      else
      {
        PHIST_CHK_IERR(private_kacz_sweep(A,sigma_r,
                  dat->xLoc_.get(),
                  rhs,omega,A->NumMyRows()-1,-1,-1,iflag),*iflag);
      }
    }

    if (dat->xLoc_.get()!=sol)
    {
      for (int i=0; i<dat->xLoc_->MyLength();i++)
      {
        for (int j=0; j< dat->xLoc_->NumVectors(); j++)
        {
          (*dat->xLoc_)[j][i]*=(*dat->invProcWeights_)[i];
          if (rc_variant)
          {
            (*dat->xLoc_i_)[j][i]*=(*dat->invProcWeights_)[i];
          }
        }
      }
      PHIST_CHK_IERR(*iflag=sol->Export(*dat->xLoc_, *(A->Importer()), Add),*iflag);
    }
    if (rc_variant && dat->xLoc_i_.get()!=sol_i)
    {
      PHIST_CHK_IERR(*iflag=sol_i->Export(*dat->xLoc_i_, *(A->Importer()), Average),*iflag);
    }
}


//! clean up data structures created by carp_setup
extern "C" void SUBR(carp_destroy)(TYPE(const_sparseMat_ptr) A,
        void* work, int *iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(carpWork_t,dat,work,*iflag);
  dat->xLoc_=Teuchos::null;
  dat->xLoc_i_=Teuchos::null;
  delete dat;
}

void private_kacz_sweep(const Epetra_CrsMatrix* A,
                        double const sigma_r[],
                              Epetra_MultiVector* xLoc, 
                        const Epetra_MultiVector* rhs,
                        double const omega[],
                        int istart, int iend, int istep,
                        int *iflag)
{
  *iflag=0;
  int* col;
  int len;
  double* val;
  
  int nvec=xLoc->NumVectors();

  for (int i=istart; i!=iend;i+=istep)
  {
    PHIST_CHK_IERR(*iflag=A->ExtractMyRowView(i,len,val,col),*iflag);

    // row 2-norm
    double nrm_ai2[nvec];
    double tmp_r[nvec];
    double d=0.0,tmp=0.0;
    int col_i=A->LCID(A->GRID64(i));
    
    for (int j=0;j<len;j++)
    {
      tmp+=val[j]*val[j];
      if (col[j]==col_i)
      {
        d=val[j]; // diagonal element
      }
    }
    for (int k=0;k<nvec;k++)
    {
      nrm_ai2[k]=tmp-2*d*sigma_r[k]+sigma_r[k]*sigma_r[k];
    }
    
    // first compute scaling factors. 
    for (int k=0;k<nvec;k++)
    {
      tmp_r[k] = -sigma_r[k]*(*xLoc)[k][col_i];
      for (int j=0;j<len;j++)
      {
        // note: xLoc lives in the column map of A, so we do not need to convert the col index
        tmp_r[k] += val[j]*(*xLoc)[k][col[j]];
      }//j
      if (rhs!=NULL) tmp_r[k]-=(*rhs)[k][i];
      tmp_r[k] *= omega[k]/nrm_ai2[k];
    }//k

    // Projection step:
    for (int k=0;k<nvec;k++)
    {
      for (int j=0;j<len;j++)
      {
        (*xLoc)[k][col[j]] -= tmp_r[k]*val[j];
      }//j
      (*xLoc)[k][col_i]+=tmp_r[k]*sigma_r[k];
    }//k
  }// i
}

// local KACZ sweep on [A-sigma*I, Q; Q' 0][xLoc; q]=[b; 0], b may be NULL
void private_kacz_sweep_aug(const Epetra_CrsMatrix* A, 
                        const double sigma_r[],
                        const Epetra_MultiVector* Q,
                              Epetra_MultiVector* xLoc, 
                              Epetra_MultiVector* q,
                        const Epetra_MultiVector* rhs,
                        double const omega[],
                        int istart, int iend, int istep,
                        int *iflag)
{
  *iflag=0;
  int* col;
  int len;
  double* val;
  
  int nvec=xLoc->NumVectors();
  int nvecQ=Q->NumVectors();

  if (istep>0)
  {
    // forward sweep, orthogonalize against Q first.    
    // Note that xLoc has a different map than Q, so we 
    // have to import it first (at most a local permutation)
    Teuchos::RCP<Epetra_MultiVector> x=Teuchos::rcp(xLoc, false);
    if (A->Importer()!=NULL)
    {
      x=Teuchos::rcp(new Epetra_MultiVector(Q->Map(),xLoc->NumVectors()));
      PHIST_CHK_IERR(*iflag=x->Import(*xLoc,*(A->Importer()),Insert),*iflag);
    }
    Epetra_MultiVector qx=*q;
    // q_tmp = Q'x
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(1.0, Q, x.get(), 0.0, &qx, iflag), *iflag);
    // scale col k by omega[k]
    for (int k=0; k<qx.NumVectors(); k++)
    {
      qx(k)->Scale(omega[k]);
    }
    // x = x - Q*qx
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-1.0,Q,&qx,1.0,x.get(),iflag),*iflag);
  }

  for (int i=istart; i!=iend;i+=istep)
  {
    PHIST_CHK_IERR(*iflag=A->ExtractMyRowView(i,len,val,col),*iflag);

    // row 2-norm
    double nrm_ai2[nvec];
    double tmp_r[nvec];
    double d=0.0,tmp=0.0;
    int col_i=A->LCID(A->GRID64(i));
    
    for (int j=0;j<len;j++)
    {
      tmp+=val[j]*val[j];
      if (col[j]==col_i)
      {
        d=val[j]; // diagonal element
      }
    }
    // adjust for augmentation mvec Q
    for (int j=0; j<nvecQ; j++)
    {
      double qij=(*Q)[j][i];
      tmp+=qij*qij;
    }
    for (int k=0;k<nvec;k++)
    {
      // adjust for diagonal shift
      nrm_ai2[k]=tmp-2*d*sigma_r[k]+sigma_r[k]*sigma_r[k];
    }
    
    // first compute scaling factors. 
    for (int k=0;k<nvec;k++)
    {
      tmp_r[k] = -sigma_r[k]*(*xLoc)[k][col_i];
      for (int j=0;j<len;j++)
      {
        // note: xLoc lives in the column map of A, so we do not need to convert the col index
        tmp_r[k] += val[j]*(*xLoc)[k][col[j]];
      }//j
      for (int j=0; j<nvecQ; j++)
      {
        tmp_r[k] += (*Q)[j][i]*(*q)[k][j];
      }
      if (rhs!=NULL) tmp_r[k]-=(*rhs)[k][i];
      tmp_r[k] *= omega[k]/nrm_ai2[k];
    }//k

    // Projection step:
    for (int k=0;k<nvec;k++)
    {
      for (int j=0;j<len;j++)
      {
        (*xLoc)[k][col[j]] -= tmp_r[k]*val[j];
      }//j
      (*xLoc)[k][col_i]+=tmp_r[k]*sigma_r[k];
      for (int j=0; j<nvecQ; j++)
      {
        (*q)[k][j]+=tmp_r[k]*(*Q)[j][i];
      }
    }//k
  }// i

  if (istep<0)
  {
    // backward sweep, orthogonalize against Q last.
    Teuchos::RCP<Epetra_MultiVector> x=Teuchos::rcp(xLoc, false);
    if (A->Importer()!=NULL)
    {
      x=Teuchos::rcp(new Epetra_MultiVector(Q->Map(),xLoc->NumVectors()));
      PHIST_CHK_IERR(*iflag=x->Import(*xLoc,*(A->Importer()),Insert),*iflag);
    }
    Epetra_MultiVector qx=*q;
    // q_tmp = Q'x
    PHIST_CHK_IERR(SUBR(mvecT_times_mvec)(1.0, Q, x.get(), 0.0, &qx, iflag), *iflag);
    // scale col k by omega[k]
    for (int k=0; k<qx.NumVectors(); k++)
    {
      qx(k)->Scale(omega[k]);
    }
    // x = x - Q*qx
    PHIST_CHK_IERR(SUBR(mvec_times_sdMat)(-1.0,Q,&qx,1.0,x.get(),iflag),*iflag);
  }
}

void private_kacz_sweep_rc(const Epetra_CrsMatrix* A,
                        double const sigma_r[],
                        double const sigma_i[],
                              Epetra_MultiVector* xLoc, 
                              Epetra_MultiVector* xLoc_i, 
                        const Epetra_MultiVector* rhs,
                        double const omega[],
                        int istart, int iend, int istep,
                        int *iflag)
{
  *iflag=0;
  int* col;
  int len;
  double* val;
  
  int nvec=xLoc->NumVectors();

  for (int i=istart; i!=iend;i+=istep)
  {
    PHIST_CHK_IERR(*iflag=A->ExtractMyRowView(i,len,val,col),*iflag);

    // row 2-norm
    double tmp_r[nvec],tmp_i[nvec];
    double nrm_ai2[nvec];
    double d=0.0,tmp=0.0;
    int col_i = A->LCID(A->GRID64(i));

    for (int j=0;j<len;j++)
    {
      tmp+=val[j]*val[j];
      if (col[j]==col_i)
      {
        d=val[j];
      }
    }
    for (int k=0;k<nvec;k++)
    {
      nrm_ai2[k]=tmp-2*d*sigma_r[k]+sigma_r[k]*sigma_r[k]
                    +sigma_i[k]*sigma_i[k];
    }
    
    // first compute scaling factors. 
    for (int k=0;k<nvec;k++)
    {
      tmp_r[k]= -sigma_r[k]*(*xLoc)[k][col_i] 
                +sigma_i[k]*(*xLoc_i)[k][col_i];
      tmp_i[k]=-sigma_i[k]*(*xLoc)[k][i]
               -sigma_r[k]*(*xLoc_i)[k][i];
      for (int j=0;j<len;j++)
      {
        // note: xLoc lives in the column map of A, so we do not need to convert the col index
        tmp_r[k] += val[j]*(*xLoc)[k][col[j]];
        tmp_i[k] += val[j]*(*xLoc_i)[k][col[j]];
      }//j
      if (rhs!=NULL) tmp_r[k]-=(*rhs)[k][i];
      tmp_r[k] *= omega[k]/nrm_ai2[k];
      tmp_i[k] *= omega[k]/nrm_ai2[k];
    }//k

    // Projection step:
    for (int k=0;k<nvec;k++)
    {
      for (int j=0;j<len;j++)
      {
        (*xLoc)[k][col[j]] -= tmp_r[k]*val[j];
        (*xLoc_i)[k][col[j]] -= tmp_i[k]*val[j];
      }//j
      (*xLoc)[k][col_i]+=tmp_r[k]*sigma_r[k]+tmp_i[k]*sigma_i[k];
      (*xLoc_i)[k][col_i]+=tmp_i[k]*sigma_r[k]-tmp_r[k]*sigma_i[k];
    }//k
  }// i
}

// local KACZ sweep on [(A-(sigma_r+i*sigma_i)I), Q; Q' 0] [(x_r+i*x_i); q_r+i*q_i]=[b; 0], where b is real and may be NULL
void private_kacz_sweep_aug_rc(const Epetra_CrsMatrix* A, 
                        const double sigma_r[],
                        const double sigma_i[],
                        const Epetra_MultiVector* Q,
                              Epetra_MultiVector* xLoc_r, 
                              Epetra_MultiVector* xLoc_i, 
                              Epetra_MultiVector* q_r, 
                              Epetra_MultiVector* q_i, 
                        const Epetra_MultiVector* rhs,
                        double const omega[],
                        int istart, int iend, int istep,
                        int *iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=-99;
  return;
}
