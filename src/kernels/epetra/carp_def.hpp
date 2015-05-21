//! kernels to implement CARP for the matrix sigma[j]I-A.

typedef struct {
Teuchos::RCP<Epetra_Map> colMap_;
Teuchos::RCP<Epetra_MultiVector> xLoc_;
} carpWork_t;

void private_kacz_sweep(const Epetra_CrsMatrix* A, 
                              Epetra_MultiVector* xLoc, 
                        const Epetra_MultiVector* rhs,
                        double omega,
                        int imin, int iend, int istep,
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
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        _MT_ **nrms_ai2i, void** work, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  bool anyComplex=false;
  for (int i=0;i<numShifts;i++)
  {
    if (sigma_i[i]!=0.0) anyComplex=true;
  }
  if (anyComplex) 
  {
    PHIST_SOUT(PHIST_ERROR,"only real variant of CARP kernel implemented in Epetra right now\n");
    *iflag=-99;
    return;
  }

  PHIST_CAST_PTR_FROM_VOID(const Epetra_CrsMatrix, A, vA, *iflag);
  
  const Epetra_Map& rowMap = A->RowMap();
  const Epetra_Map& colMap = A->ColMap();
  
  if (!rowMap.SameAs(A->RangeMap()) || 
      !rowMap.SameAs(A->DomainMap()))
      {
        PHIST_SOUT(PHIST_ERROR,"Epetra CARP implementation currently needs a simple square matrix.\n");
        *iflag=-99;
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

  *work=(void*)dat;
}

//! forward/backward sweep of Kaczmarz/CARP algorithm (SSOR sweep on the normal equations),
//! with matrix sigma[j]*I-A applied to the columns of mvec X[j] with a single rhs B. The
//! input arguments nrms_ai2i and work must be unchanged from the _setup routine. For each
//! shift sigma[j], a separate relaxation parameter omega[j] must be provided.
extern "C" void SUBR(carp_sweep)(TYPE(const_sparseMat_ptr) vA, int numShifts, 
        _MT_ const sigma_r[], _MT_ const sigma_i[],
        TYPE(const_mvec_ptr) vrhs,
        TYPE(mvec_ptr) X_r[], TYPE(mvec_ptr) X_i[],
        _MT_ const* nrm_ai2i, void* const work,
        _MT_ const * omega, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Epetra_CrsMatrix, A, vA, *iflag);
  const Epetra_MultiVector* rhs=(const Epetra_MultiVector*)vrhs;
  
    carpWork_t* dat=(carpWork_t*)work;


  for (int isys=0;isys<numShifts;isys++)
  {
    if (sigma_r[isys]!=0.0||sigma_i[isys]!=0.0)
    {
      PHIST_SOUT(PHIST_WARNING,"CARP: ignoring shift (not implemented in Epetra)\n");
    }
  
    PHIST_CAST_PTR_FROM_VOID(Epetra_MultiVector, sol, X_r[isys], *iflag);

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

    if (sol->Map().Comm().NumProc()==1)
    {
      dat->xLoc_=Teuchos::rcp(sol,false);
    }
    else if (dat->xLoc_!=Teuchos::null)
    {
      if (sol->NumVectors()!=dat->xLoc_->NumVectors())
      {
        dat->xLoc_=Teuchos::null;
      }
    }
      
  
    if (dat->xLoc_==Teuchos::null)
    {
      PHIST_SOUT(PHIST_DEBUG,"(re-)allocate temporary vector");
      dat->xLoc_ = Teuchos::rcp(new Epetra_MultiVector(A->ColMap(),sol->NumVectors()));
    }
  
    // import the sol vector into the column map of A
    if (dat->xLoc_.get()!=sol)
    {
      PHIST_CHK_IERR(*iflag=dat->xLoc_->Import(*sol, *(A->Importer()),Insert),*iflag);
    }

    ///////////////////////////////////////////////////////////
    // local forward Kaczmarz sweep                         //
    ///////////////////////////////////////////////////////////

    PHIST_CHK_IERR(private_kacz_sweep(A, dat->xLoc_.get(),rhs,omega[isys],
                0,A->NumMyRows(),+1,iflag),*iflag);

    ///////////////////////////////////////////////////////////
    // component averaging                                   //
    ///////////////////////////////////////////////////////////

    // average overlapping nodes into X. As explained by Gordon & Gordon (SISC 27, 2005), this
    // parallel algorithm (which they call CARP) is equivalent to Kaczmarz in a superspace of 
    // R^n in which the overlapping elements appear multiple times.
    if (dat->xLoc_.get()!=sol)
    {
      PHIST_CHK_IERR(*iflag=sol->Export(*dat->xLoc_, *(A->Importer()), Average),*iflag);
      PHIST_CHK_IERR(*iflag=dat->xLoc_->Import(*sol, *A->Importer(),Insert),*iflag);
    }

    ///////////////////////////////////////////////////////////
    // local backward Kaczmarz sweep                         //
    ///////////////////////////////////////////////////////////

    PHIST_CHK_IERR(private_kacz_sweep(A, dat->xLoc_.get(),rhs,omega[isys],
                A->NumMyRows()-1,-1,-1,iflag),*iflag);

    if (dat->xLoc_.get()!=sol)
    {
      PHIST_CHK_IERR(*iflag=sol->Export(*dat->xLoc_, *(A->Importer()), Average),*iflag);
    }
  }
}

//! clean up data structures created by carp_setup
extern "C" void SUBR(carp_destroy)(TYPE(const_sparseMat_ptr) A, int numShifts,
        _MT_* nrms_ai2i, void* work, int *iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);
  *iflag=0;
  PHIST_CAST_PTR_FROM_VOID(carpWork_t,dat,work,*iflag);
  dat->xLoc_=Teuchos::null;
  delete dat;
}

void private_kacz_sweep(const Epetra_CrsMatrix* A, 
                              Epetra_MultiVector* xLoc, 
                        const Epetra_MultiVector* rhs,
                        double omega,
                        int istart, int iend, int istep,
                        int *iflag)
{
  *iflag=0;
  int* col;
  int len;
  double* val;
  
  int nvecs=xLoc->NumVectors();

  for (int i=istart; i!=iend;i+=istep)
  {
    PHIST_CHK_IERR(*iflag=A->ExtractMyRowView(i,len,val,col),*iflag);

    // row 2-norm
    double nrm_ai2 = 0.0;
    for (int j=0;j<len;j++)
    {
      nrm_ai2+=val[j]*val[j];
    }
    
    // first compute scaling factors. Note that
    // the interface tells us to use sI-A as operator,
    // see also issue #119.
    double tmp_r[nvecs];
    for (int k=0;k<nvecs;k++)
    {
      tmp_r[k]=0.0;
      for (int j=0;j<len;j++)
      {
        // note: xLoc lives in the column map of A, so we do not need to convert the col index
        tmp_r[k] -= val[j]*(*xLoc)[k][col[j]];
      }//j
      if (rhs!=NULL) tmp_r[k]-=(*rhs)[k][i];
      tmp_r[k] *= omega/nrm_ai2;
    }//k

    // Projection step:
    for (int k=0;k<nvecs;k++)
    {
      for (int j=0;j<len;j++)
      {
        (*xLoc)[k][col[j]] += tmp_r[k]*val[j];
      }//j
    }//k
  }// i
}
