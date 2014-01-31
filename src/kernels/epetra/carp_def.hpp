extern "C" {

void SUBR(carp_create)(TYPE(const_crsMat_ptr) vA, TYPE(carpData)** dat_ptr, int* ierr)
{
#include "phist_std_typedefs.hpp"
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_CrsMatrix, A, vA, *ierr);
  const Epetra_Map& rowMap = A->RowMap();
  const Epetra_Map& colMap = A->ColMap();
  TYPE(carpData)* dat = new TYPE(carpData);
  
  Epetra_Vector* diagA = new Epetra_Vector(colMap);
  Epetra_Vector* rowScaling = new Epetra_Vector(colMap);
  
  PHIST_CHK_IERR(*ierr=A->ExtractDiagonalCopy(*diagA),*ierr);
  
  // compute squared two-norm of each row of A
  int* col;
  int len;
  double* val;
#pragma omp parallel for private(len,col,val) schedule(static)
  for (int i=0;i<A->NumMyRows();i++)
  {
    *ierr=A->ExtractMyRowView(i,len,val,col);
    (*rowScaling)[i]=0.0;
    for (int j=0;j<len;j++)
    {
      rowScaling[0][i] += val[j]*val[j];
    }
  }
  
  *dat_ptr=dat;
  dat->omega_=st::one();
  dat->rowScaling_=(TYPE(mvec_ptr))rowScaling;
  dat->diagA_=(TYPE(mvec_ptr))diagA;
  dat->xLoc_=NULL;
}

void SUBR(carp_delete)(TYPE(carpData)* dat, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  if (dat==NULL) return;
  *ierr=-99; //TODO - this can probably  be implemented in a kernel-lib independent way
}

void SUBR(carp_fb)(TYPE(carpData)* dat, TYPE(const_crsMat_ptr) vA, 
        TYPE(const_mvec_ptr) vB, _ST_ const * sigma, 
        TYPE(const_mvec_ptr) vrhs, TYPE(mvec_ptr) vsol, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  CAST_PTR_FROM_VOID(const Epetra_CrsMatrix, A, vA, *ierr);
  const Epetra_Vector* B = (const Epetra_Vector*)(vB);
  CAST_PTR_FROM_VOID(const Epetra_MultiVector, rhs, vrhs, *ierr);
  CAST_PTR_FROM_VOID(Epetra_MultiVector, sol, vsol, *ierr);
  
#ifdef TESTING
if ( (sol->Map().SameAs(rhs->Map())==false) ||
     (A->RowMap().SameAs(sol->Map())==false) )
     {
       PHIST_SOUT(PHIST_ERROR,"CARP: vectors and matrix do not match or case not "
       "implemented\n(file %s, line %d)\n",__FILE__,__LINE__);
     }
#endif  

  bool needVecs=false;
  bool needImport=true;
  
  if (sol->Map().Comm().NumProc()>1)
  {
    if (dat->xLoc_==NULL) 
    {
      needVecs=true;
    }
    else
    {
      CAST_PTR_FROM_VOID(Epetra_MultiVector,xLoc,dat->xLoc_,*ierr);
      if (xLoc->NumVectors()!=sol->NumVectors())
      {
        delete xLoc;
        needVecs=true;
      }
    }
  }
  else
  {
    // create a view that will be re-created in the next call and deleted in carp_delete
    PHIST_CHK_IERR(SUBR(mvec_view_block)(vsol,&dat->xLoc_,0,sol->NumVectors()-1,ierr),*ierr);
    needImport=false;
  }

  if (needVecs)
  {
    PHIST_SOUT(PHIST_DEBUG,"(re-)allocate temporary vector");
    dat->xLoc_ = (TYPE(mvec_ptr))(new Epetra_MultiVector(A->ColMap(),sol->NumVectors()));
  }
  
  // import the sol vector into the column map of A
  CAST_PTR_FROM_VOID(Epetra_MultiVector,xLoc,dat->xLoc_,*ierr);
  if (needImport)
  {
    PHIST_CHK_IERR(*ierr=xLoc->Import(*sol, *A->Importer(),Insert),*ierr);
  }
  CAST_PTR_FROM_VOID(const Epetra_Vector,diagA,dat->diagA_,*ierr);
  CAST_PTR_FROM_VOID(const Epetra_Vector,rowScaling,dat->rowScaling_,*ierr);

  int* col;
  int len;
  double* val;

  ///////////////////////////////////////////////////////////
  // forward Kaczmarz sweep                                //
  ///////////////////////////////////////////////////////////
  for (int i=0;i<A->NumMyRows();i++)
  {
    double Aii = (*diagA)[i];
    double Bii = 1.0;
    if (B!=NULL) Bii=(*B)[i];
    double nrm_ai = (*rowScaling)[i];
    PHIST_CHK_IERR(*ierr=A->ExtractMyRowView(i,len,val,col),*ierr);
    //TODO - switch loop order
    for (int k=0;k<sol->NumVectors();k++)
    {
      double factor=(*rhs)[k][i];
      // scale by 1/||A(i,:)||_2^2, correct row norm for shifted diagonal entry:
      // ||A(i,:) + sigmaB(i,i)||_2^2 = ||A(i,:)||_2^2 - A(i,i)^2 + (A(i,i)-sigma*B(i,i))^2
      //                              = ||A(i,:)||_2^2 -2*sigma*A(i,i)*B(i,i)
      double scal=dat->omega_/(nrm_ai-2.0*sigma[k]*Aii*Bii+sigma[k]*sigma[k]*Bii*Bii);
      for (int j=0;j<len;j++)
      {
        // note: xLoc lives in the column map of A, so we do not need to convert the col index
        factor += val[j]*(*xLoc)[k][col[j]];
      }//j

      for (int j=0;j<len;j++)
      {
        // Projection step: 
        // update all elements j in one step (this prevents straight-forward OpenMP usage,
        // we need a distance-2 coloring for intra-node parallelization here).
        (*xLoc)[k][col[j]] += factor*scal*val[j];
      }//j
    }//k
  }// i

  ///////////////////////////////////////////////////////////
  // communication/averaging                               //
  ///////////////////////////////////////////////////////////

  // average overlapping nodes into X. As explained by Gordon & Gordon (SISC 27, 2005), this
  // parallel algorithm (which they call CARP) is equivalent to Kaczmarz in a superspace of 
  // R^n in which the overlapping elements appear multiple times.
  if (needImport)
  {
    PHIST_CHK_IERR(*ierr=sol->Export(*xLoc, *(A->Importer()), Average),*ierr);
    PHIST_CHK_IERR(*ierr=xLoc->Import(*sol, *A->Importer(),Insert),*ierr);
  }

  ///////////////////////////////////////////////////////////
  // forward Kaczmarz sweep                                //
  ///////////////////////////////////////////////////////////
  for (int i=A->NumMyRows()-1;i>=0; i--)
  {
    double Aii = (*diagA)[i];
    double Bii = 1.0;
    if (B!=NULL) Bii=(*B)[i];
    double nrm_ai = (*rowScaling)[i];
    PHIST_CHK_IERR(*ierr=A->ExtractMyRowView(i,len,val,col),*ierr);
    //TODO - switch loop order
    for (int k=0;k<sol->NumVectors();k++)
    {
      double factor=(*rhs)[k][i];
      // scale by 1/||A(i,:)||_2^2, correct row norm for shifted diagonal entry:
      // ||A(i,:) + sigmaB(i,i)||_2^2 = ||A(i,:)||_2^2 - A(i,i)^2 + (A(i,i)-sigma*B(i,i))^2
      //                              = ||A(i,:)||_2^2 -2*sigma*A(i,i)*B(i,i)
      double scal=dat->omega_/(nrm_ai-2.0*sigma[k]*Aii*Bii+sigma[k]*sigma[k]*Bii*Bii);
      for (int j=0;j<len;j++)
      {
        // note: xLoc lives in the column map of A, so we do not need to convert the col index
        factor += val[j]*(*xLoc)[k][col[j]];
      }//j

      for (int j=0;j<len;j++)
      {
        // Projection step: 
        // update all elements j in one step (this prevents straight-forward OpenMP usage,
        // we need a distance-2 coloring for intra-node parallelization here).
        (*xLoc)[k][col[j]] += factor*scal*val[j];
      }//j
    }//k
  }// i

  if (needImport)
  {
    PHIST_CHK_IERR(*ierr=sol->Export(*xLoc, *(A->Importer()), Average),*ierr);
  }

}//carp_fb

} // extern "C"
