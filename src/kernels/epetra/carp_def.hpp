extern "C" {

SUBR(carp_create)(TYPE(const_crsMat_ptr) vA, TYPE(carpData)** dat_ptr, int* ierr)
{
  ENTER_FCN(__FUNCTION__);
  *ierr=0;
  CAST_PTR_FROM_VOID(const Epetra_CrsMatrix, A, vA, *ierr);
  const Epetra_Map& rowMap = A->RowMap();
  const Epetra_Map& colMap = A->ColMap();
  TYPE(carpData)* dat = new TYPE(carpData);
  
  Epetra_Vector* diagA = new Epetra_Vector(colMap);
  Epetra_Vector* rowScaling = new Epetra_Vector(colMap);
  
  PHIST_CHK_IERR(*ierr=A.ExtractDiagonalCopy(*diagA),*ierr);
  
  // compute squared two-norm of each row of A
  int* col;
  int len;
  double* val
#pragma omp parallel private(len,col,val) schedule(static)  
  for (int i=0;i<A.NumMyRows();i++)
  {
    PHIST_CHK_IERR(*ierr=A.ExtractMyRowView(i,len,val,col),*ierr);
    (*rowScaling)[i]=0.0;
    for (int j=0;j<len;j++)
    {
      rowScaling[i] += val[j]*val[j];
    }
  }
  
  *dat_ptr=dat;
  dat->rowScaling_=(TYPE(mvec_ptr))rowScaling;
  dat->diagA_=(TYPE(mvec_ptr))diagA;
  dat->xLoc_=NULL;
}

SUBR(carp_delete)(TYPE(carpData)* dat, int* ierr)
{
  *ierr=0;
  if (dat==NULL) return;
  *ierr=-99;
}

SUBR(carp_fb)(TYPE(carpData)* dat, TYPE(const_crsMat_ptr) vA, 
        TYPE(const_mvec_ptr) vB, _ST_ const * sigma, 
        TYPE(const_mvec_ptr) vrhs, TYPE(mvec_ptr) vsol, int* ierr)
{
  CAST_PTR_FROM_VOID(const Epetra_CrsMatrix, A, vA, *ierr);
  const Epetra_Vector* B = (const Epetra_Vector*)(vB);
  CAST_PTR_FROM_VOID(const Epetra_MultiVector, rhs, vrhs, *ierr);
  CAST_PTR_FROM_VOID(Epetra_MultiVector, sol, vsol, *ierr);
  
#ifdef TESTING
if ( (X->Map().SameAs(B->Map())==false) ||
     (A->RowMap().SameAs(X->Map())==false) )
     {
       PHIST_SOUT(PHIST_ERROR("CARP: vectors and matrix do not match or case not "
       "implemented\n(file %s, line %d)\n",__FILE__,__LINE__);
     }
#endif  

  bool needVecs=false;
  if (sol->Map().Comm().NumProc()>1)
  {
    if (carp->xLoc==NULL || carp->bLoc_==NULL) 
    {
      needVecs=true;
    }
    else
    {
      if (dat->Xloc_->NumVectors()!=sol->NumVectors())
      {
        needVecs=true;
      }
    }
  }

  if (needVecs)
  {
    dat->Xloc_ = new Epetra_MultiVector(A->ColMap(),sol->NumVectors());
  }
  
  // import the sol vector into the column map of A
  PHIST_CHK_IERR(*ierr=dat->xLoc->Import(*sol, *A->Importer(),Insert),ierr);

  int* col;
  int len;
  double* val

  ///////////////////////////////////////////////////////////
  // forward Kaczmarz sweep                                //
  ///////////////////////////////////////////////////////////
  for (int i=0;i<A.NumMyRows();i++)
  {
    double Aii = *(dat->diagA_)[i];
    double Bii = 1.0;
    if (B!=NULL) Bii=(*B)[i];
    double nrm_ai = (*(dat->rowScaling_))[i];
    PHIST_CHK_IERR(*ierr=A.ExtractMyRowView(i,len,val,col),*ierr);
    for (int k=0;k<X->NumVectors();k++)
    {
      double factor=(*rhs)[k][i];
      // scale by 1/||A(i,:)||_2^2, correct row norm for shifted diagonal entry:
      // ||A(i,:) + sigmaB(i,i)||_2^2 = ||A(i,:)||_2^2 - A(i,i)^2 + (A(i,i)-sigma*B(i,i))^2
      //                              = ||A(i,:)||_2^2 -2*sigma*A(i,i)*B(i,i)
      double cal=omega/(nrm_ai-2.0*sigma[k]*Aii*Bii+sigma[k]*sigma[k]*Bii*Bii);
      for (int j=0;j<len;j++)
      {
        // note: xLoc lives in the column map of A, so we do not need to convert the col index
        factor += val[j]*(*(dat->xLoc_))[k][cols[j]];

        // Projection step: 
        // update all elements j in one step (this prevents straight-forward OpenMP usage,
        // we need a distance-2 coloring for intra-node parallelization here).
        (*(dat->xLoc_))[k][cols[j]] += factor*scal*val[j];
      }//j
    }//k
  }// i

  ///////////////////////////////////////////////////////////
  // communication/averaging                               //
  ///////////////////////////////////////////////////////////

  // average overlapping nodes into X. As explained by Gordon & Gordon (SISC 27, 2005), this
  // parallel algorithm (which they call CARP) is equivalent to Kaczmarz in a superspace of 
  // R^n in which the overlapping elements appear multiple times.
  PHIST_CHK_IERR(*ierr=X->Export(*(dat->xLoc_), *(A->Importer()), Average),*ierr);
  PHIST_CHK_IERR(*ierr=dat->xLoc->Import(*sol, *A->Importer(),Insert),*ierr);

  ///////////////////////////////////////////////////////////
  // forward Kaczmarz sweep                                //
  ///////////////////////////////////////////////////////////
  for (int i=A.NumMyRows()-1;i>=0; i--)
  {
    double Aii = *(dat->diagA_)[i];
    double Bii = 1.0;
    if (B!=NULL) Bii=(*B)[i];
    double nrm_ai = (*(dat->rowScaling_))[i];
    PHIST_CHK_IERR(*ierr=A.ExtractMyRowView(i,len,val,col),*ierr);
    for (int k=0;k<X->NumVectors();k++)
    {
      double factor=(*rhs)[k][i];
      // scale by 1/||A(i,:)||_2^2, correct row norm for shifted diagonal entry:
      // ||A(i,:) + sigmaB(i,i)||_2^2 = ||A(i,:)||_2^2 - A(i,i)^2 + (A(i,i)-sigma*B(i,i))^2
      //                              = ||A(i,:)||_2^2 -2*sigma*A(i,i)*B(i,i)
      double cal=omega/(nrm_ai-2.0*sigma[k]*Aii*Bii+sigma[k]*sigma[k]*Bii*Bii);
      for (int j=0;j<len;j++)
      {
        // note: xLoc lives in the column map of A, so we do not need to convert the col index
        factor += val[j]*(*(dat->xLoc_))[k][cols[j]];

        // Projection step: 
        // update all elements j in one step (this prevents straight-forward OpenMP usage,
        // we need a distance-2 coloring for intra-node parallelization here).
        (*(dat->xLoc_))[k][cols[j]] += factor*scal*val[j];
      }//j
    }//k
  }// i


}//carp_fb

} // extern "C"
