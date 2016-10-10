/*! \file petsc/kernels_def.hpp
 * included by petsc/kernels.cpp
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 * \author "Jonas Thies <Jonas.Thies@DLR.de>
 *
*/

extern "C" void SUBR(type_avail)(int *iflag)
{
  *iflag=0;
}

extern "C" void SUBR(sparseMat_read_mm)(TYPE(sparseMat_ptr)* A, phist_const_comm_ptr vcomm,
        const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  std::string line;
  std::ifstream infile(filename);
  getline(infile, line);
  if( line != "%%MatrixMarket matrix coordinate real general" &&  line != "%%MatrixMarket matrix coordinate real symmetric" )
  {
    *iflag = PHIST_NOT_IMPLEMENTED;
    return;
  }

  while(infile.peek() == '%')
    getline(infile, line);

  phist_gidx globalRows, globalCols, globalLines;
  infile >> globalRows >> globalCols >> globalLines;
  phist_map_ptr map;
  PHIST_CHK_IERR(phist_map_create(&map,vcomm,globalRows,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(sparseMat_read_mm_with_map)(A, map, filename, iflag), *iflag);

  *iflag=PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_read_bin)(TYPE(sparseMat_ptr)* A, phist_const_comm_ptr vcomm,
const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_hb)(TYPE(sparseMat_ptr)* A, phist_const_comm_ptr vcomm,
const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_mm_with_map)(TYPE(sparseMat_ptr)* vA, phist_const_map_ptr map,
        const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sparseMat_t*,A,vA,*iflag);

  std::string line;
  std::ifstream infile(filename);
  getline(infile, line);
  bool symm;
  if( line == "%%MatrixMarket matrix coordinate real general" )
    symm = false;
  else if( line == "%%MatrixMarket matrix coordinate real symmetric" )
    symm = true;
  else
  {
    *iflag = PHIST_NOT_IMPLEMENTED;
    return;
  }

  while(infile.peek() == '%')
    getline(infile, line);

  phist_gidx globalRows, globalCols, globalLines;
  infile >> globalRows >> globalCols >> globalLines;
  phist_gidx globalEntries = globalLines;
  if( symm )
    globalEntries += globalLines - globalRows;
  phist_lidx avg_nne = globalEntries/globalRows;

  phist_lidx nlocal;
  phist_gidx ilower, iupper, nglob;
  phist_const_comm_ptr comm;
  PHIST_CHK_IERR(phist_map_get_local_length(map, &nlocal, iflag), *iflag);
  PHIST_CHK_IERR(phist_map_get_ilower(map, &ilower, iflag), *iflag);
  PHIST_CHK_IERR(phist_map_get_iupper(map, &iupper, iflag), *iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(map, &nglob, iflag), *iflag);
  PHIST_CHK_IERR(phist_map_get_comm(map, &comm, iflag), *iflag);

  PHIST_CHK_IERR( *iflag = PetscNew(A), *iflag);
  (*A)->map = map;
  PHIST_SOUT(PHIST_INFO, "reading mm mat: %d x %d\n", nglob, nglob);
  //PHIST_CHK_IERR( *iflag = MatCreateAIJ(*(MPI_Comm*)comm, nlocal, nlocal, nglob, nglob, avg_nne, NULL, avg_nne, NULL, &((*A)->m)), *iflag);
  PHIST_CHK_IERR( *iflag = MatCreate(*(MPI_Comm*)comm, &((*A)->m)), *iflag);
  PHIST_CHK_IERR( *iflag = MatSetType((*A)->m, MATAIJ), *iflag);
  PHIST_CHK_IERR( *iflag = MatSetSizes((*A)->m, nlocal, nlocal, nglob, nglob), *iflag);
  MatType matType;
  PHIST_CHK_IERR( *iflag = MatGetType((*A)->m, &matType), *iflag);
  if( std::string(matType) == std::string(MATMPIAIJ) )
  {
    PHIST_CHK_IERR( *iflag = MatMPIAIJSetPreallocation((*A)->m, avg_nne, NULL, avg_nne, NULL), *iflag);
  }
  else if( std::string(matType) == std::string(MATSEQAIJ) )
  {
    PHIST_CHK_IERR( *iflag = MatSeqAIJSetPreallocation((*A)->m, avg_nne, NULL), *iflag);
  }
  else
  {
    PHIST_SOUT(PHIST_ERROR, "strange PETSc matrix type %s\n", matType);
    *iflag = PHIST_NOT_IMPLEMENTED;
    return;
  }
  PHIST_CHK_IERR( *iflag = MatSetOption((*A)->m, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE), *iflag);


  for (phist_gidx i = 0; i < globalLines; i++)
  {
    _ST_ val;
    phist_gidx row, col;
    infile >> row >> col >> val;
    // count from zero...
    row--, col--;
    if( ilower <= row && row <= iupper )
    {
      PHIST_CHK_IERR( *iflag = MatSetValue((*A)->m,row,col,val,INSERT_VALUES), *iflag);
    }
    if( symm && col != row && ilower <= col && col <= iupper )
    {
      PHIST_CHK_IERR( *iflag = MatSetValue((*A)->m,col,row,val,INSERT_VALUES), *iflag);
    }
  }
  PHIST_CHK_IERR( *iflag = MatAssemblyBegin((*A)->m,MAT_FINAL_ASSEMBLY), *iflag);
  PHIST_CHK_IERR( *iflag = MatAssemblyEnd((*A)->m,MAT_FINAL_ASSEMBLY), *iflag);

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_read_bin_with_map)(TYPE(sparseMat_ptr)* A, phist_const_map_ptr map,
        const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_read_hb_with_map)(TYPE(sparseMat_ptr)* A, phist_const_map_ptr map,
        const char* filename,int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}
extern "C" void SUBR(sparseMat_get_row_map)(TYPE(const_sparseMat_ptr) vA, phist_const_map_ptr* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t,A,vA,*iflag);

  *map = A->map;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_get_col_map)(TYPE(const_sparseMat_ptr) vA, phist_const_map_ptr* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t,A,vA,*iflag);

  *map = A->map;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_get_domain_map)(TYPE(const_sparseMat_ptr) vA, phist_const_map_ptr* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t,A,vA,*iflag);

  *map = A->map;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_get_range_map)(TYPE(const_sparseMat_ptr) vA, phist_const_map_ptr* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t,A,vA,*iflag);

  *map = A->map;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_create)(TYPE(mvec_ptr)* vV, 
    phist_const_map_ptr map, phist_lidx nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t*,V,vV,*iflag);

  phist_lidx nlocal;
  phist_gidx nglob;
  phist_const_comm_ptr comm;
  PHIST_CHK_IERR(phist_map_get_local_length(map,&nlocal,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(map,&nglob,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_comm(map,&comm,iflag),*iflag);
  PHIST_CHK_IERR( *iflag = PetscNew(V), *iflag);
  (*V)->map = map;
  PHIST_CHK_IERR( *iflag = PetscCalloc1(nlocal*nvec, &((*V)->rawData)), *iflag);
  PHIST_CHK_IERR( *iflag = MatCreate(*(MPI_Comm*)comm, &((*V)->v)), *iflag);
  PHIST_CHK_IERR( *iflag = MatSetType((*V)->v, MATDENSE), *iflag);
  PHIST_CHK_IERR( *iflag = MatSetSizes((*V)->v, nlocal, PETSC_DECIDE, nglob, nvec), *iflag);
  MatType matType;
  PHIST_CHK_IERR( *iflag = MatGetType((*V)->v, &matType), *iflag);
  if( std::string(matType) == std::string(MATMPIDENSE) )
  {
    PHIST_CHK_IERR( *iflag = MatMPIDenseSetPreallocation((*V)->v, (PetscScalar*)(*V)->rawData), *iflag);
  }
  else if( std::string(matType) == std::string(MATSEQDENSE) )
  {
    PHIST_CHK_IERR( *iflag = MatSeqDenseSetPreallocation((*V)->v, (PetscScalar*)(*V)->rawData), *iflag);
  }
  else
  {
    PHIST_SOUT(PHIST_ERROR, "strange PETSc matrix type %s\n", matType);
    *iflag = PHIST_NOT_IMPLEMENTED;
    return;
  }
  PHIST_CHK_IERR( *iflag = MatAssemblyBegin((*V)->v, MAT_FINAL_ASSEMBLY), *iflag);
  PHIST_CHK_IERR( *iflag = MatAssemblyEnd((*V)->v, MAT_FINAL_ASSEMBLY), *iflag);
  (*V)->is_view = false;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_create_view)(TYPE(mvec_ptr)* vV, phist_const_map_ptr map, 
    _ST_* values, phist_lidx lda, int nvec,
    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t*,V,vV,*iflag);

  phist_lidx nlocal;
  phist_gidx nglob;
  phist_const_comm_ptr comm;
  PHIST_CHK_IERR(phist_map_get_local_length(map,&nlocal,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(map,&nglob,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_comm(map,&comm,iflag),*iflag);
  PHIST_CHK_IERR( *iflag = (lda==nlocal) ? PHIST_SUCCESS : PHIST_NOT_IMPLEMENTED, *iflag);
  PHIST_CHK_IERR( *iflag = PetscNew(V), *iflag);
  (*V)->map = map;
  (*V)->rawData = values;
  PHIST_CHK_IERR( *iflag = MatCreate(*(MPI_Comm*)comm, &((*V)->v)), *iflag);
  PHIST_CHK_IERR( *iflag = MatSetType((*V)->v, MATDENSE), *iflag);
  PHIST_CHK_IERR( *iflag = MatSetSizes((*V)->v, nlocal, PETSC_DECIDE, nglob, nvec), *iflag);
  MatType matType;
  PHIST_CHK_IERR( *iflag = MatGetType((*V)->v, &matType), *iflag);
  if( std::string(matType) == std::string(MATMPIDENSE) )
  {
    PHIST_CHK_IERR( *iflag = MatMPIDenseSetPreallocation((*V)->v, (PetscScalar*)(*V)->rawData), *iflag);
  }
  else if( std::string(matType) == std::string(MATSEQDENSE) )
  {
    PHIST_CHK_IERR( *iflag = MatSeqDenseSetPreallocation((*V)->v, (PetscScalar*)(*V)->rawData), *iflag);
  }
  else
  {
    PHIST_SOUT(PHIST_ERROR, "strange PETSc matrix type %s\n", matType);
    *iflag = PHIST_NOT_IMPLEMENTED;
    return;
  }

  PHIST_CHK_IERR( *iflag = MatAssemblyBegin((*V)->v, MAT_FINAL_ASSEMBLY), *iflag);
  PHIST_CHK_IERR( *iflag = MatAssemblyEnd((*V)->v, MAT_FINAL_ASSEMBLY), *iflag);
  (*V)->is_view = true;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_create)(TYPE(sdMat_ptr)* vM, 
    int nrows, int ncols, phist_const_comm_ptr comm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t*,M,vM,*iflag);

  PHIST_CHK_IERR( *iflag = PetscNew(M), *iflag);
  (*M)->comm = comm;
  PHIST_CHK_IERR( *iflag = PetscCalloc1(nrows*ncols, &((*M)->rawData)), *iflag);
  (*M)->lda = nrows;
  PHIST_CHK_IERR( *iflag = MatCreate(PETSC_COMM_SELF, &((*M)->m)), *iflag);
  PHIST_CHK_IERR( *iflag = MatSetType((*M)->m, MATSEQDENSE), *iflag);
  PHIST_CHK_IERR( *iflag = MatSetSizes((*M)->m, nrows, ncols, nrows, ncols), *iflag);
  PHIST_CHK_IERR( *iflag = MatSeqDenseSetLDA((*M)->m, (*M)->lda), *iflag);
  PHIST_CHK_IERR( *iflag = MatSeqDenseSetPreallocation((*M)->m, (PetscScalar*)(*M)->rawData), *iflag);
  PHIST_CHK_IERR( *iflag = MatAssemblyBegin((*M)->m, MAT_FINAL_ASSEMBLY), *iflag);
  PHIST_CHK_IERR( *iflag = MatAssemblyEnd((*M)->m, MAT_FINAL_ASSEMBLY), *iflag);
  (*M)->is_view = false;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_create_view)(TYPE(sdMat_ptr)* vM, phist_const_comm_ptr comm,
        _ST_* values, phist_lidx lda, int nrows, int ncols,
        int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t*,M,vM,*iflag);

  PHIST_CHK_IERR( *iflag = PetscNew(M), *iflag);
  (*M)->comm = comm;
  (*M)->rawData = values;
  (*M)->lda = lda;
  PHIST_CHK_IERR( *iflag = MatCreate(PETSC_COMM_SELF, &((*M)->m)), *iflag);
  PHIST_CHK_IERR( *iflag = MatSetType((*M)->m, MATSEQDENSE), *iflag);
  PHIST_CHK_IERR( *iflag = MatSetSizes((*M)->m, nrows, ncols, nrows, ncols), *iflag);
  PHIST_CHK_IERR( *iflag = MatSeqDenseSetLDA((*M)->m, (*M)->lda), *iflag);
  PHIST_CHK_IERR( *iflag = MatSeqDenseSetPreallocation((*M)->m, (PetscScalar*)(*M)->rawData), *iflag);
  PHIST_CHK_IERR( *iflag = MatAssemblyBegin((*M)->m, MAT_FINAL_ASSEMBLY), *iflag);
  PHIST_CHK_IERR( *iflag = MatAssemblyEnd((*M)->m, MAT_FINAL_ASSEMBLY), *iflag);
  (*M)->is_view = true;
  *iflag = PHIST_SUCCESS;
}
                  

extern "C" void SUBR(mvec_get_map)(TYPE(const_mvec_ptr) vV, phist_const_map_ptr* map, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  *map = V->map;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_num_vectors)(TYPE(const_mvec_ptr) vV, int* nvec, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);

  phist_lidx nrows;
  PHIST_CHK_IERR( *iflag = MatGetSize(V->v, &nrows, nvec), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_get_nrows)(TYPE(const_sdMat_ptr) vM, int* nrows, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,M,vM,*iflag);

  phist_lidx ncols;
  PHIST_CHK_IERR( *iflag = MatGetSize(M->m, nrows, &ncols), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_get_ncols)(TYPE(const_sdMat_ptr) vM, int* ncols, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,M,vM,*iflag);

  phist_lidx nrows;
  PHIST_CHK_IERR( *iflag = MatGetSize(M->m, &nrows, ncols), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_extract_view)(TYPE(mvec_ptr) vV, _ST_** val, phist_lidx* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);

  *val = V->rawData;
  phist_lidx nvec;
  PHIST_CHK_IERR( *iflag = MatGetLocalSize(V->v, lda, &nvec), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_extract_view)(TYPE(sdMat_ptr) vM, _ST_** val, phist_lidx* lda, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,M,vM,*iflag);
  *val = M->rawData;
  *lda = M->lda;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_to_mvec)(TYPE(const_mvec_ptr) vV, TYPE(mvec_ptr) vW, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,W,vW,*iflag);

  PHIST_CHK_IERR( *iflag = MatCopy(V->v, W->v, SAME_NONZERO_PATTERN), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_view_block)(TYPE(mvec_ptr) vV,
    TYPE(mvec_ptr)* vVblock,
    int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t*,Vblock,vVblock,*iflag);

  phist_lidx nvec = jmax-jmin+1;
  phist_lidx nlocal;
  phist_gidx nglob;
  phist_const_comm_ptr comm;
  PHIST_CHK_IERR(phist_map_get_local_length(V->map,&nlocal,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(V->map,&nglob,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_comm(V->map,&comm,iflag),*iflag);
  if( *vVblock )
  {
    PHIST_CHK_IERR(SUBR(mvec_delete)(*vVblock, iflag), *iflag);
  }
  PHIST_CHK_IERR( *iflag = PetscNew(Vblock), *iflag);
  (*Vblock)->map = V->map;
  (*Vblock)->rawData = V->rawData + jmin*nlocal;
  PHIST_CHK_IERR( *iflag = MatCreate(*(MPI_Comm*)comm, &((*Vblock)->v)), *iflag);
  PHIST_CHK_IERR( *iflag = MatSetType((*Vblock)->v, MATDENSE), *iflag);
  PHIST_CHK_IERR( *iflag = MatSetSizes((*Vblock)->v, nlocal, PETSC_DECIDE, nglob, nvec), *iflag);
  MatType matType;
  PHIST_CHK_IERR( *iflag = MatGetType((*Vblock)->v, &matType), *iflag);
  if( std::string(matType) == std::string(MATMPIDENSE) )
  {
    PHIST_CHK_IERR( *iflag = MatMPIDenseSetPreallocation((*Vblock)->v, (PetscScalar*)(*Vblock)->rawData), *iflag);
  }
  else if( std::string(matType) == std::string(MATSEQDENSE) )
  {
    PHIST_CHK_IERR( *iflag = MatSeqDenseSetPreallocation((*Vblock)->v, (PetscScalar*)(*Vblock)->rawData), *iflag);
  }
  else
  {
    PHIST_SOUT(PHIST_ERROR, "strange PETSc matrix type %s\n", matType);
    *iflag = PHIST_NOT_IMPLEMENTED;
    return;
  }

  PHIST_CHK_IERR( *iflag = MatAssemblyBegin((*Vblock)->v, MAT_FINAL_ASSEMBLY), *iflag);
  PHIST_CHK_IERR( *iflag = MatAssemblyEnd((*Vblock)->v, MAT_FINAL_ASSEMBLY), *iflag);
  (*Vblock)->is_view = true;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_get_block)(TYPE(const_mvec_ptr) vV,
    TYPE(mvec_ptr) vVblock,
    int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  TYPE(mvec_ptr) tmpV = NULL;
  PHIST_CHK_IERR(SUBR(mvec_view_block)((TYPE(mvec_ptr))vV, &tmpV, jmin, jmax, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_to_mvec)(tmpV, vVblock, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(tmpV, iflag), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_set_block)(TYPE(mvec_ptr) vV,
    TYPE(const_mvec_ptr) vVblock,
    int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  TYPE(mvec_ptr) tmpV = NULL;
  PHIST_CHK_IERR(SUBR(mvec_view_block)(vV, &tmpV, jmin, jmax, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_to_mvec)(vVblock, tmpV, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(tmpV, iflag), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_view_block)(TYPE(sdMat_ptr) vM, 
    TYPE(sdMat_ptr)* vMblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t*,Mblock,vMblock,*iflag);

  phist_lidx nrows = imax-imin+1;
  phist_lidx ncols = jmax-jmin+1;
  if( *vMblock )
  {
    PHIST_CHK_IERR(SUBR(sdMat_delete)(*vMblock, iflag), *iflag);
  }
  PHIST_CHK_IERR( *iflag = PetscNew(Mblock), *iflag);
  (*Mblock)->comm = M->comm;
  (*Mblock)->rawData = M->rawData + jmin*M->lda + imin;
  (*Mblock)->lda = M->lda;
  PHIST_CHK_IERR( *iflag = MatCreate(PETSC_COMM_SELF, &((*Mblock)->m)), *iflag);
  PHIST_CHK_IERR( *iflag = MatSetType((*Mblock)->m, MATSEQDENSE), *iflag);
  PHIST_CHK_IERR( *iflag = MatSetSizes((*Mblock)->m, nrows, ncols, nrows, ncols), *iflag);
  PHIST_CHK_IERR( *iflag = MatSeqDenseSetLDA((*Mblock)->m, (*Mblock)->lda), *iflag);
  PHIST_CHK_IERR( *iflag = MatSeqDenseSetPreallocation((*Mblock)->m, (PetscScalar*)(*Mblock)->rawData), *iflag);
  PHIST_CHK_IERR( *iflag = MatAssemblyBegin((*Mblock)->m, MAT_FINAL_ASSEMBLY), *iflag);
  PHIST_CHK_IERR( *iflag = MatAssemblyEnd((*Mblock)->m, MAT_FINAL_ASSEMBLY), *iflag);
  (*Mblock)->is_view = true;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_get_block)(TYPE(const_sdMat_ptr) vM, 
    TYPE(sdMat_ptr) vMblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,M,vM,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,Mblock,vMblock,*iflag);

  TYPE(sdMat_ptr) vtmpM = NULL;
  PHIST_CHK_IERR(SUBR(sdMat_view_block)((TYPE(sdMat_ptr))vM, &vtmpM, imin, imax, jmin, jmax, iflag), *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,tmpM,vtmpM,*iflag);
  PHIST_CHK_IERR( *iflag = MatCopy(tmpM->m, Mblock->m, SAME_NONZERO_PATTERN), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(vtmpM, iflag), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_set_block)(TYPE(sdMat_ptr) vM, 
    TYPE(const_sdMat_ptr) vMblock,
    int imin, int imax, int jmin, int jmax, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,Mblock,vMblock,*iflag);

  TYPE(sdMat_ptr) vtmpM = NULL;
  PHIST_CHK_IERR(SUBR(sdMat_view_block)(vM, &vtmpM, imin, imax, jmin, jmax, iflag), *iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,tmpM,vtmpM,*iflag);
  PHIST_CHK_IERR( *iflag = MatCopy(Mblock->m, tmpM->m, SAME_NONZERO_PATTERN), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_delete)(vtmpM, iflag), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_delete)(TYPE(sparseMat_ptr) vA, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  if( vA == NULL )
  {
    *iflag = PHIST_SUCCESS;
    return;
  }
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sparseMat_t,A,vA,*iflag);

  PHIST_CHK_IERR( *iflag = MatDestroy(&(A->m)), *iflag);
  PHIST_CHK_IERR( *iflag = PetscFree(A), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_delete)(TYPE(mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  if( vV == NULL )
  {
    *iflag = PHIST_SUCCESS;
    return;
  }
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);

  PHIST_CHK_IERR( *iflag = MatDestroy(&(V->v)), *iflag);
  if( !V->is_view )
  {
    PHIST_CHK_IERR( *iflag = PetscFree(V->rawData), *iflag);
  }
  PHIST_CHK_IERR( *iflag = PetscFree(V), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_delete)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  if( vM == NULL )
  {
    *iflag = PHIST_SUCCESS;
    return;
  }
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);

  PHIST_CHK_IERR( *iflag = MatDestroy(&(M->m)), *iflag);
  if( !M->is_view )
  {
    PHIST_CHK_IERR( *iflag = PetscFree(M->rawData), *iflag);
  }
  PHIST_CHK_IERR( *iflag = PetscFree(M), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_put_value)(TYPE(mvec_ptr) vV, _ST_ value, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);
  
  // no real support for this in Petsc...
  phist_lidx nvec, nlocal;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vV, &nvec, iflag), *iflag);
  PHIST_CHK_IERR(phist_map_get_local_length(V->map,&nlocal,iflag),*iflag);
  for(phist_lidx j = 0; j < nvec; j++)
    for(phist_lidx i = 0; i < nlocal; i++)
      V->rawData[j*nlocal+i] = value;
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_put_func)(TYPE(mvec_ptr) vV,
        phist_mvec_elemFunc funPtr,void* last_arg, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);

  // no real support for this in Petsc...
  phist_lidx nvec, nlocal;
  phist_gidx ilower;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vV, &nvec, iflag), *iflag);
  PHIST_CHK_IERR(phist_map_get_local_length(V->map,&nlocal,iflag),*iflag);
  PHIST_CHK_IERR(phist_map_get_ilower(V->map,&ilower,iflag),*iflag);
  for(phist_lidx j = 0; j < nvec; j++)
  {
    for(phist_lidx i = 0; i < nlocal; i++)
    {
      PHIST_CHK_IERR(*iflag=funPtr(ilower+i,j,(void*)&(V->rawData[j*nlocal+i]),last_arg),*iflag);
    }
  }
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_put_value)(TYPE(sdMat_ptr) vM, _ST_ value, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);

  // no real support for this in Petsc...
  phist_lidx nrows, ncols;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(vM, &nrows, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(vM, &ncols, iflag), *iflag);
  for(phist_lidx j = 0; j < ncols; j++)
    for(phist_lidx i = 0; i < nrows; i++)
      M->rawData[j*M->lda+i] = value;
  *iflag = PHIST_SUCCESS;
}

#ifndef PHIST_BUILTIN_RNG
extern "C" void SUBR(mvec_random)(TYPE(mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);

  PHIST_CHK_IERR( *iflag = MatSetRandom(V->m,NULL), *iflag);
  *iflag = PHIST_SUCCESS;
}
#endif

extern "C" void SUBR(mvec_print)(TYPE(const_mvec_ptr) vV, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  
  PHIST_CHK_IERR( *iflag = MatView(V->v,PETSC_VIEWER_STDOUT_WORLD), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_print)(TYPE(const_sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,M,vM,*iflag);

  PHIST_CHK_IERR( *iflag = MatView(M->m,PETSC_VIEWER_STDOUT_SELF), *iflag);
  *iflag = PHIST_SUCCESS;
}

#ifndef PHIST_BUILTIN_RNG
extern "C" void SUBR(sdMat_random)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);

  PHIST_CHK_IERR( *iflag = MatSetRandom(M->m,NULL), *iflag);
  *iflag = PHIST_SUCCESS;
}
#endif

extern "C" void SUBR(sdMat_identity)(TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);

  // no real support for this in Petsc...
  phist_lidx nrows, ncols;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(vM, &nrows, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(vM, &ncols, iflag), *iflag);
  for(phist_lidx j = 0; j < ncols; j++)
    for(phist_lidx i = 0; i < nrows; i++)
      M->rawData[j*M->lda+i] = (i==j) ? st::one() : st::zero();
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_norm2)(TYPE(const_mvec_ptr) vV,
    _MT_* vnrm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);

  // invalidate Petsc cached norms
  PHIST_CHK_IERR( *iflag = PetscObjectStateIncrease((PetscObject)V->v), *iflag);

  PHIST_CHK_IERR( *iflag = MatGetColumnNorms(V->v, NORM_2, (PetscReal*) vnrm), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_normalize)(TYPE(mvec_ptr) vV,
    _MT_* vnrm, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);

  phist_lidx nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vV, &nvec, iflag), *iflag);
  _ST_ * vnrm_inv = NULL;
  PHIST_CHK_IERR( *iflag = PetscMalloc1(nvec, &vnrm_inv), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_norm2)(vV, vnrm, iflag), *iflag);
  for(phist_lidx i = 0; i < nvec; i++)
    vnrm_inv[i] = st::one() / vnrm[i];
  PHIST_CHK_IERR(SUBR(mvec_vscale)(vV, vnrm_inv, iflag), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_scale)(TYPE(mvec_ptr) vV, 
    _ST_ scalar, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);

  PHIST_CHK_IERR( *iflag = MatScale(V->v, scalar), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_vscale)(TYPE(mvec_ptr) vV, 
    const _ST_* scalar, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  // no real support for this in petsc, scale columns one by one
  phist_lidx nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vV, &nvec, iflag), *iflag);
  TYPE(mvec_ptr) col = NULL;
  for(phist_lidx i = 0; i < nvec; i++)
  {
    PHIST_CHK_IERR(SUBR(mvec_view_block)(vV, &col, i, i, iflag), *iflag);
    PHIST_CHK_IERR(SUBR(mvec_scale)(col, scalar[i], iflag), *iflag);
  }
  PHIST_CHK_IERR(SUBR(mvec_delete)(col, iflag), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_add_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) vV,
                                    _ST_ beta,  TYPE(mvec_ptr)       vW,
                                    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,W,vW,*iflag);

  PHIST_CHK_IERR( *iflag = MatScale(W->v, beta), *iflag);
  PHIST_CHK_IERR( *iflag = MatAXPY(W->v, alpha, V->v, SAME_NONZERO_PATTERN), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_vadd_mvec)(const _ST_ alpha[], TYPE(const_mvec_ptr) vV,
                                    _ST_ beta,  TYPE(mvec_ptr)       vW,
                                    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);

  // no real support for this in petsc, scale columns one by one
  phist_lidx nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vV, &nvec, iflag), *iflag);
  TYPE(mvec_ptr) colV = NULL, colW = NULL;
  for(phist_lidx i = 0; i < nvec; i++)
  {
    PHIST_CHK_IERR(SUBR(mvec_view_block)((TYPE(mvec_ptr))vV, &colV, i, i, iflag), *iflag);
    PHIST_CHK_IERR(SUBR(mvec_view_block)((TYPE(mvec_ptr))vW, &colW, i, i, iflag), *iflag);
    PHIST_CHK_IERR(SUBR(mvec_add_mvec)(alpha[i], colV, beta, colW, iflag), *iflag);
  }
  PHIST_CHK_IERR(SUBR(mvec_delete)(colV, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(colW, iflag), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA,
                                      _ST_ beta,  TYPE(sdMat_ptr)       vB, 
                                      int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,B,vB,*iflag);

  PHIST_CHK_IERR( *iflag = MatScale(B->m, beta), *iflag);
  PHIST_CHK_IERR( *iflag = MatAXPY(B->m, alpha, A->m, SAME_NONZERO_PATTERN), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMatT_add_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA,
    _ST_ beta,  TYPE(sdMat_ptr)       vB, 
    int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,B,vB,*iflag);

  PHIST_CHK_IERR( *iflag = MatScale(B->m, beta), *iflag);
  Mat AT;
#ifdef IS_COMPLEX
  PHIST_CHK_IERR( *iflag = MatHermitianTranspose(A->m, MAT_INITIAL_MATRIX, &AT), *iflag);
#else
  PHIST_CHK_IERR( *iflag = MatTranspose(A->m, MAT_INITIAL_MATRIX, &AT), *iflag);
#endif
  PHIST_CHK_IERR( *iflag = MatAXPY(B->m, alpha, AT, SAME_NONZERO_PATTERN), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMat_times_mvec_communicate)(TYPE(const_sparseMat_ptr) vA, TYPE(const_mvec_ptr) vx, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  *iflag = 0;
}

extern "C" void SUBR(sparseMat_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) vA, 
    TYPE(const_mvec_ptr) vV, _ST_ beta, TYPE(mvec_ptr) vW, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,W,vW,*iflag);

  Mat tmpW;
  PHIST_CHK_IERR( *iflag = MatMatMult(A->m, V->v, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &tmpW), *iflag);
  PHIST_CHK_IERR( *iflag = MatScale(W->v, beta), *iflag);
  PHIST_CHK_IERR( *iflag = MatAXPY(W->v, alpha, tmpW, SAME_NONZERO_PATTERN), *iflag);
  PHIST_CHK_IERR( *iflag = MatDestroy(&tmpW), *iflag);

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sparseMatT_times_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) vA, 
    TYPE(const_mvec_ptr) vV, _ST_ beta, TYPE(mvec_ptr) vW, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sparseMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,W,vW,*iflag);

  Mat tmpW;
#ifdef IS_COMPLEX
  PHIST_CHK_IERR( *iflag = MatConjugate(V->v), *iflag);
#endif
  PHIST_CHK_IERR( *iflag = MatTransposeMatMult(A->m, V->v, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &tmpW), *iflag);
#ifdef IS_COMPLEX
  PHIST_CHK_IERR( *iflag = MatConjugate(V->v), *iflag);
#endif
  PHIST_CHK_IERR( *iflag = MatScale(W->v, beta), *iflag);
  PHIST_CHK_IERR( *iflag = MatAXPY(W->v, alpha, tmpW, SAME_NONZERO_PATTERN), *iflag);
  PHIST_CHK_IERR( *iflag = MatDestroy(&tmpW), *iflag);

  *iflag = PHIST_SUCCESS;
}

//! y[i]=alpha*(A*x[i]+shifts[i]*x[i]) + beta*y[i]
extern "C" void SUBR(sparseMat_times_mvec_vadd_mvec)(_ST_ alpha, TYPE(const_sparseMat_ptr) A,
        const _ST_ shifts[], TYPE(const_mvec_ptr) x, _ST_ beta, TYPE(mvec_ptr) y, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"

  PHIST_CHK_IERR(SUBR(sparseMat_times_mvec)(alpha,A,x,beta,y,iflag),*iflag);
  phist_lidx nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(x,&nvec,iflag),*iflag);
  {
    _ST_ alphashifts[nvec];
    for(phist_lidx i = 0; i < nvec; i++)
      alphashifts[i] = alpha*shifts[i];
    PHIST_CHK_IERR(SUBR(mvec_vadd_mvec)(alphashifts,x,st::one(),y,iflag),*iflag);
  }

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_dot_mvec)(TYPE(const_mvec_ptr) vV, 
                                    TYPE(const_mvec_ptr) vW, 
                                    _ST_* s, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,W,vW,*iflag);

  // quite a workaround...
  phist_lidx nvec;
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vV, &nvec, iflag), *iflag);
  Mat m;
  TYPE(mvec_ptr) vcolV = NULL, vcolW = NULL;
  for(phist_lidx i = 0; i < nvec; i++)
  {
    PHIST_CHK_IERR(SUBR(mvec_view_block)((TYPE(mvec_ptr))vV, &vcolV, i, i, iflag), *iflag);
    PHIST_CHK_IERR(SUBR(mvec_view_block)((TYPE(mvec_ptr))vW, &vcolW, i, i, iflag), *iflag);
    PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,colV,vcolV,*iflag);
    PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,colW,vcolW,*iflag);
#ifdef IS_COMPLEX
    PHIST_CHK_IERR( *iflag = MatConjugate(colV->v), *iflag);
#endif
    PHIST_CHK_IERR( *iflag = MatTransposeMatMult(colV->v, colW->v, (i==0) ? MAT_INITIAL_MATRIX : MAT_REUSE_MATRIX, PETSC_DEFAULT, &m), *iflag);
#ifdef IS_COMPLEX
    PHIST_CHK_IERR( *iflag = MatConjugate(colV->v), *iflag);
#endif
    PHIST_CHK_IERR( *iflag = MatGetTrace(m, &s[i]), *iflag);
  }
  if( nvec > 0 )
  {
    PHIST_CHK_IERR( *iflag = MatDestroy(&m), *iflag);
  }
  PHIST_CHK_IERR(SUBR(mvec_delete)(vcolV, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_delete)(vcolW, iflag), *iflag);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_times_sdMat)(_ST_ alpha, TYPE(const_mvec_ptr) vV, 
    TYPE(const_sdMat_ptr) vM, 
    _ST_ beta, TYPE(mvec_ptr) vW, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,M,vM,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,W,vW,*iflag);

  // not really easily available with petsc... lets just call BLAS
  phist_lidx nlocal, nvecv, nvecw;
  PHIST_CHK_IERR(phist_map_get_local_length(V->map,&nlocal,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vV, &nvecv, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(vW, &nvecw, iflag), *iflag);
  PHIST_TG_PREFIX(GEMM)("N", "N", &nlocal, &nvecw, &nvecv, &alpha, V->rawData, &nlocal, M->rawData, &M->lda, &beta, W->rawData, &nlocal);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA, 
    TYPE(const_sdMat_ptr) vB, 
    _ST_ beta, TYPE(sdMat_ptr) vC, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,B,vB,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,C,vC,*iflag);

  // not really easily available with petsc... lets just call BLAS
  phist_lidx m, n, k;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(vC, &m, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(vC, &n, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(vA, &k, iflag), *iflag);
  PHIST_TG_PREFIX(GEMM)("N", "N", &m, &n, &k, &alpha, A->rawData, &A->lda, B->rawData, &B->lda, &beta, C->rawData, &C->lda);
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMatT_times_sdMat)(_ST_ alpha, TYPE(const_sdMat_ptr) vA, 
    TYPE(const_sdMat_ptr) vB, 
    _ST_ beta, TYPE(sdMat_ptr) vC, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,B,vB,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,C,vC,*iflag);

  // not really easily available with petsc... lets just call BLAS
  phist_lidx m, n, k;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(vC, &m, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(vC, &n, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(vA, &k, iflag), *iflag);
#ifdef IS_COMPLEX
  PHIST_TG_PREFIX(GEMM)("H", "N", &m, &n, &k, &alpha, A->rawData, &A->lda, B->rawData, &B->lda, &beta, C->rawData, &C->lda);
#else
  PHIST_TG_PREFIX(GEMM)("T", "N", &m, &n, &k, &alpha, A->rawData, &A->lda, B->rawData, &B->lda, &beta, C->rawData, &C->lda);
#endif
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(sdMat_times_sdMatT)(_ST_ alpha, TYPE(const_sdMat_ptr) vA,
    TYPE(const_sdMat_ptr) vB,
    _ST_ beta, TYPE(sdMat_ptr) vC, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,A,vA,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::sdMat_t,B,vB,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,C,vC,*iflag);

  // not really easily available with petsc... lets just call BLAS
  phist_lidx m, n, k;
  PHIST_CHK_IERR(SUBR(sdMat_get_nrows)(vC, &m, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(vC, &n, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sdMat_get_ncols)(vA, &k, iflag), *iflag);
#ifdef IS_COMPLEX
  PHIST_TG_PREFIX(GEMM)("N", "H", &m, &n, &k, &alpha, A->rawData, &A->lda, B->rawData, &B->lda, &beta, C->rawData, &C->lda);
#else
  PHIST_TG_PREFIX(GEMM)("N", "T", &m, &n, &k, &alpha, A->rawData, &A->lda, B->rawData, &B->lda, &beta, C->rawData, &C->lda);
#endif
  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvecT_times_mvec)(_ST_ alpha, TYPE(const_mvec_ptr) vV, 
                                       TYPE(const_mvec_ptr) vW, 
                                       _ST_ beta, TYPE(sdMat_ptr) vM, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(const Traits<_ST_>::mvec_t,W,vW,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,M,vM,*iflag);

#ifdef IS_COMPLEX
    PHIST_CHK_IERR( *iflag = MatConjugate(V->v, &colVT), *iflag);
#endif
    Mat tmpM;
    PHIST_CHK_IERR( *iflag = MatTransposeMatMult(V->v, W->v, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &tmpM), *iflag);
    PHIST_CHK_IERR( *iflag = MatScale(M->m, beta), *iflag);
    // we obtained a parallel distributed matrix in tmpM, now "allgather" it
    // (required even on one process to convert it to a seqdense matrix
    Mat *tmpM_ = NULL;
    IS allRows, allCols;
    phist_lidx nrows, ncols;
    PHIST_CHK_IERR( *iflag = MatGetSize(tmpM, &nrows, &ncols), *iflag);
    PHIST_CHK_IERR( *iflag = ISCreateStride(*(MPI_Comm*)M->comm, nrows, 0, 1, &allRows), *iflag);
    PHIST_CHK_IERR( *iflag = ISCreateStride(*(MPI_Comm*)M->comm, ncols, 0, 1, &allCols), *iflag);
    PHIST_CHK_IERR( *iflag = MatGetSubMatrices(tmpM, 1, &allRows, &allCols, MAT_INITIAL_MATRIX, &tmpM_), *iflag);
    PHIST_CHK_IERR( *iflag = MatAXPY(M->m,alpha,tmpM_[0],SAME_NONZERO_PATTERN), *iflag);
    PHIST_CHK_IERR( *iflag = MatDestroyMatrices(1, &tmpM_), *iflag);
    PHIST_CHK_IERR( *iflag = ISDestroy(&allRows), *iflag);
    PHIST_CHK_IERR( *iflag = ISDestroy(&allCols), *iflag);
    PHIST_CHK_IERR( *iflag = MatDestroy(&tmpM), *iflag);
#ifdef IS_COMPLEX
    PHIST_CHK_IERR( *iflag = MatConjugate(V->v, &colVT), *iflag);
#endif

  *iflag = PHIST_SUCCESS;
}

extern "C" void SUBR(mvec_QR)(TYPE(mvec_ptr) vV, TYPE(sdMat_ptr) vR, int* iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
#include "phist_std_typedefs.hpp"
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::mvec_t,V,vV,*iflag);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sdMat_t,R,vR,*iflag);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_gather_mvecs)(TYPE(mvec_ptr) V, TYPE(const_mvec_ptr) W[], int nblocks, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_MARK_AS_EXPERIMENTAL(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(mvec_scatter_mvecs)(TYPE(const_mvec_ptr) V, TYPE(mvec_ptr) W[], int nblocks, int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_MARK_AS_EXPERIMENTAL(__FUNCTION__);
  *iflag=PHIST_NOT_IMPLEMENTED;
}

extern "C" void SUBR(sparseMat_create_fromRowFuncAndMap)(TYPE(sparseMat_ptr) *vA,
        phist_const_map_ptr map,
        phist_lidx maxnne,phist_sparseMat_rowFunc rowFunPtr,void* last_arg,
        int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CAST_PTR_FROM_VOID(Traits<_ST_>::sparseMat_t*,A,vA,*iflag);

  phist_lidx nlocal;
  phist_gidx ilower, nglob;
  phist_const_comm_ptr comm;
  PHIST_CHK_IERR(phist_map_get_local_length(map, &nlocal, iflag), *iflag);
  PHIST_CHK_IERR(phist_map_get_ilower(map, &ilower, iflag), *iflag);
  PHIST_CHK_IERR(phist_map_get_global_length(map, &nglob, iflag), *iflag);
  PHIST_CHK_IERR(phist_map_get_comm(map, &comm, iflag), *iflag);

  PHIST_CHK_IERR( *iflag = PetscNew(A), *iflag);
  (*A)->map = map;
  //PHIST_CHK_IERR( *iflag = MatCreateAIJ(*(MPI_Comm*)comm, nlocal, nlocal, nglob, nglob, maxnne, NULL, maxnne, NULL, &((*A)->m)), *iflag);
  PHIST_CHK_IERR( *iflag = MatCreate(*(MPI_Comm*)comm, &((*A)->m)), *iflag);
  PHIST_CHK_IERR( *iflag = MatSetType((*A)->m, MATAIJ), *iflag);
  PHIST_CHK_IERR( *iflag = MatSetSizes((*A)->m, nlocal, nlocal, nglob, nglob), *iflag);
  MatType matType;
  PHIST_CHK_IERR( *iflag = MatGetType((*A)->m, &matType), *iflag);
  if( std::string(matType) == std::string(MATMPIAIJ) )
  {
    PHIST_CHK_IERR( *iflag = MatMPIAIJSetPreallocation((*A)->m, maxnne, NULL, maxnne, NULL), *iflag);
  }
  else if( std::string(matType) == std::string(MATSEQAIJ) )
  {
    PHIST_CHK_IERR( *iflag = MatSeqAIJSetPreallocation((*A)->m, maxnne, NULL), *iflag);
  }
  else
  {
    PHIST_SOUT(PHIST_ERROR, "strange PETSc matrix type %s\n", matType);
    *iflag = PHIST_NOT_IMPLEMENTED;
    return;
  }

  _ST_ *vals = NULL;
  phist_gidx *cols = NULL;
  PHIST_CHK_IERR( *iflag = PetscMalloc1(maxnne, &vals), *iflag);
  PHIST_CHK_IERR( *iflag = PetscMalloc1(maxnne, &cols), *iflag);

  for (phist_lidx i = 0; i < nlocal; i++)
  {
    phist_gidx row = ilower+i;
    phist_lidx row_nnz = 0;
    PHIST_CHK_IERR( *iflag = rowFunPtr(row,&row_nnz,cols,vals,last_arg), *iflag);
    PHIST_CHK_IERR( *iflag = MatSetValues((*A)->m,1,&row,row_nnz,cols,vals,INSERT_VALUES), *iflag);
  }
  PHIST_CHK_IERR( *iflag = PetscFree(vals), *iflag);
  PHIST_CHK_IERR( *iflag = PetscFree(cols), *iflag);
  PHIST_CHK_IERR( *iflag = MatAssemblyBegin((*A)->m,MAT_FINAL_ASSEMBLY), *iflag);
  PHIST_CHK_IERR( *iflag = MatAssemblyEnd((*A)->m,MAT_FINAL_ASSEMBLY), *iflag);

  *iflag = PHIST_SUCCESS;
}


// NOTE: see the description of sparseMat_read_mm on how we treat input flags for this function
extern "C" void SUBR(sparseMat_create_fromRowFunc)(TYPE(sparseMat_ptr) *A, phist_const_comm_ptr vcomm,
        phist_gidx nrows, phist_gidx ncols, phist_lidx maxnne,
                phist_sparseMat_rowFunc rowFunPtr, void* last_arg,
                int *iflag)
{
  PHIST_ENTER_KERNEL_FCN(__FUNCTION__);
  PHIST_CHK_IERR( *iflag = (nrows == ncols) ? PHIST_SUCCESS : PHIST_NOT_IMPLEMENTED, *iflag);

  phist_map_ptr map = NULL;
  PHIST_CHK_IERR(phist_map_create(&map, vcomm, nrows, iflag), *iflag);
  PHIST_CHK_IERR(SUBR(sparseMat_create_fromRowFuncAndMap)(A, map, maxnne, rowFunPtr, last_arg, iflag), *iflag);

  *iflag = PHIST_SUCCESS;
}

