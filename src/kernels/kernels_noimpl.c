// these implementations can be used if a kernel package
// does not implement all four data types (cf. epetra/ for
// an example, which supports only double (D))

// \name Matrix input from a file

//@{

//! read a matrix from a MatrixMarket (ASCII) file
_SUBROUTINE_(read_crsMat_mm)(void** A, char* filename,int* ierr)
  {
  *ierr=-99;
  }

//! read a matrix from a Ghost CRS (binary) file.
_SUBROUTINE_(read_crsMat_bin)(void** A, char* filename,int* ierr)
  {
  *ierr=-99;
  }

//! read a matrix from a Harwell-Boeing (HB) file
_SUBROUTINE_(read_crsMat_hb)(void** vA, char* filename,int* ierr)
  {
  *ierr=-99;
  }
//!@}

//! \name get information about the data distribution in a matrix (maps)

//!@{
//! get the row distribution of the matrix
_SUBROUTINE_(get_row_map)(void* vA, void** vmap, int* ierr)
  {
  *ierr=-99;
  }

//! get column distribution of a matrix
_SUBROUTINE_(get_col_map)(void* A, void** map, int* ierr)
  {
  *ierr=-99;
  }

//! get the map for vectors x in y=A*x
_SUBROUTINE_(get_domain_map)(void* A, void** map, int* ierr)
  {
  *ierr=-99;
  }

//! get the map for vectors y in y=A*x
_SUBROUTINE_(get_range_map)(void* A, void** map, int* ierr)
  {
  *ierr=-99;
  }
//@}

//! \name constructors

//@{
//! create a block-vector. The entries are stored contiguously
//! at val in column major ordering.
_SUBROUTINE_(create_mvec)(void* vmap, int nvec, void** vV, _ST_** val, int* ierr)
  {
  *ierr=-99;
  }

//! create a serial dense n x m matrix on all procs, with column major
//! ordering.
_SUBROUTINE_(create_sdMat)(int n, int m, void** M, _ST_** val, int* ierr)
  {
  *ierr=-99;
  }

//@}

//! \name destructors

//@{

//!
_SUBROUTINE_(delete_crsMat)(void* vA, int* ierr)
  {
  *ierr=-99;
  }

//!
_SUBROUTINE_(delete_mvec)(void* vV, int* ierr)
  {
  *ierr=-99;
  }

//!
_SUBROUTINE_(delete_sdMat)(void* vM, int* ierr)
  {
  *ierr=-99;
  }

//@}

//! \name Numerical functions
//!@{

//! y=alpha*A*x+beta*y.
_SUBROUTINE_(crsMat_X_mvec)(_ST_ alpha, void* vA, void* vx, _ST_ beta, void* vy, int* ierr)
  {
  *ierr=-99;
  }

//! dot product of vectors v_i and w_i, i=1..numvecs
_SUBROUTINE_(mvec_dot_mvec)(void* vv, void* vw, _ST_* s, int* ierr)
  {
  *ierr=-99;
  }

//! dense tall skinny matrix-matrix product yielding a serial dense matrix
//! C=V'*W. C is replicated on all MPI processes sharing V and W.
_SUBROUTINE_(mvecT_X_mvec)(void* vV, void* vW, void* vC, int* ierr)
  {
  *ierr=-99;
  }

//! 'tall skinny' QR decomposition, V=Q*R, Q'Q=I, R upper triangular.
_SUBROUTINE_(mvec_QR)(void* V, void* Q, void* R, int* ierr)
  {
  *ierr=-99;
  }

//!@}

