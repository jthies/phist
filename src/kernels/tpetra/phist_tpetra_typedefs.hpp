#ifndef KERNELS_TPETRA_TYPEDEFS_HPP
#define KERNELS_TPETRA_TYPEDEFS_HPP

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "phist_typedefs.h"

#include "Tpetra_Map_decl.hpp"
#include "Tpetra_MultiVector_decl.hpp"
#include "Tpetra_CrsMatrix_decl.hpp"
#include "Tpetra_Import_decl.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace phist
{
namespace tpetra
{
using map_type = Tpetra::Map<phist_lidx, phist_gidx>;
using sdMat_map_type = Tpetra::Map<phist_lidx, phist_gidx>;
using comm_type =  Teuchos::Comm<int>;
using import_type = Tpetra::Import<phist_lidx, phist_gidx>;

template<typename ST>
class Traits
{
  public:
  
  //! multi vectors
  using mvec_t = Tpetra::MultiVector<ST, phist_lidx, phist_gidx>;

  //! serial dense matrix - just a multivector with a serial map.
  using sdMat_t = Tpetra::MultiVector<ST, phist_lidx, phist_gidx>;

  //! serial dense matrix from Teuchos, we need this for e.g. the BLAS interface.
  using Teuchos_sdMat_t = Teuchos::SerialDenseMatrix<phist_lidx, ST>;

  //! CRS matrices
  using sparseMat_t = Tpetra::CrsMatrix<ST, phist_lidx, phist_gidx>;

  //! for performing the MVM
//  typedef Tpetra::CrsMatrixMultiplyOp<ST,ST,phist_lidx,phist_gidx,node_type> crsMVM_t;

  //! scalar 1
  static inline ST one()
  {
    return Teuchos::ScalarTraits<ST>::one();
  }

  //! scalar 0
  static inline ST zero()
  {
      return Teuchos::ScalarTraits<ST>::zero();
  }
/*
  //! create a Teuchos' view of a local mvec/sdMat
  static Teuchos::RCP<const Teuchos_sdMat_t> CreateTeuchosView(Teuchos::RCP<const sdMat_t> M, int* iflag)
  {
    *iflag=0;
    phist_lidx stride = M->getStride();
    phist_lidx nrows = M->getLocalLength();
    phist_lidx ncols = M->getNumVectors();
    
    Teuchos::ArrayRCP<const ST> M_tmp;
    bool status=true;
    try {
    M_tmp=M->get1dView();
    } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,status);
    if (!status) {*iflag=-1; return Teuchos::null;}
    const ST *M_val = M_tmp.getRawPtr();
    Teuchos::RCP<const Teuchos_sdMat_t> M_view
                  = Teuchos::rcp(new Teuchos_sdMat_t(Teuchos::View,M_val,stride,nrows,ncols));
    return M_view;     
  }

  //! create a non-const Teuchos' view of a local mvec/sdMat
  static Teuchos::RCP<Teuchos_sdMat_t> CreateTeuchosViewNonConst(Teuchos::RCP<sdMat_t> M, int* iflag)
  {
    *iflag=0;
    int stride = M->getStride();
    int nrows = M->getLocalLength();
    int ncols = M->getNumVectors();
    Teuchos::ArrayRCP<ST> M_tmp=Teuchos::null;
    bool status=true;
    try {
    M_tmp=M->get1dViewNonConst();
    } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,status);
    if (!status) {*iflag=-1; return Teuchos::null;}
    ST *M_val = M_tmp.getRawPtr();
    Teuchos::RCP<Teuchos_sdMat_t> M_view
                  = Teuchos::rcp(new Teuchos_sdMat_t(Teuchos::View,M_val,stride,nrows,ncols));
    return M_view;                  
  }
*/
  
};




} //namespace tpetra
} //namespace phist
#endif
