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
};
} //namespace tpetra
} //namespace phist
#endif
