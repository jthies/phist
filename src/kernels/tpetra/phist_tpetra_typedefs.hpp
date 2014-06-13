#ifndef KERNELS_TPETRA_TYPEDEFS_HPP
#define KERNELS_TPETRA_TYPEDEFS_HPP

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DataAccess.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_DefaultKernels.hpp"

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_CrsMatrixMultiplyOp.hpp"

#include "phist_typedefs.h"

namespace phist {
namespace tpetra {

// currently typedef'd in phist_typedefs.hpp (TODO)
//typedef Kokkos::DefaultNode::DefaultNodeType node_t; // from the Kokkos node API
typedef Tpetra::Map<lidx_t,gidx_t,node_t> map_t;
typedef Tpetra::Import<lidx_t,gidx_t,node_t> import_t;
typedef Tpetra::Export<lidx_t,gidx_t,node_t> export_t;
typedef Teuchos::Comm<int> comm_t;

template<typename ST>
class Traits
  {
  public:
  
  //!
  typedef typename Kokkos::DefaultKernels<ST,lidx_t,node_t>::SparseOps localOps_t;

  //! multi vectors
  typedef Tpetra::MultiVector<ST,lidx_t,gidx_t,node_t> mvec_t;

  //! serial dense matrix - just a multivector with a serial map.
  typedef Tpetra::MultiVector<ST,lidx_t,gidx_t,node_t> sdMat_t;

  //! serial dense matrix from Teuchos, we need this for e.g. the BLAS interface.
  typedef Teuchos::SerialDenseMatrix<lidx_t,ST> Teuchos_sdMat_t;

  //! CRS matrices
  typedef Tpetra::CrsMatrix<ST,lidx_t,gidx_t,node_t,localOps_t> crsMat_t;

  //! for performing the MVM
  typedef Tpetra::CrsMatrixMultiplyOp<ST,ST,lidx_t,gidx_t,node_t,localOps_t> crsMVM_t;

  //! scalar 1
  static inline ST one(){return Teuchos::ScalarTraits<ST>::one();}

  //! scalar 0
  static inline ST zero(){return Teuchos::ScalarTraits<ST>::zero();}

  //! create a Teuchos' view of a local mvec/sdMat
  static Teuchos::RCP<const Teuchos_sdMat_t> CreateTeuchosView(Teuchos::RCP<const sdMat_t> M, int* ierr)
    {
    *ierr=0;
    lidx_t stride = M->getStride();
    lidx_t nrows = M->getLocalLength();
    lidx_t ncols = M->getNumVectors();
    
    Teuchos::ArrayRCP<const ST> M_tmp;
    bool status=true;
    try {
    M_tmp=M->get1dView();
    } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,status);
    if (!status) {*ierr=-1; return Teuchos::null;}
    const ST *M_val = M_tmp.getRawPtr();
    Teuchos::RCP<const Teuchos_sdMat_t> M_view
                  = Teuchos::rcp(new Teuchos_sdMat_t(Teuchos::View,M_val,stride,nrows,ncols));
    return M_view;     
    }

  //! create a non-const Teuchos' view of a local mvec/sdMat
  static Teuchos::RCP<Teuchos_sdMat_t> CreateTeuchosViewNonConst(Teuchos::RCP<sdMat_t> M, int* ierr)
    {
    *ierr=0;
    int stride = M->getStride();
    int nrows = M->getLocalLength();
    int ncols = M->getNumVectors();
    Teuchos::ArrayRCP<ST> M_tmp=Teuchos::null;
    bool status=true;
    try {
    M_tmp=M->get1dViewNonConst();
    } TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,status);
    if (!status) {*ierr=-1; return Teuchos::null;}
    ST *M_val = M_tmp.getRawPtr();
    Teuchos::RCP<Teuchos_sdMat_t> M_view
                  = Teuchos::rcp(new Teuchos_sdMat_t(Teuchos::View,M_val,stride,nrows,ncols));
    return M_view;                  
    }

  
  };

}//namespace tpetra
}//namespace phist

#endif
