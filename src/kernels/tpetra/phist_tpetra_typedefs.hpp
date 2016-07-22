#ifndef KERNELS_TPETRA_TYPEDEFS_HPP
#define KERNELS_TPETRA_TYPEDEFS_HPP

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "phist_typedefs.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DataAccess.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "KokkosClassic_config.h"
#ifdef HAVE_KOKKOSCLASSIC_OPENMP
#include "Kokkos_OpenMPNode.hpp"
#endif
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
#include "Kokkos_TPINode.hpp"
#endif
#ifdef HAVE_KOKKOSCLASSIC_TBB
#include "Kokkos_TBBNode.hpp"
#endif
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_DefaultKernels.hpp"

#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
//#include "Tpetra_CrsMatrixMultiplyOp.hpp"

#include "phist_typedefs.h"

namespace phist {
namespace tpetra {

// notte: with the OpenMP node and Trilinos 11.12, TSQR doesn't work
//#ifdef HAVE_KOKKOSCLASSIC_OPENMP
//typedef KokkosClassic::OpenMPNode node_type; // from the Kokkos node API
#if defined(HAVE_KOKKOSCLASSIC_TBB)
typedef KokkosClassic::TBBNode node_type; // from the Kokkos node API
#elif defined(HAVE_KOKKOSCLASSIC_THREADPOOL)
typedef KokkosClassic::TPINode node_type; // from the Kokkos node API
#else
typedef Kokkos::DefaultNode::DefaultNodeType node_type; // from the Kokkos node API
#endif
typedef Tpetra::Map<phist_lidx,phist_gidx,node_type> map_type;
typedef Tpetra::Import<phist_lidx,phist_gidx,node_type> import_type;
typedef Tpetra::Export<phist_lidx,phist_gidx,node_type> export_type;
typedef Teuchos::Comm<int> comm_type;

template<typename ST>
class Traits
{
  public:
  
  //! multi vectors
  typedef Tpetra::MultiVector<ST,phist_lidx,phist_gidx,node_type> mvec_t;

  //! serial dense matrix - just a multivector with a serial map.
  typedef Tpetra::MultiVector<ST,phist_lidx,phist_gidx,node_type> sdMat_t;

  //! serial dense matrix from Teuchos, we need this for e.g. the BLAS interface.
  typedef Teuchos::SerialDenseMatrix<phist_lidx,ST> Teuchos_sdMat_t;

  //! CRS matrices
  typedef Tpetra::CrsMatrix<ST,phist_lidx,phist_gidx,node_type> sparseMat_t;

  //! for performing the MVM
//  typedef Tpetra::CrsMatrixMultiplyOp<ST,ST,phist_lidx,phist_gidx,node_type> crsMVM_t;

  //! scalar 1
  static inline ST one(){return Teuchos::ScalarTraits<ST>::one();}

  //! scalar 0
  static inline ST zero(){return Teuchos::ScalarTraits<ST>::zero();}

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

  
};

}//namespace tpetra
}//namespace phist

#endif
