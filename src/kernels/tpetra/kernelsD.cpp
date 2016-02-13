#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_macros.h"
#include "phist_typedefs.h"
#include "phist_tpetra_typedefs.hpp"
#include "../phist_kernels.h"
#include "phist_trilinos_macros.h"
#include "phist_ScalarTraits.hpp"

#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_DefaultComm.hpp"
#ifdef PHIST_HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_OpaqueWrapper.hpp"
#endif
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "Tpetra_DefaultPlatform.hpp"

#include "./Tpetra_TsqrAdaptor.hpp"
#include "./BelosTpetraAdapter.hpp"
#include "BelosTsqrOrthoManager.hpp"

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

#include <fstream>
#include <sstream>


using namespace phist::tpetra;

extern "C" void phist_tpetra_node_create(node_t** node, const_comm_ptr_t vcomm, int* iflag);

#include "phist_gen_d.h"
#include "kernels_def.hpp"
#include "carp_def.hpp"
#include "../common/kernels_no_inplace_VC.cpp"
#include "../common/kernels_no_fused.cpp"
#include "../common/kernels_no_gpu.cpp"

