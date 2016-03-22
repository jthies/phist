#ifndef PHIST_EPETRA_HG_REPART_HPP
#define PHIST_EPETRA_HG_REPART_HPP

#include "phist_config.h"

#ifndef DOXYGEN
#include <mpi.h>
#include "Teuchos_RCP.hpp"
#endif

class Epetra_CrsMatrix;
class Epetra_RowMatrix;
class Epetra_Map;


namespace phist 
{
  namespace epetra_internal 
  {

    // on input, rowmatrix must be Filled(). on output, bal_map points to the new (repartitioned) map,
    // and bal_mat contains the repartitioned matrix *if redist==true*. Otherwise, bal_mat=Teuchos::null
    // on output.
    void repartition(Teuchos::RCP<const Epetra_RowMatrix> rowmatrix,
                     Teuchos::RCP<      Epetra_Map>& bal_map,
                     Teuchos::RCP<      Epetra_CrsMatrix>& bal_mat, bool redist,
                     int *iflag);
                                                   

  }
}

#endif
