/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  Write sparse matrices to file
 *  \ingroup toolbox
 */

#ifndef SCAMAC_SPARSEMAT_IO_H
#define SCAMAC_SPARSEMAT_IO_H

#include "scamac_include.h"
#include "scamac_sparsemat.h"

/* output in Matrix Market format */
ScamacErrorCode scamac_sparsemat_io_write_mm(const scamac_sparsemat_st *sm, char * fname);

/* output in Harwell-Boeing format */
ScamacErrorCode scamac_sparsemat_io_write_hb(const scamac_sparsemat_st *sm, char * fname);

/* output in MATLAB binary format */
ScamacErrorCode scamac_sparsemat_io_write_matlab(const scamac_sparsemat_st *sm, char * fname);

/* output in GHOST binary format */
ScamacErrorCode scamac_sparsemat_io_write_ghost(const scamac_sparsemat_st *sm, char * fname);


#endif /* SCAMAC_SPARSEMAT_IO_H */
