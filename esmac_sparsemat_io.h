#ifndef ESMAC_SPARSEMAT_IO_H
#define ESMAC_SPARSEMAT_IO_H

#include "esmac_sparsemat.h"

/* output in Matrix Market format */
int esmac_sparsemat_io_write_mm(const esmac_sparsemat_t *sm, char * fname);

/* output in Harwell-Boeing format */
int esmac_sparsemat_io_write_hb(const esmac_sparsemat_t *sm, char * fname);

/* output in MATLAB binary format */
int esmac_sparsemat_io_write_matlab(const esmac_sparsemat_t *sm, char * fname);

/* output in GHOST binary format */
int esmac_sparsemat_io_write_ghost(const esmac_sparsemat_t *sm, char * fname);


#endif /* ESMAC_SPARSEMAT_IO_H */
