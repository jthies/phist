#ifndef ESMAC_SPECTRUM_H
#define ESMAC_SPECTRUM_H

#include "esmac_generator.h"
#include "esmac_sparsemat.h"

// Compute entire spectrum of (real, symmetric) matrix. Use only for *really* small matrices
int esmac_spectrum(const esmac_generator_t * gen, double **spec);
int esmac_spectrum_mat(const esmac_sparsemat_t * sm, double **spec);

#endif /* ESMAC_SPECTRUM_H */
