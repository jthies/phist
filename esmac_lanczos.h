#ifndef ESMAC_LANCZOS_H
#define ESMAC_LANCZOS_H

#include "esmac_sparsemat.h"
#include "esmac_generator.h"

int esmac_lanczos_ev_mat(const esmac_sparsemat_t *sm, double tol, double *ev1, double *ev2, double *eps1, double *eps2);

int esmac_lanczos_ev(const esmac_generator_t *gen, double tol, double *ev1, double *ev2, double *eps1, double *eps2);


#endif /* ESMAC_LANCZOS_H */
