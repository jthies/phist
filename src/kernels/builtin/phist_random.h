/*! \file random.h
 * a variant of George Marsaglia's KISS random number generator with jump ahead
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 * inspired by https://github.com/cmcqueen/simplerandom and http://web.mst.edu/~vojtat/class_5403/kiss05/rkiss05.f90
 *
*/

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <stdint.h>


//! initialize the random number generator
void phist_random_init();

//! get a set of random numbers between -1 and +1
void phist_Drandom_number(int n, double*x);

// generate random numbers in parallel for a dense block
void drandom_1(int nrows, double *restrict y, int64_t pre_skip, int64_t post_skip);

// generate random numbers in parallel for strided block
void drandom_general(int nvec, int nrows, double *restrict v, int ldv, int64_t pre_skip, int64_t post_skip);
