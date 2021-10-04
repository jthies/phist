/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
/*! \file stream_bench.h
 * header file for stream_bench.c
 * \author "Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
 *
*/
#ifndef PHIST_COMMON_STREAM_BENCH_H
#define PHIST_COMMON_STREAM_BENCH_H

#ifdef __cplusplus
extern "C" {
#endif

//! allocate memory for bench_stream_load_run
void dbench_stream_load_create(double** x, int* ierr);

//! purely load dominated micro benchmark (for large n), actually calculates a sum and determines the bandwidth
//! \warning trust the compiler for now to do useful things
void dbench_stream_load_run(const double* x, double* res, double* bw, int* ierr);

//! delete memory for bench_stream_load_run
void dbench_stream_load_destroy(double* x, int* ierr);


//! allocate memory for bench_stream_store_run
void dbench_stream_store_create(double** x, int* ierr);

//! purely store dominated micro benchmark (for large n), actually calculates some values and determines the bandwidth
void dbench_stream_store_run(double* x, const double* res, double* bw, int* ierr);

//! delete memory for bench_stream_store_run
void dbench_stream_store_destroy(double* x, int* ierr);


//! allocate memory for bench_stream_triad_run
void dbench_stream_triad_create(double** x, double** y, double** z, int* ierr);

//! stream triad micro benchmark, determines the bandwidth
void dbench_stream_triad_run(const double* x, const double* y, double* z, const double* res, double* bw, int* ierr);

//! delete memory for bench_stream_triad_run
void dbench_stream_triad_destroy(double* x, double* y, double* z, int* ierr);

#ifdef __cplusplus
} // extern "C"
#endif

#endif
