/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#ifndef PHIST_BENCH_KERNELS_H
#define PHIST_BENCH_KERNELS_H

#ifdef __cplusplus
extern "C" {
#endif

void phist_bench_stream_load(double* mean_bw, double* max_bw, int* iflag);

void phist_bench_stream_store(double* mean_bw, double* max_bw, int* iflag);

void phist_bench_stream_triad(double* mean_bw, double* max_bw, int* iflag);


#ifdef __cplusplus
} //extern "C"
#endif

#endif
