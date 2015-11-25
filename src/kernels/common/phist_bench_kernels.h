#ifndef PHIST_BENCH_KERNELS_H
#define PHIST_BENCH_KERNELS_H

#ifdef __cplusplus
extern "C" {
#endif

void phist_bench_stream_load(double* max_bw, int* iflag);

void phist_bench_stream_store(double* max_bw, int* iflag);

void phist_bench_stream_triad(double* max_bw, int* iflag);


#ifdef __cplusplus
} //extern "C"
#endif

#endif
