#ifndef ESMAC_PLOT_H
#define ESMAC_PLOT_H

#include "esmac_generator.h"

// plot sparsity pattern "one pixel per entry". Use only for small matrices
int esmac_plot_pattern_onetoone(const esmac_generator_t * gen, const char * filename);

// plot sparsity pattern with npx x npx pixels
int esmac_plot_pattern_from_gen(const esmac_generator_t * gen, int downscale, const char * filename);


int esmac_plot_pattern(int pattern_px, int pattern_py, const int * pattern, int downscale, const char * filename);


#endif /* ESMAC_PLOT_H */
