/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  ...
 *  \ingroup internal
 */

#ifndef SCAMAC_OMP_H
#define SCAMAC_OMP_H

#include "scamac_config.h"

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#include "time.h"
#define omp_get_wtime() ( (double) clock() / (double) CLOCKS_PER_SEC )
#endif

#endif /* SCAMAC_OMP_H */
