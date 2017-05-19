#ifndef ESMAC_GENERATOR_H
#define ESMAC_GENERATOR_H

#include "esmac_constants.h"
#include "esmac_types.h"


typedef void * (*esmac_gen_alloc_t)(const void * par, int * info);
typedef void  (*esmac_gen_free_t)(void * ws);
//typedef esmac_idx_t (*esmac_gen_ns_t)(void * ws);
//typedef int (*esmac_gen_maxrowlength_t)(void * ws);
typedef int (*esmac_gen_set_info_t)(const void * par, void *ws, esmac_matrix_info_t *info);
typedef int (*esmac_gen_row_t)(const void * par, void * ws, esmac_idx_t irow, esmac_idx_t *cind, double *val);

typedef struct {
  char name[ESMAC_NAME_LENGTH];

  void * par;
 
  esmac_gen_alloc_t gen_alloc;
  esmac_gen_free_t gen_free;
 // esmac_gen_ns_t gen_ns;
 // esmac_gen_maxrowlength_t gen_maxrowlength;
  esmac_gen_set_info_t gen_set_info;
  esmac_gen_row_t gen_row;

} esmac_generator_t;

typedef struct {
  void *ws;
  esmac_matrix_info_t info;
  esmac_gen_free_t work_free;
} esmac_generator_work_t;


esmac_generator_work_t * esmac_generator_alloc(const esmac_generator_t * gen, int * info) ;
void esmac_generator_free(esmac_generator_work_t * ws);

esmac_idx_t esmac_generator_query(const esmac_generator_t * gen, const esmac_generator_work_t * ws, const char * query);

/*
esmac_idx_t esmac_generator_nrow(const esmac_generator_t * gen, const esmac_generator_work_t * ws);
esmac_idx_t esmac_generator_ncol(const esmac_generator_t * gen, const esmac_generator_work_t * ws);
int esmac_generator_maxnzr(const esmac_generator_t * gen, const esmac_generator_work_t * ws);
*/

int esmac_generator_row(const esmac_generator_t * gen, esmac_generator_work_t * ws, esmac_idx_t irow, esmac_idx_t *cind, double *val);


#endif /* ESMAC_GENERATOR_H */
