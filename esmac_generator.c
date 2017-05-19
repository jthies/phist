#include <stdlib.h>
#include <string.h>
#include "esmac_constants.h"
#include "esmac_generator.h"

esmac_generator_work_t * esmac_generator_alloc(const esmac_generator_t * gen, int * info) {
  esmac_generator_work_t * ws;
  ws = malloc(sizeof * ws);
  if (gen->gen_alloc) {
    ws->ws = gen -> gen_alloc(gen->par, info);
  } else {
    ws->ws = NULL;
  }
  if (gen->gen_set_info) {
    gen -> gen_set_info(gen->par, ws->ws, &(ws->info));
  }
  ws -> work_free = gen->gen_free;
  if (info) {
    *info=0;
  }  
  return ws;
}

void esmac_generator_free(esmac_generator_work_t * ws) {
  if (ws) {
    if (ws->ws && ws->work_free) {
      ws->work_free(ws->ws);
    }
    free(ws);
  }
}

esmac_idx_t esmac_generator_query(const esmac_generator_t * gen, const esmac_generator_work_t * ws, const char * query) {
  if (!strcmp(query,"nrow")) {
    return ws->info.nrow;
  } else if (!strcmp(query,"ncol")) {
    return ws->info.nrow;
  } else if (!strcmp(query,"maxnzrow")) {
    return ws->info.maxrowlength;
  } else if (!strcmp(query,"valtype")) {
    return ws->info.valtype;
  } else {
    return -1;
  }
}

/*
//esmac_idx_t esmac_generator_ns(const esmac_generator_t * gen, const esmac_generator_work_t * ws) {
//  return gen->gen_ns(ws->ws);
//}
esmac_idx_t esmac_generator_nrow(const esmac_generator_t * gen, const esmac_generator_work_t * ws) {
  return ws->info.nrow;
}
esmac_idx_t esmac_generator_ncol(const esmac_generator_t * gen, const esmac_generator_work_t * ws) {
  return ws->info.ncol;
}
int esmac_generator_maxnzr(const esmac_generator_t * gen, const esmac_generator_work_t * ws) {
  return ws->info.maxrowlength;
}
*/

int esmac_generator_row(const esmac_generator_t * gen, esmac_generator_work_t * ws, esmac_idx_t irow, esmac_idx_t *cind, double *val) {
  return gen->gen_row(gen->par, ws->ws, irow, cind, val);
}

