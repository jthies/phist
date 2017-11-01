/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  ...
 *  \ingroup internal
 */

#ifndef SCAMAC_INTERNAL_H
#define SCAMAC_INTERNAL_H

#include "scamac_include.h"
#include "scamac_sparserow.h"

typedef ScamacErrorCode (*scamac_gen_work_alloc_ft)(const void * par, const void * tab, void ** ws);
typedef ScamacErrorCode (*scamac_gen_work_free_ft)(void * ws);
typedef ScamacErrorCode (*scamac_gen_tables_create_ft)(const void * par, void ** tables, ScamacInfo * info);
typedef ScamacErrorCode (*scamac_gen_tables_destroy_ft)(void * tables);
typedef ScamacErrorCode (*scamac_gen_row_ft)(const void * par, const void * table, void * ws, ScamacIdx irow, ScamacFlag flag, void *row);
typedef ScamacErrorCode (*scamac_gen_check_ft)(const void * par, char ** desc);

typedef ScamacErrorCode (*scamac_gen_unwrap_ft)(const void * wrapped_par, void *par);

struct scamac_info_st {
  ScamacIdx nrow;
  ScamacIdx ncol;
  ScamacIdx maxnzrow;
  ScamacIdx maxnzcol;
  ScamacIdx maxnz;
  int valtype; // SCAMAC_NONE or SCAMAC_VAL_REAL or SCAMAC_VAL_COMPLEX
  int symmetry; // SCAMAC_NONE or SCAMAC_SYMMETRIC or SCAMAC_HERMITIAN
};

struct scamac_generator_st {
  char name[SCAMAC_NAME_LENGTH];

  void * par;

  scamac_gen_work_alloc_ft fct_work_alloc;
  scamac_gen_work_free_ft fct_work_free;
  scamac_gen_tables_create_ft fct_tables_create;
  scamac_gen_tables_destroy_ft fct_tables_destroy;
  scamac_gen_row_ft fct_gen_row;
  scamac_gen_check_ft fct_check;

  // for "wrappers"
  bool needs_unwrapping;
  void * wrapped_par;
  scamac_gen_unwrap_ft fct_unwrap_par;
  scamac_gen_check_ft fct_check_wrapped;

  // fully set after "finalize"
  bool needs_finalization;
  void * tables;
  ScamacInfo info;

};

struct scamac_workspace_st {
  scamac_sparserow_real_st * row_real;
  scamac_sparserow_cplx_st * row_cplx;
  scamac_gen_work_free_ft fct_work_free;
  void * ws;
};

// macro for scamac_matrix_***_check routines
#define SCAMAC_DESC_ERR(cond,conddesc)  \
    do { \
      if (cond) { \
        err=SCAMAC_EINVAL; \
        if (desc) {scamac_string_append(&str,conddesc"\n");} \
      } \
    } while (0)
#define SCAMAC_DESC_WARN(cond,conddesc)  \
    do { \
      if (cond) { \
        err=SCAMAC_EWARNING; \
        if (desc) {scamac_string_append(&str,"[WARNING] "conddesc"\n");} \
      } \
    } while (0)

#endif /* SCAMAC_INTERNAL_H */
