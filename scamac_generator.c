#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "scamac_include.h"
#include "scamac_generator.h"

#include "scamac_internal.h"

ScamacErrorCode scamac_workspace_alloc(const ScamacGenerator * gen, ScamacWorkspace ** ws) {

  if (!gen) {
    return SCAMAC_ENULL | 1 << SCAMAC_ESHIFT;
  }
  if (!ws) {
    return SCAMAC_ENULL | 2 << SCAMAC_ESHIFT;
  }
  
  if (gen->needs_finalization) {
    return SCAMAC_EINVALID | 1 << SCAMAC_ESHIFT; // forgot to call "scamac_generator_finalize" !
  }

  ScamacErrorCode err;

  ScamacWorkspace * myws;
  myws = malloc(sizeof * myws);
  if (!myws) {
    return SCAMAC_EMALLOCFAIL;
  }

  if (gen->fct_work_alloc) {
    err = gen -> fct_work_alloc(gen->par, gen->tables, &(myws->ws) );
    if (err) {
      return (err | SCAMAC_EINTERNAL);
    }
    myws->fct_work_free = gen->fct_work_free;
  } else {
    myws->fct_work_free = NULL;
    myws->ws = NULL;
  }

  // allocate the respective sparse row
  if (gen->info.valtype == SCAMAC_VAL_REAL) {
    err = scamac_sparserow_real_alloc(&(myws->row_real));
    myws->row_cplx = NULL;
    if (err) {
      return  (err | SCAMAC_EINTERNAL);
    }
  } else if (gen->info.valtype == SCAMAC_VAL_COMPLEX) {
    myws->row_real = NULL;
    err = scamac_sparserow_cplx_alloc(&(myws->row_cplx));
    if (err) {
      return  (err | SCAMAC_EINTERNAL);
    }
  } else {
    return SCAMAC_EFAIL | SCAMAC_EINTERNAL;
  }

  *ws = myws;
  return SCAMAC_EOK;
}

ScamacErrorCode scamac_workspace_free(ScamacWorkspace * ws) {
  ScamacErrorCode err;
  if (ws) {
    if (ws->fct_work_free && ws->ws) {
      err = ws->fct_work_free(ws->ws);
      if (err) {
        return (err | SCAMAC_EINTERNAL);
      }
    }
    if (ws->row_real) {
      err = scamac_sparserow_real_free(ws->row_real);
      if (err) {
        return (err | SCAMAC_EINTERNAL);
      }
    }
    if (ws->row_cplx) {
      err = scamac_sparserow_cplx_free(ws->row_cplx);
      if (err) {
        return (err | SCAMAC_EINTERNAL);
      }
    }
    free(ws);
  }
  return SCAMAC_EOK;
}

ScamacErrorCode scamac_generator_check(const ScamacGenerator * gen, char ** desc) {
  ScamacErrorCode err;

  if (!gen) {
    return SCAMAC_ENULL | 1 << SCAMAC_ESHIFT;
  }

  /* for wrappers */
  if (gen->wrapped_par) {
    if (gen->fct_check_wrapped) {
      if (desc) {
        *desc = NULL;
      }
      err = gen->fct_check_wrapped(gen->wrapped_par, desc);
      if (err) {
        return err;
      }
    }
    if (!gen->needs_unwrapping) {
      // check unwrapped parameters. No error should occur here.
      err = gen->fct_check(gen->par, NULL);
      if (err) {
        return (err | SCAMAC_EINTERNAL);
      }
    }
    return SCAMAC_EOK;
  }

  if (gen->par) {
    if (gen->fct_check) {
      return gen->fct_check(gen->par, desc);
    } else {
      if (*desc) {
        *desc = malloc(100 * sizeof **desc);
        snprintf(*desc, 100, "Oops - internal error in: scamac_generator_check");
      }
      return SCAMAC_EFAIL | SCAMAC_EINTERNAL;
    }
  } else {
    if (desc) {
      *desc = NULL;
    }
    return SCAMAC_EOK;
  }

}

ScamacErrorCode scamac_generator_destroy(ScamacGenerator * gen) {
  ScamacErrorCode err;
  if (gen) {
    if (gen->fct_tables_destroy && gen->tables) {
      err = gen->fct_tables_destroy(gen->tables);
      if (err) {
        return err; // TODO
      }
    }
    if (gen->par) {
      free(gen->par);
    }
    if (gen->wrapped_par) {
      free(gen->wrapped_par);
    }
    free(gen);
  }
  return SCAMAC_EOK;
}


ScamacErrorCode scamac_generator_finalize(ScamacGenerator * gen) {
  ScamacErrorCode err;

  if (!gen) {
    return SCAMAC_ENULL | 1 << SCAMAC_ESHIFT;
  }

  if (gen -> needs_finalization) {
    err = scamac_generator_check(gen, NULL);
    if (SCAMAC_DISWARN(err)) {
      return err;
    }

    /* for wrappers: unwrap parameters  */
    if (gen -> fct_unwrap_par && gen -> wrapped_par) {
      err = gen->fct_unwrap_par(gen->wrapped_par, gen->par);
      if (SCAMAC_DISWARN(err)) {
        return err;
      }
      gen->needs_unwrapping=false;
      // check unwrapped parameters: No error should occur now.
      err = scamac_generator_check(gen, NULL);
      if (SCAMAC_DISWARN(err)) {
        return (err | SCAMAC_EINTERNAL);
      }
    }

    // have tables been allocated previously? Then, destroy them first.
    if (gen->fct_tables_destroy && gen->tables) {
      err = gen->fct_tables_destroy(gen->tables);
      if (err) {
        return err; // TODO
      }
      gen->tables = NULL;
    }
    if (gen->fct_tables_create) {
      err = gen->fct_tables_create(gen->par, &(gen->tables), &(gen->info) );
      if (err) {
        return err;
      }
    }

    gen->needs_finalization=false;
    return SCAMAC_EOK;
  } else {
    return SCAMAC_EOK;
  }
}


ScamacIdx scamac_generator_query_nrow    (const ScamacGenerator * gen) {
  if (gen) {
    return gen->info.nrow;
  } else {
    return 0;
  }
}

ScamacIdx scamac_generator_query_ncol    (const ScamacGenerator * gen) {
  if (gen) {
    return gen->info.ncol;
  } else {
    return 0;
  }
}

ScamacIdx scamac_generator_query_maxnzrow(const ScamacGenerator * gen) {
  if (gen) {
    return gen->info.maxnzrow;
  } else {
    return 0;
  }
}

ScamacIdx scamac_generator_query_maxnzcol(const ScamacGenerator * gen) {
  if (gen) {
    return gen->info.maxnzcol;
  } else {
    return 0;
  }
}

ScamacIdx scamac_generator_query_maxnz   (const ScamacGenerator * gen) {
  if (gen) {
    return gen->info.maxnz;
  } else {
    return 0;
  }
}

ScamacIdx scamac_generator_query_valtype (const ScamacGenerator * gen) {
  if (gen) {
    return gen->info.valtype;
  } else {
    return SCAMAC_NONE;
  }
}

ScamacIdx scamac_generator_query_symmetry(const ScamacGenerator * gen) {
  if (gen) {
    return gen->info.symmetry;
  } else {
    return SCAMAC_NONE;
  }
}

const char * scamac_generator_query_name(const ScamacGenerator * gen) {
  static const char msg[]="Unknown example";
  if (gen) {
    if (gen->name) {
      //my_string = malloc((strlen(gen->name) + 1) * sizeof *my_string);
      //strncpy(my_string,gen->name,strlen(gen->name)+1);
      return gen->name;
    } else {
      //my_string = malloc(16 * sizeof *my_string);
      //strncpy(my_string,"Unknown example",16);
      return msg;
    }
  } else {
    return msg;
  }
}

ScamacErrorCode scamac_generate_row(const ScamacGenerator * gen, ScamacWorkspace * ws, ScamacIdx irow, ScamacFlag flag, ScamacIdx * nzr, ScamacIdx * cind, double * val) {
  ScamacErrorCode err;

  if (flag & ~SCAMAC_TRANSPOSE & ~SCAMAC_CONJUGATE & ~SCAMAC_KEEPZEROS) {
    return SCAMAC_ERANGE;
  }
  bool fl_transpose = (flag & SCAMAC_TRANSPOSE) != 0;
  bool fl_keepzeros = (flag & SCAMAC_KEEPZEROS) != 0;

  if (!gen) {
    return SCAMAC_ENULL;
  }
  if (gen->needs_finalization) {
    return SCAMAC_EFAIL;
  }
  if (!gen->fct_gen_row) {
    return (SCAMAC_ENULL | SCAMAC_EINTERNAL);
  }

  if (!nzr) {
    return SCAMAC_ENULL;
  }
  if (val && (!cind)) {
    return SCAMAC_ENULL;
  }
  if (irow < 0) {
    return SCAMAC_ERANGE;
  }

  if (!ws) {
    return SCAMAC_ENULL;
  }
  if ( (!ws->row_real) && (!ws->row_cplx) ) {
    return (SCAMAC_ENULL     | SCAMAC_EINTERNAL);
  }
  if (   ws->row_real  &&   ws->row_cplx  ) {
    return (SCAMAC_EINVALID | SCAMAC_EINTERNAL);
  }

  ScamacIdx maxnzr;
  if (fl_transpose) {
    if ( (irow < 0) || (irow >= gen->info.ncol) ) {
      return SCAMAC_ERANGE;
    }
    maxnzr = gen->info.maxnzcol;
  } else {
    if ( (irow < 0) || (irow >= gen->info.nrow) ) {
      return SCAMAC_ERANGE;
    }
    maxnzr = gen->info.maxnzrow;
  }

  if (ws->row_real) {
    err = scamac_sparserow_real_zero(ws->row_real);
    if (err) {
      return (err | SCAMAC_EINTERNAL);
    }

    err = gen->fct_gen_row(gen->par, gen->tables, ws->ws, irow, flag, ws->row_real);
    if (err) {
      return (err | SCAMAC_EINTERNAL);
    }

    err = scamac_sparserow_real_normalize(ws->row_real, fl_keepzeros);
    if (err) {
      return (err | SCAMAC_EINTERNAL);
    }

    err = scamac_sparserow_real_to_idxval(ws->row_real, maxnzr, nzr, cind, val);
    if (err) {
      return (err | SCAMAC_EINTERNAL);
    }

  } else if (ws->row_cplx) {
    err = scamac_sparserow_cplx_zero(ws->row_cplx);
    if (err) {
      return (err | SCAMAC_EINTERNAL);
    }

    err = gen->fct_gen_row(gen->par, gen->tables, ws->ws, irow, flag, ws->row_cplx);
    if (err) {
      return (err | SCAMAC_EINTERNAL);
    }

    err = scamac_sparserow_cplx_normalize(ws->row_cplx, fl_keepzeros);
    if (err) {
      return (err | SCAMAC_EINTERNAL);
    }

    err = scamac_sparserow_cplx_to_idxval(ws->row_cplx, maxnzr, nzr, cind, (double complex *) val);
    if (err) {
      return (err | SCAMAC_EINTERNAL);
    }
  }

  return SCAMAC_EOK;
}

ScamacErrorCode scamac_alloc_cind_val(const ScamacGenerator * gen, ScamacFlag flag, ScamacIdx ** cind, double ** val) {
  if (!gen) {
    return SCAMAC_ENULL;
  }
  if (gen->needs_finalization) {
    return SCAMAC_EFAIL;
  }

  if (flag & ~SCAMAC_TRANSPOSE & ~SCAMAC_CONJUGATE & ~SCAMAC_KEEPZEROS) {
    return SCAMAC_ERANGE;
  }
  bool fl_transpose = (flag & SCAMAC_TRANSPOSE) != 0;

  ScamacIdx maxnzr;
  if (fl_transpose) {
    maxnzr = gen->info.maxnzcol;
  } else {
    maxnzr = gen->info.maxnzrow;
  }

  if (cind) {
    if (maxnzr>0) {
      *cind = malloc(maxnzr * sizeof **cind);
      if (! *cind) {
        return SCAMAC_EMALLOCFAIL;
      }
    } else {
      *cind = NULL;
    }
  }
  if (val) {
    if (maxnzr>0) {
      if (gen->info.valtype == SCAMAC_VAL_REAL) {
        *val = malloc(maxnzr * sizeof **val);
      } else if (gen->info.valtype == SCAMAC_VAL_COMPLEX) {
        *val = malloc(2 * maxnzr * sizeof **val);
      } else {
        return SCAMAC_EFAIL | SCAMAC_EINTERNAL;
      }
    } else {
      *val = NULL;
    }
    if (! *val) {
      return SCAMAC_EMALLOCFAIL;
    }
  }
  return SCAMAC_EOK;
}

ScamacErrorCode scamac_generate_row_real    (const ScamacGenerator * gen, ScamacWorkspace * ws, ScamacIdx irow, ScamacFlag flag, ScamacIdx * nzr, ScamacIdx * cind, double * val) {
  if (!gen) {
    return SCAMAC_ENULL;
  }
  if (gen->needs_finalization) {
    return SCAMAC_EFAIL;
  }
  if (gen->info.valtype != SCAMAC_VAL_REAL) {
    return SCAMAC_EFAIL;
  }
  ScamacErrorCode err;
  err = scamac_generate_row(gen,ws,irow,flag,nzr,cind,val);
  return err;
}

ScamacErrorCode scamac_generate_row_cplx    (const ScamacGenerator * gen, ScamacWorkspace * ws, ScamacIdx irow, ScamacFlag flag, ScamacIdx * nzr, ScamacIdx * cind, double complex * val) {
  if (!gen) {
    return SCAMAC_ENULL;
  }
  if (gen->needs_finalization) {
    return SCAMAC_EFAIL;
  }
  if (gen->info.valtype != SCAMAC_VAL_COMPLEX) {
    return SCAMAC_EFAIL;
  }
  ScamacErrorCode err;
  err = scamac_generate_row(gen,ws,irow,flag,nzr,cind,(double *) val);
  return err;
}


ScamacErrorCode scamac_generate_row_int     (const ScamacGenerator * gen, ScamacWorkspace * ws, ScamacIdx irow, ScamacFlag flag, int * nzr, int * cind, double * val) {
  /* we duplicate the code of scamac_generate_row, with an additional check on maxnzr, and a call to scamac_sparserow_..._intval
   * Probably, both fcts. ( scamac_generate_row and scamac_generate_row_int) should be merged.
   */

  ScamacErrorCode err;

  if (flag & ~SCAMAC_TRANSPOSE & ~SCAMAC_CONJUGATE & ~SCAMAC_KEEPZEROS) {
    return SCAMAC_ERANGE;
  }
  bool fl_transpose = (flag & SCAMAC_TRANSPOSE) != 0;
  bool fl_keepzeros = (flag & SCAMAC_KEEPZEROS) != 0;

  if (!gen) {
    return SCAMAC_ENULL;
  }
  if (gen->needs_finalization) {
    return SCAMAC_EFAIL;
  }
  if (!gen->fct_gen_row) {
    return (SCAMAC_ENULL | SCAMAC_EINTERNAL);
  }

  if (!nzr) {
    return SCAMAC_ENULL;
  }
  if (val && (!cind)) {
    return SCAMAC_ENULL;
  }
  if (irow < 0) {
    return SCAMAC_ERANGE;
  }

  if (!ws) {
    return SCAMAC_ENULL;
  }
  if ( (!ws->row_real) && (!ws->row_cplx) ) {
    return (SCAMAC_ENULL     | SCAMAC_EINTERNAL);
  }
  if (   ws->row_real  &&   ws->row_cplx  ) {
    return (SCAMAC_EINVALID | SCAMAC_EINTERNAL);
  }

  if ( (gen->info.ncol > INT_MAX) || (gen->info.nrow > INT_MAX) ) {
    return SCAMAC_EOVERFLOW;
  }

  ScamacIdx maxnzr;
  if (fl_transpose) {
    if ( (irow < 0) || (irow >= gen->info.ncol) ) {
      return SCAMAC_ERANGE;
    }
    maxnzr = gen->info.maxnzcol;
  } else {
    if ( (irow < 0) || (irow >= gen->info.nrow) ) {
      return SCAMAC_ERANGE;
    }
    maxnzr = gen->info.maxnzrow;
  }

  if (ws->row_real) {
    err = scamac_sparserow_real_zero(ws->row_real);
    if (err) {
      return (err | SCAMAC_EINTERNAL);
    }

    err = gen->fct_gen_row(gen->par, gen->tables, ws->ws, irow, flag, ws->row_real);
    if (err) {
      return (err | SCAMAC_EINTERNAL);
    }

    err = scamac_sparserow_real_normalize(ws->row_real, fl_keepzeros);
    if (err) {
      return (err | SCAMAC_EINTERNAL);
    }

    err = scamac_sparserow_real_to_idxval_int(ws->row_real, maxnzr, nzr, cind, val);
    if (err) {
      return (err | SCAMAC_EINTERNAL);
    }

  } else if (ws->row_cplx) {
    err = scamac_sparserow_cplx_zero(ws->row_cplx);
    if (err) {
      return (err | SCAMAC_EINTERNAL);
    }

    err = gen->fct_gen_row(gen->par, gen->tables, ws->ws, irow, flag, ws->row_cplx);
    if (err) {
      return (err | SCAMAC_EINTERNAL);
    }

    err = scamac_sparserow_cplx_normalize(ws->row_cplx, fl_keepzeros);
    if (err) {
      return (err | SCAMAC_EINTERNAL);
    }

    err = scamac_sparserow_cplx_to_idxval_int(ws->row_cplx, maxnzr, nzr, cind, (double complex *) val);
    if (err) {
      return (err | SCAMAC_EINTERNAL);
    }
  }

  return SCAMAC_EOK;
}

ScamacErrorCode scamac_generate_row_int_real(const ScamacGenerator * gen, ScamacWorkspace * ws, ScamacIdx irow, ScamacFlag flag, int * nzr, int * cind, double * val) {
  if (!gen) {
    return SCAMAC_ENULL;
  }
  if (gen->needs_finalization) {
    return SCAMAC_EFAIL;
  }
  if (gen->info.valtype != SCAMAC_VAL_REAL) {
    return SCAMAC_EFAIL;
  }
  ScamacErrorCode err;
  err = scamac_generate_row_int(gen,ws,irow,flag,nzr,cind,val);
  return err;
}

ScamacErrorCode scamac_generate_row_int_cplx(const ScamacGenerator * gen, ScamacWorkspace * ws, ScamacIdx irow, ScamacFlag flag, int * nzr, int * cind, double complex * val) {
  if (!gen) {
    return SCAMAC_ENULL;
  }
  if (gen->needs_finalization) {
    return SCAMAC_EFAIL;
  }
  if (gen->info.valtype != SCAMAC_VAL_COMPLEX) {
    return SCAMAC_EFAIL;
  }
  ScamacErrorCode err;
  err = scamac_generate_row_int(gen,ws,irow,flag,nzr,cind,(double *) val);
  return err;
}

