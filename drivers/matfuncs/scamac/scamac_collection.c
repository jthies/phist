#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "scamac_include.h"
#include "scamac_collection.h"
#include "scamac_internal.h"

// includes headers of all matrix examples. Autogenerated.
#include "scamac_collection_inc.h"

ScamacErrorCode scamac_generator_obtain(const char * matname, ScamacGenerator ** gen) {

  //include code for all matrix examples. Autogenerated.
#include "scamac_collection_example_inc.c"

  return SCAMAC_EFAIL;
}

ScamacErrorCode scamac_generator_set_int(ScamacGenerator * gen, const char * parname, int val) {
  if (!gen) {
    return SCAMAC_ENULL | 1 << SCAMAC_ESHIFT;
  }
  if (!parname) {
    return SCAMAC_ENULL | 2 << SCAMAC_ESHIFT;
  }
  
  gen->needs_finalization = true; // mark as "tainted"

#include "scamac_collection_set_int_inc.c"

  // if we made it until here, something has gone wrong
  return SCAMAC_EFAIL | SCAMAC_EINTERNAL;
}

ScamacErrorCode scamac_generator_set_double(ScamacGenerator * gen, const char * parname, double val) {
  if (!gen) {
    return SCAMAC_ENULL | 1 << SCAMAC_ESHIFT;
  }
  if (!parname) {
    return SCAMAC_ENULL | 2 << SCAMAC_ESHIFT;
  }
 
  gen->needs_finalization = true; // mark as "tainted"

#include "scamac_collection_set_double_inc.c"

  // if we made it until here, something has gone wrong
  return SCAMAC_EFAIL | SCAMAC_EINTERNAL;
}

ScamacErrorCode scamac_generator_set_bool(ScamacGenerator * gen, const char * parname, bool val) {
  if (!gen) {
    return SCAMAC_ENULL | 1 << SCAMAC_ESHIFT;
  }
  if (!parname) {
    return SCAMAC_ENULL | 2 << SCAMAC_ESHIFT;
  } 
 
  gen->needs_finalization = true; // mark as "tainted"

#include "scamac_collection_set_bool_inc.c"

  // if we made it until here, something has gone wrong
  return SCAMAC_EFAIL | SCAMAC_EINTERNAL;
}


ScamacErrorCode scamac_generator_set_seed(ScamacGenerator * gen, const char * parname, uint64_t seed) {
  if (!gen) {
    return SCAMAC_ENULL | 1 << SCAMAC_ESHIFT;
  }
  if (!parname) {
    return SCAMAC_ENULL | 2 << SCAMAC_ESHIFT;
  }
  
  gen->needs_finalization = true; // mark as "tainted"

#include "scamac_collection_set_rngseed_inc.c"

  // if we made it until here, something has gone wrong
  return SCAMAC_EFAIL | SCAMAC_EINTERNAL;
}

ScamacErrorCode scamac_generator_set_seed_str(ScamacGenerator * gen, const char * parname, const char * seedstr) {
  if (!gen) {
    return SCAMAC_ENULL | 1 << SCAMAC_ESHIFT;
  }
  if (!parname) {
    return SCAMAC_ENULL | 2 << SCAMAC_ESHIFT;
  }
  
  gen->needs_finalization = true; // mark as "tainted"
  
#include "scamac_collection_set_rngseed_str_inc.c"

  // if we made it until here, something has gone wrong
  return SCAMAC_EFAIL | SCAMAC_EINTERNAL;
}

ScamacErrorCode scamac_generator_set_option(ScamacGenerator * gen, const char * parname, const char * option) {
  if (!gen) {
    return SCAMAC_ENULL | 1 << SCAMAC_ESHIFT;
  }
  if (!parname) {
    return SCAMAC_ENULL | 2 << SCAMAC_ESHIFT;
  }
  if (!option) {
    return SCAMAC_ENULL | 3 << SCAMAC_ESHIFT;
  }
  
  gen->needs_finalization = true; // mark as "tainted"

#include "scamac_collection_set_option_inc.c"

  // if we made it until here, something has gone wrong
  return SCAMAC_EFAIL | SCAMAC_EINTERNAL;
}

ScamacErrorCode scamac_list_examples(char ** desc) {
  if (desc) {
    char * my_string;
#include "scamac_collection_list_examples_inc.c"
    *desc = my_string;
    return SCAMAC_EOK;
  } else {
    return SCAMAC_ENULL;
  }
}

ScamacErrorCode scamac_example_desc(const char * matname, int * valtype, int * symmetry, char ** desc) {
  if (!matname) {
    return SCAMAC_ENULL | 1 << SCAMAC_ESHIFT;
  }

#include "scamac_collection_example_desc_inc.c"

  return SCAMAC_EINVALID | 1 << SCAMAC_ESHIFT;
}

ScamacErrorCode scamac_example_parameters(const char * matname, char ** desc) {
  if (!matname) {
    return SCAMAC_ENULL | 1 << SCAMAC_ESHIFT;
  }
   if (!desc) {
    return SCAMAC_ENULL | 2 << SCAMAC_ESHIFT;
  }
 
  char * my_string;
#include "scamac_collection_list_parameters_inc.c"

  // unknown example
  if (*desc) {
    *desc=NULL;
  }
  return SCAMAC_EINVALID | 1 << SCAMAC_ESHIFT;
}

ScamacPar scamac_identify_parameter(const char * matname, const char * parname) {
  if (matname && parname) {
#include "scamac_collection_identify_parameter_inc.c"
    return SCAMAC_PAR_NONE; // not found
  } else {
    return SCAMAC_PAR_NONE; // default fall back
  }
}

static char * strip_string(char * s) {
  if (s) {
    while (*s == ' ') {
      s++;
    }
    int n = strlen(s);
    while (n>0 && s[n-1] == ' ') {
      s[n-1]=0;
      n--;
    }
    return s;
  } else {
    return NULL;
  }
}

static void delete_last_char(const char * ch, char *s) {
  int l = strlen(s);
  if (l>0) {
    if (s[l-1]==ch[0]) {
      s[l-1]=0;
    }
  }
}

ScamacErrorCode scamac_generator_parameter_desc(const ScamacGenerator * gen, const char * format, char ** desc) {
  if (!gen) {
    return SCAMAC_ENULL | 1 << SCAMAC_ESHIFT;
  }
  if (!format) {
    return SCAMAC_ENULL | 2 << SCAMAC_ESHIFT;
  }
  if (!desc) {
    return SCAMAC_ENULL | 3 << SCAMAC_ESHIFT;
  }
  
  bool print_name = false;
  bool print_as_argstr = false;
  bool print_double_hex = false;
      
  if (!strcmp(format,"argstr")) {
    print_name=true;
    print_as_argstr = true;
    print_double_hex = false;
    } else if (!strcmp(format,"desc")) {
    print_name=true;
    print_as_argstr = false;
    print_double_hex = false;
  } else {
    return SCAMAC_EINVALID | 2 << SCAMAC_ESHIFT;
  }
    
  
  // separation character
  char sepchar[2];
  if (print_as_argstr) {
    snprintf(sepchar, 2, "%s", ",");
  } else {
    snprintf(sepchar, 2, "%s", "\n");
  }
  
  // format strings
  char format_name[30],format_int[30],format_double[30], format_bool[30], format_rngseed[30], format_option[30];

  snprintf(format_name,   30, "%%s%s", sepchar);
  snprintf(format_int,    30, "%%s=%%d%s",sepchar);
  if (print_double_hex) {
    snprintf(format_double,   30, "%%s=%%a%s",sepchar);
  } else {
    if (print_as_argstr) {
      snprintf(format_double, 30, "%%s=%%1.16e%s",sepchar);
    } else {
      snprintf(format_double, 30, "%%s=%%1.16g%s",sepchar);
    }
  }
  snprintf(format_bool,   30, "%%s=%%s%s",sepchar);
  snprintf(format_rngseed,30, "%%s=%%"PRIu64"%s",sepchar);
  snprintf(format_option, 30, "%%s=%%s%s",sepchar);
  
  *desc = NULL;
  
  char * my_string;
#include "scamac_collection_example_data_inc.c"

  return SCAMAC_ECORRUPTED | 1 << SCAMAC_ESHIFT;
}

static int split_string(const char *s, char **split, char ***sfst) {
  if (s) {

    // strdup
    char * s2;
    s2 = malloc((strlen(s)+1) * sizeof *s2);
    memcpy(s2, s, strlen(s) * sizeof *s2);
    s2[strlen(s)]=0;

    // first: count

    int n=0;
    char * s1 = s2;

    while (s1) {
      n++;
      s1 = strpbrk(s1, ",");
      if (s1) {
        s1++;
      }
    }

    // second: split

    char ** my_list = malloc( 2* n * sizeof *my_list);

    // split at ","
    s1 = s2;
    n=0;
    while (s1) {
      my_list[2*n]=s1;
      my_list[2*n+1]=NULL;
      n++;
      s1 = strpbrk(s1, ",");
      if (s1) {
        *s1=0; // NULL terminate string
        s1++;
      }
    }

    // split each at "="
    int i;
    for (i=0; i<n; i++) {
      s1 = strpbrk(my_list[2*i], "=");
      if (s1) {
        *s1=0;
        my_list[2*i+1]=s1+1;
      }
    }

    // strip
    for (i=0; i<2*n; i++) {
      if (my_list[i]) {
        my_list[i]=strip_string(my_list[i]);
      }
    }

    *split = s2;
    *sfst = my_list;
    return n;

  } else {
    *split=NULL;
    *sfst=NULL;
    return 0;
  }
}

static ScamacErrorCode atobool(const char * bs, bool * val) {
  if (!bs) {
    return SCAMAC_ENULL | (1 << SCAMAC_ESHIFT);
  }
  if (!val) {
    return SCAMAC_ENULL | (2 << SCAMAC_ESHIFT);
  }
  if (!strcmp(bs,"true") || !strcmp(bs,"True") || !strcmp(bs,"TRUE")) {
    *val = true;
    return SCAMAC_EOK;
  } else if (!strcmp(bs,"false") || !strcmp(bs,"False") || !strcmp(bs,"FALSE")) {
    *val = false;
    return SCAMAC_EOK;
  } else {
    return SCAMAC_EINVALID | (2 << SCAMAC_ESHIFT);
  }
}

ScamacErrorCode scamac_parse_argstr(const char * argstr, ScamacGenerator **gen, char ** errdesc) {

  ScamacErrorCode err=SCAMAC_EOK;

  if (!argstr) {
    if (errdesc) {
      char *str = malloc(SCAMAC_NAME_LENGTH * sizeof *str);
      snprintf(str,SCAMAC_NAME_LENGTH,"missing argstr");
      *errdesc = str;
    }
    return SCAMAC_ENULL;
  }

  if (!gen) {
    if (errdesc) {
      char *str = malloc(SCAMAC_NAME_LENGTH * sizeof *str);
      snprintf(str,SCAMAC_NAME_LENGTH,"gen. pointer = NULL");
      *errdesc = str;
    }
    return SCAMAC_ENULL;
  }

  ScamacGenerator * my_gen;

  char * my_s=NULL;
  char ** my_slst=NULL;

  // TODO : free memory allocated here
  int n = split_string(argstr, &my_s, &my_slst);


  if (n<1) {
    if (errdesc) {
      char *str = malloc(SCAMAC_NAME_LENGTH * sizeof *str);
      snprintf(str,SCAMAC_NAME_LENGTH,"few parameters");
      *errdesc = str;
    }
    return SCAMAC_EFAIL;
  }

  // matrix name
  if (my_slst[1]) {
    if (errdesc) {
      char *str = malloc(SCAMAC_NAME_LENGTH * sizeof *str);
      snprintf(str,SCAMAC_NAME_LENGTH,"matrix name accepts no values");
      *errdesc = str;
    }
    return SCAMAC_EFAIL;
  }

  if (!my_slst[0]) {
    if (errdesc) {
      char *str = malloc(SCAMAC_NAME_LENGTH * sizeof *str);
      snprintf(str,SCAMAC_NAME_LENGTH,"matrix name required");
      *errdesc = str;
    }
    return SCAMAC_EFAIL;
  }

  err =  scamac_generator_obtain(my_slst[0], &my_gen);

  if (err) {
    if (errdesc) {
      char *str = malloc(SCAMAC_NAME_LENGTH * sizeof *str);
      snprintf(str,SCAMAC_NAME_LENGTH,"matrix >%s< does not exist",my_slst[0]);
      *errdesc = str;
    }
    return SCAMAC_EFAIL;
  }

  char *my_par = NULL;
  err = scamac_example_parameters(my_slst[0], &my_par);

  if (err) {
    if (errdesc) {
      char *str = malloc(SCAMAC_NAME_LENGTH * sizeof *str);
      snprintf(str,SCAMAC_NAME_LENGTH,"matrix >%s< has no parameters",my_slst[0]);
      *errdesc = str;
    }
    return SCAMAC_EFAIL;
  }


  // set parameters
  int i;
  int n_read;
  char dummy[2];
  bool err_read=false;
  for (i=1; i<n; i++) {
    ScamacPar par_type = scamac_identify_parameter(my_slst[0],my_slst[2*i]);
    if (par_type == SCAMAC_PAR_NONE) {
      if (errdesc) {
        char *str = malloc(SCAMAC_NAME_LENGTH * sizeof *str);
        snprintf(str,SCAMAC_NAME_LENGTH,"Parameter %s does not exist",my_slst[2*i]);
        *errdesc = str;
      }
      err=SCAMAC_EFAIL;
    } else if (par_type == SCAMAC_PAR_INT) {
      if (!my_slst[2*i+1]) {
        err_read=true;
      } else {
        int par;
        n_read = sscanf(my_slst[2*i+1],"%9d %1s",&par,dummy); // read at most 999999999 < 2^31
        if (n_read != 1) {
          err_read = true;
        } else {
          scamac_generator_set_int(my_gen, my_slst[2*i], par);  
        }
      }
      if (err_read) {
        if (errdesc) {
          char *str = malloc(SCAMAC_NAME_LENGTH * sizeof *str);
          snprintf(str,SCAMAC_NAME_LENGTH,"Parameter %s requires INT value",my_slst[2*i]);
          *errdesc = str;
        }
        err=SCAMAC_EFAIL;
      }
    } else if (par_type == SCAMAC_PAR_DOUBLE) {
      if (!my_slst[2*i+1]) {
        err_read=true;
      } else {
        double par;
        n_read = sscanf(my_slst[2*i+1],"%lf %1s",&par,dummy);
        if (n_read != 1) {
          err_read = true;
        } else {
          if (isnan(par) || isinf(par)) {
            err_read = true;
          } else {
            scamac_generator_set_double(my_gen, my_slst[2*i], par);  
          }
        }
      }
      if (err_read) {
        if (errdesc) {
          char *str = malloc(SCAMAC_NAME_LENGTH * sizeof *str);
          snprintf(str,SCAMAC_NAME_LENGTH,"Parameter %s requires DOUBLE value",my_slst[2*i]);
          *errdesc = str;
        }
        err=SCAMAC_EFAIL;
      } 
    } else if (par_type == SCAMAC_PAR_BOOL) {
      if (!my_slst[2*i+1]) {
        err_read=true;
      } else {
        bool par;
        err = atobool(my_slst[2*i+1],&par);
        if (err) {
          err_read = true;
        } else {
          scamac_generator_set_bool(my_gen, my_slst[2*i], par);  
        }
      }
      if (err_read) {
        if (errdesc) {
          char *str = malloc(SCAMAC_NAME_LENGTH * sizeof *str);
          snprintf(str,SCAMAC_NAME_LENGTH,"Parameter %s requires BOOL value",my_slst[2*i]);
          *errdesc = str;
        }
        err=SCAMAC_EFAIL;
      }
    } else if (par_type == SCAMAC_PAR_RNGSEED) {
      if (!my_slst[2*i+1]) {
        err_read=true;
      } else {
        if (strlen(my_slst[2*i+1])<1) {
          err_read=true;
        } else {
        scamac_generator_set_seed_str(my_gen, my_slst[2*i], my_slst[2*i+1] ); // every string is accepted
        }
      }
      if (err_read) {
        if (errdesc) {
          char *str = malloc(SCAMAC_NAME_LENGTH * sizeof *str);
          snprintf(str,SCAMAC_NAME_LENGTH,"Parameter %s requires RNGSEED (dec,hex,string) value",my_slst[2*i]);
          *errdesc = str;
        }
        err=SCAMAC_EFAIL;
      } 
    } else if (par_type == SCAMAC_PAR_OPTION) {
      if (!my_slst[2*i+1]) {
         err_read=true;
      } else {
        if (strlen(my_slst[2*i+1])<1) {
          err_read=true;
        } else {
          err = scamac_generator_set_option(my_gen, my_slst[2*i], my_slst[2*i+1]);
          if (err) {
            if (errdesc) {
              char *str = malloc(SCAMAC_NAME_LENGTH * sizeof *str);
              snprintf(str,SCAMAC_NAME_LENGTH,"Unknown option for parameter %s : %s",my_slst[2*i],my_slst[2*i+1]);
              *errdesc = str;
            }
            err=SCAMAC_EFAIL;
          }
        }
      }
      if (err_read) {
        if (errdesc) {
          char *str = malloc(SCAMAC_NAME_LENGTH * sizeof *str);
          snprintf(str,SCAMAC_NAME_LENGTH,"Parameter %s requires OPTION (string) value",my_slst[2*i]);
          *errdesc = str;
        }
        err=SCAMAC_EFAIL;
      }
    } else {
      if (errdesc) {
        char *str = malloc(SCAMAC_NAME_LENGTH * sizeof *str);
        snprintf(str,SCAMAC_NAME_LENGTH,"Unknown parameter type");
        *errdesc = str;
      }
      err=SCAMAC_EFAIL;
    }
    if (err) { break; }
  }

  free(my_s);
  free(my_slst);
  free(my_par);

  if (!err) {// SCAMAC_EOK
    *gen =  my_gen;
    if (errdesc) {
      *errdesc = NULL;
    }
  }
  
  return err;

}

const char * scamac_error_desc(ScamacErrorCode err) {
  const int MaxErrorN = 14;
  /* This array must agree with the def. of ScamacErrorCode in scamac_include.h */
  static const char *ErrorStrings[] = {
    /* OK */           "== OK: Hurrah! ==",
    /* FAIL */         "== error FAIL (Failed for no reason) ==",
    /* NULL */         "== error NULL (Unexpected NULL pointer) ==",
    /* INVALID */      "== error INVALID (Invalid parameter) ==",
    /* RANGE */        "== error RANGE (Parameter out of range) ==",
    /* CORRUPTED */    "== error CORRUPTED (Object is corrupted) ==",
    /* SCOPE */        "== error SCOPE (Requested omputation not covered) ==",
    /* INPUT */        "== error INPUT (Invalid set of input parameters) ==",
    /* OVERFLOW */     "== error OVERFLOW (Matrix size exceeds INT_MAX) ==",
    /* MALLOCFAIL */   "== error MALLOCFAIL (Memory allocation failed) ==",
    /* HUGEINT */      "== error HUGEINT (An integer parameters is too large) ==",
    /* HUGECOMP */     "== error HUGECOMP (A computed integer is too large) ==",
    /* WARNING */      "== error WARNING (A warning -- might be ignored) ==",
    /* NOTCONVERGED */ "== error NOTCONVERGED (Failed to converged) ==",
    /* SHORTROW */     "== error SHORTROW (Short row) =="
  };
  static const char *InternalErrorStrings[] = {
    /* OK */           "== OK: Hurrah! ==",
    /* FAIL */         "== internal error FAIL (Failed for no reason) ==",
    /* NULL */         "== internal error NULL (Unexpected NULL pointer) ==",
    /* INVALID */      "== internal error INVALID (Invalid parameter) ==",
    /* RANGE */        "== internal error RANGE (Parameter out of range) ==",
    /* CORRUPTED */    "== internal error CORRUPTED (Object is corrupted) ==",
    /* SCOPE */        "== internal error SCOPE (Requested omputation not covered) ==",
    /* INPUT */        "== internal error INPUT (Invalid set of input parameters) ==",
    /* OVERFLOW */     "== internal error OVERFLOW (Matrix size exceeds INT_MAX) ==",
    /* MALLOCFAIL */   "== internal error MALLOCFAIL (Memory allocation failed) ==",
    /* HUGEINT */      "== internal error HUGEINT (An integer parameters is too large) ==",
    /* HUGECOMP */     "== internal error HUGECOMP (A computed integer is too large) ==",
    /* WARNING */      "== internal error WARNING (A warning -- might be ignored) ==",
    /* NOTCONVERGED */ "== internal error NOTCONVERGED (Failed to converged) ==",
    /* SHORTROW */     "== internal error SHORTROW (Short row) =="
  };

  int errno = err & SCAMAC_EMASK;
  int isinternal = err & SCAMAC_EINTERNAL;
  if (errno>MaxErrorN) {
    errno = 1; // (int) cast to avoid "tautological compare warning"
  }
  if (isinternal) {
    return InternalErrorStrings[errno];
  } else {
    return ErrorStrings[errno];
  }
}


int scamac_error_par (ScamacErrorCode err) {
  int parno = (err >> SCAMAC_ESHIFT) & SCAMAC_EMASK;
  return parno;
}

const char * scamac_error_dpar(ScamacErrorCode err) {
  static const char *ErrorParStrings[] = {
    "",
    " at parameter (1) ==",
    " at parameter (2) ==",
    " at parameter (3) ==",
    " at parameter (4) ==",
    " at parameter (5) ==",
    " at parameter (6) ==",
    " at parameter (7) ==",
    " at parameter (8) ==",
    " at parameter (9) ==",
    " at parameter (10) ==",
    " at parameter (11) ==",
    " at parameter (12) ==",
    " at parameter (13) ==",
    " at parameter (14) ==",
    " at parameter (15) ==",
    " at parameter (16) ==",
    " at parameter (17) ==",
    " at parameter (18) ==",
    " at parameter (19) ==",
    " at parameter (20) ==",
    " at parameter (21) ==",
    " at parameter (22) ==",
    " at parameter (23) ==",
    " at parameter (24) ==",
    " at parameter (25) ==",
    " at parameter (26) ==",
    " at parameter (27) ==",
    " at parameter (28) ==",
    " at parameter (29) ==",
    " at parameter (30) ==",
    " at parameter (31) ==",
    " at parameter (32) ==",
    " at parameter (33) ==",
    " at parameter (34) ==",
    " at parameter (35) ==",
    " at parameter (36) ==",
    " at parameter (37) ==",
    " at parameter (38) ==",
    " at parameter (39) ==",
    " at parameter (40) ==",
    " at parameter (41) ==",
    " at parameter (42) ==",
    " at parameter (43) ==",
    " at parameter (44) ==",
    " at parameter (45) ==",
    " at parameter (46) ==",
    " at parameter (47) ==",
    " at parameter (48) ==",
    " at parameter (49) ==",
    " at parameter (50) ==",
    " at parameter (51) ==",
    " at parameter (52) ==",
    " at parameter (53) ==",
    " at parameter (54) ==",
    " at parameter (55) ==",
    " at parameter (56) ==",
    " at parameter (57) ==",
    " at parameter (58) ==",
    " at parameter (59) ==",
    " at parameter (60) ==",
    " at parameter (61) ==",
    " at parameter (62) ==",
    " at parameter (63) =="
  };
  int parno = scamac_error_par(err);
  int errno = err & SCAMAC_EMASK;
  if (errno == SCAMAC_ENULL || errno == SCAMAC_EINVALID || errno == SCAMAC_ERANGE ||
      errno == SCAMAC_ECORRUPTED || errno == SCAMAC_EWARNING) {
    return ErrorParStrings[parno];
  } else {
    return ErrorParStrings[0];
  }
}