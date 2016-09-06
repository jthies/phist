#ifdef __cplusplus
#include <cstring>
#else
#include <string.h>
#endif

#include "phist_config.h"
#include "phist_kernels.h"

// some helper functions

int phist_sizeof_lidx()
{
  return sizeof(phist_lidx);
}

int phist_sizeof_gidx()
{
  return sizeof(phist_gidx);
}

const char* phist_file_extension(const char* filename)
{
  const char *c=filename+strlen(filename)-1;
  while (c!=filename)
  {
    if (*c=='.') break;
    c--;
  }
  c++;
  if (c==filename) c=NULL;
  return c;
}

int phist_filename_isMM(const char* filename)
{
  const char *c = phist_file_extension(filename);
  return c? (!strcmp(c,"mm"))||(!strcmp(c,"mtx")): 0;
}

int phist_filename_isCRS(const char* filename)
{
  const char *c = phist_file_extension(filename);
  return c? (!strcmp(c,"crs"))||(!strcmp(c,"bin")): 0;
}

int phist_filename_isHB(const char* filename)
{
  const char *c = phist_file_extension(filename);
  return c? (!strcmp(c,"rua"))||(!strcmp(c,"cua"))||(!strcmp(c,"rsa"))||(!strcmp(c,"csa")): 0;
}

#ifdef PHIST_HAVE_SP
#include "phist_gen_s.h"
#include "phist_driver_utils_def.h"
#include "phist_gen_c.h"
#include "phist_driver_utils_def.h"
#endif
#include "phist_gen_d.h"
#include "phist_driver_utils_def.h"
#include "phist_gen_z.h"
#include "phist_driver_utils_def.h"

