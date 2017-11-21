/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include <string.h>

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

#define PHIST_CLASSFILE_DEF "phist_driver_utils_def.h"
#include "phist_gen_all.h"
