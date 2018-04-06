/** \file
 *  \author Andreas Alvermann (<alvermann@physik.uni-greifswald.de>)
 *  \date   October 2017
 *  \brief  ...
 *  \ingroup internal
 */

#ifndef SCAMAC_STRING_H
#define SCAMAC_STRING_H

typedef struct {
  int nalloc;
  int nstr;
  char * str;
} scamac_string_st;

/* extendable string routines */

void scamac_string_empty(scamac_string_st * estr);
void scamac_string_append(scamac_string_st * estr, const char * str);
char * scamac_string_get(const scamac_string_st * estr);

#endif /* SCAMAC_STRING_H */
