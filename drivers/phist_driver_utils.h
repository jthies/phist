#ifndef PHIST_DRIVER_UTILS_H
#define PHIST_DRIVER_UTILS_H
#ifdef __cplusplus
#include <cstring>
#else
#include <string.h>
#endif
typedef _ST_ ST;
typedef _MT_ MT;
typedef TYPE(mvec_ptr) mvec_ptr_t;
typedef TYPE(const_mvec_ptr) const_mvec_ptr_t;

typedef TYPE(sdMat_ptr) sdMat_ptr_t;
typedef TYPE(const_sdMat_ptr) const_sdMat_ptr_t;

typedef TYPE(crsMat_ptr) crsMat_ptr_t;
typedef TYPE(const_crsMat_ptr) const_crsMat_ptr_t;

typedef TYPE(op_ptr) op_ptr_t;
typedef TYPE(const_op_ptr) const_op_ptr_t;

// auto-detects the file type by looking at the file extension
void SUBR(crsMat_read)(TYPE(crsMat_ptr)* A, char* filename, int* ierr)
  {
  char *c=&filename[0];
  int isMM=0;
  int isCRS=0;
  int isHB=0;
  
  while (*c)
    {
    if (*c=='.') break;
    c++;
    }
  c++;
  isMM=!strcmp(c,"mm");
  isCRS=!strcmp(c,"crs");
  isCRS=isCRS||!strcmp(c,"bin");
  isHB=!strcmp(c,"rua")||!strcmp(c,"cua");
  if (isMM)
    {
    PHIST_CHK_IERR(SUBR(crsMat_read_mm)(A,filename,ierr),*ierr);
    }
  else if (isHB)
    {
    PHIST_CHK_IERR(SUBR(crsMat_read_hb)(A,filename,ierr),*ierr);
    }
  else if (isCRS)
    {
    PHIST_CHK_IERR(SUBR(crsMat_read_bin)(A,filename,ierr),*ierr);
    }
  else
    {
    *ierr=-1;
    }
  return;
  }
#endif
