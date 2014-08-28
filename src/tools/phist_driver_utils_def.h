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
