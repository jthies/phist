void SUBR(read_mat)(const char* filebase,const_comm_ptr_t comm, int nglob,TYPE(sparseMat_ptr) *ptr, int* iflag)
{
  int _iflag;
  *ptr = NULL;
  char tpc = ::phist::ScalarTraits< _ST_ >::type_char();
  char mmfile[256],hbfile[256],binfile[256],crsfile[256];
  sprintf(mmfile,"%c%s%d.mm",tpc,filebase,nglob);
  if (::phist::ScalarTraits< _ST_ >::is_complex())
  {
    sprintf(hbfile,"%c%s%d.cua",tpc,filebase,nglob);
  } 
  else
  {
    sprintf(hbfile,"%c%s%d.rua",tpc,filebase,nglob);
  }

  sprintf(binfile,"%c%s%d.bin",tpc,filebase,nglob);
  sprintf(crsfile,"%c%s%d.crs",tpc,filebase,nglob);
  
  PHIST_SOUT(PHIST_DEBUG, "Looking for matrix \'%s\'..\n", filebase);
  *iflag=0;

  PHIST_SOUT(PHIST_DEBUG, "... try \'%s\'\n", binfile);
  SUBR(sparseMat_read_bin)(ptr,comm,binfile,&_iflag);
  if (_iflag==PHIST_SUCCESS) return;

  // try same format, different extension
  PHIST_SOUT(PHIST_DEBUG, "... try \'%s\'\n", crsfile);
  SUBR(sparseMat_read_bin)(ptr,comm,crsfile,&_iflag);
  if (_iflag==PHIST_SUCCESS) return;
  
  PHIST_SOUT(PHIST_DEBUG, "... try \'%s\'\n", mmfile);
  SUBR(sparseMat_read_mm)(ptr,comm,mmfile,&_iflag);
  if (_iflag==PHIST_SUCCESS) return;

  PHIST_SOUT(PHIST_DEBUG, "... try \'%s\'\n", hbfile);
  SUBR(sparseMat_read_hb)(ptr,comm,hbfile,&_iflag);
  *iflag=_iflag;
  return;
}//read_mat
