void SUBR(read_mat)(const char* filebase,int nglob,TYPE(crsMat_ptr) *ptr, int* ierr)
{
  int _ierr;
  *ptr = NULL;
  comm_ptr_t comm;
  phist_comm_create(&comm,&_ierr);
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
  PHIST_SOUT(PHIST_DEBUG, "... try \'%s\'\n", mmfile);
  SUBR(crsMat_read_mm)(ptr,comm,mmfile,&_ierr);
  if (_ierr!=PHIST_SUCCESS) // kernel lib can't read MatrixMarket format or file not found
  {
    PHIST_SOUT(PHIST_DEBUG, "... try \'%s\'\n", hbfile);
    SUBR(crsMat_read_hb)(ptr,comm,hbfile,&_ierr);
    if (_ierr!=PHIST_SUCCESS) // kernel lib can't read Harwell-Boeing or file not found
    {
      PHIST_SOUT(PHIST_DEBUG, "... try \'%s\'\n", binfile);
      SUBR(crsMat_read_bin)(ptr,comm,binfile,&_ierr);
      if (_ierr!=PHIST_SUCCESS) // kernel lib can't read binCRS or file not found
      {
        // try same format, different extension
        PHIST_SOUT(PHIST_DEBUG, "... try \'%s\'\n", crsfile);
        SUBR(crsMat_read_bin)(ptr,comm,crsfile,&_ierr);
      }
    }
  }
*ierr=_ierr;
return;
}//read_mat
