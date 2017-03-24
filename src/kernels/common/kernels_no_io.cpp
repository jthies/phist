/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
void SUBR(mvec_write_bin)(TYPE(const_mvec_ptr) V, const char* filename, int* iflag)
{
  *iflag=-99;
}

void SUBR(mvec_read_bin)(TYPE(mvec_ptr) M, const char* filename, int* iflag)
{
  *iflag=-99;
}

void SUBR(sdMat_write_bin)(TYPE(const_sdMat_ptr) M, const char* filename, int* iflag)
{
  *iflag=-99;
}

void SUBR(sdMat_read_bin)(TYPE(sdMat_ptr) M, const char* filename, int* iflag)
{
  *iflag=-99;
}
