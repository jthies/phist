/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
extern "C" void SUBR(mvec_to_device)(TYPE(mvec_ptr) V, int* iflag)
{
  *iflag=0;
}

extern "C" void SUBR(mvec_from_device)(TYPE(mvec_ptr) V, int* iflag)
{
  *iflag=0;
}

extern "C" void SUBR(sdMat_to_device)(TYPE(sdMat_ptr) M, int* iflag)
{
  *iflag=0;
}

extern "C" void SUBR(sdMat_from_device)(TYPE(sdMat_ptr) M, int* iflag)
{
  *iflag=0;
}

