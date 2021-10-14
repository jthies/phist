/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
// this file can be included if a kernel library does not implement
// the rather uncommon 'fused kernels', which compute several independent
// operations at once to exploit memory locality.
#include "kernels_no_fused_mvec.cpp"
#include "kernels_no_fused_spmv.cpp"
