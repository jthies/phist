/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

// "fill" an mvec from a user-provided array.
void SUBR(mvec_get_data)(TYPE(const_mvec_ptr) V, _ST_* data_out, phist_lidx lda_out, int output_row_major, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  phist_lidx lda_in;
  phist_lidx lnrows;
  int ncols;
  _ST_* V_raw;
#ifdef PHIST_MVECS_ROW_MAJOR
  int input_row_major=1;
#else
  int input_row_major=0;
#endif
  // note: we allow a device sync here, strictly speaking V is const, but we trick it here
  TYPE(mvec_ptr) V_nonconst=const_cast<TYPE(mvec_ptr)>(V);
  PHIST_CHK_IERR(SUBR(mvec_from_device)(V_nonconst, iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&ncols,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&lnrows,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_extract_view)(V_nonconst,&V_raw, &lda_in, iflag), *iflag);
  for (phist_lidx i=0; i<lnrows; i++)
  {
    for (int j=0; j<ncols; j++)
    {
      phist_lidx idx_in  =  input_row_major? i*lda_in +j: j*lda_in +i;
      phist_lidx idx_out = output_row_major? i*lda_out+j: j*lda_out+i;
      data_out[idx_out]=V_raw[idx_in];
    }
  }
}
