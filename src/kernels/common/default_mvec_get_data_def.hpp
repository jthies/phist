// "fill" an mvec from a user-provided array.
void SUBR(mvec_get_data)(TYPE(const_mvec_ptr) V, _ST_* data_out, phist_lidx lda_out, int output_row_major, int* iflag)
{
#include "phist_std_typedefs.hpp"
  PHIST_ENTER_FCN(__FUNCTION__);
  phist_lidx lda;
  phist_lidx lnrows;
  int ncols;
  _ST_* V_raw;
#ifdef PHIST_MVECS_ROW_MAJOR
  int input_row_major=1;
#else
  int input_row_major=0;
#endif
  PHIST_CHK_IERR(SUBR(mvec_from_device)(V,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_num_vectors)(V,&ncols,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_my_length)(V,&lnrows,iflag),*iflag);
  PHIST_CHK_IERR(SUBR(mvec_extract_view)(V,&V_raw, &lda, iflag)
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


