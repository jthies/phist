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

