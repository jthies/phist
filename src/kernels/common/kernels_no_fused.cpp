// this file can be included if a kernel library does not implement
// the rather uncommon 'fused kernels', which compute several independent
// operations at once to exploit memory locality.
#include "kernels_no_fused_mvec.cpp"
#include "kernels_no_fused_spmv.cpp"
