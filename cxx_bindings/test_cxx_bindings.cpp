/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/

#include "phist_config.h"
#include "phist_types.hpp"
#include "phist_kernels.hpp"
#include "phist_core.hpp"
#include "phist_krylov.hpp"
#include "phist_precon.hpp"
#include "phist_jada.hpp"

template<typename ST>
bool phist_foo()
{
  int type_avail=0;
  phist::kernels<ST>::type_avail(&type_avail);
  phist_gidx n=42;
  int m=5;
  int iflag=0;
  phist_map_ptr map=NULL;
  phist_comm_ptr comm=NULL;
  phist_comm_create(&comm,&iflag);
  phist_map_create(&map,comm,n,&iflag);
  typename phist::types<ST>::mvec_ptr V=NULL;
  try {
    phist::kernels<ST>::mvec_create(&V,map,m,&iflag);
    phist::kernels<ST>::mvec_delete(&V,map,m,&iflag);
  } catch (phist::Exception e)
  {
  int iflag=e.iflag();
    if (type_avail==PHIST_NOT_IMPLEMENTED &&iflag==PHIST_NOT_IMPLEMENTED) return true; // correct behavior
    std::cerr << "caught phist exception: '" << e.what() << "'"<<std::endl
              << "(phist iflag="<<iflag<<")"<<std::endl;
    return false;
  }
}

int main(int argc, char* argv[])
{
  int iflag=0;
  phist_kernels_init(&argc,&argv,&iflag);

  int error=0;

  if (!phist_foo<double>()) error|=1;
  if (!phist_foo<phist_d_complex>()) error|=2;

#ifdef PHIST_HAVE_SP
  if (!phist_foo<float>()) error|=4;
  if (!phist_foo<phist_s_complex>()) error|=8;
#endif

  phist_kernels_finalize(&iflag);
  return -error;
}
