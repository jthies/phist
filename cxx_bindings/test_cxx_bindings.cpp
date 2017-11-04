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
  phist_gidx n=42;
  int iflag=0;
  phist_map_ptr map=NULL;
  phist_comm_ptr comm=NULL;
  phist_comm_create(&comm,&iflag);
  phist_map_create(&map,comm,n,&iflag);
  phist::types<ST>::mvec_ptr V=NULL;
  try {
    phist::kernels<ST>::mvec_create(&V,map,m,&iflag);
    phist::kernels<ST>::mvec_delete(&V,map,m,&iflag);
  } catch (phist::Exception e)
  {
    if (phist::kernels<ST>>::type_avail()==PHIST_NOT_IMPLEMENTED) return true;
    return false;
  }
}

int main(int argc, char* argv[])
{
  phist_kernels_init(&argc,&argv);
  if (
  phist_foo<float>() &&
  phist_foo<double>() &&
  phist_foo<phist_Scomplex>() &&
  phist_foo<phist_Scomplex>()) error=0;
  else error=-1;
  phist_kernels_finalize();
  return error;
}
