/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/

#ifndef PHIST_DISALLOW_PRECON_TRAITS
#define PHIST_DISABLED_PRECON_TRAITS

// this macro allows us to make certain preconditioners for certain data types act as
// identity operators, we use it to e.g. disale Ifpack2<float> if the Tpetra installation
// does not support float as a data type.
#define PHIST_DISALLOW_PRECON_TRAITS(ST,PT) \
template<> \
class PreconTraits< ST , PT > \
{ \
\
  static void NotImplemented(int* iflag) \
  { \
    PHIST_SOUT(PHIST_ERROR,"This preconditioner/dat type combination is explicitly disabled in your PHIST installation. \n" \
                           "The most likely cause is that the installation of the kernel library does not allow us to\n" \
                           "instantiate the preconditioner/data type combination.");\
    *iflag=PHIST_NOT_IMPLEMENTED;\
  }\
\
public: \
\
  static void Usage()\
  { \
    int iflag=-0; NotImplemented(&iflag);\
  }\
\
  static void Create(void** P, \
        const void* A, ST sigma, const void* B, \
        const void* Vkern, const void* BVkern,\
        const char* options, void* last_arg, int* iflag)\
  {\
    NotImplemented(iflag);\
    return;\
  }\
\
  static void Wrap(void** P, \
        const void* A, ST sigma, const void* B, \
        const void* Vkern, const void* BVkern,\
        void* last_arg, int* iflag)\
  {\
    NotImplemented(iflag);\
    return;\
  }\
\
  static void Update(void* P, const void* A, ST sigma, const void* B,\
        const void* Vkern, const void* BVkern,\
        int* iflag)\
  {\
    NotImplemented(iflag);\
    return;\
  }\
  static void Delete(void* P, int* iflag)\
  {\
    NotImplemented(iflag);\
    return;\
  }\
  static void Apply(ST alpha, void const* P, void const* X, ST beta, void* Y, int* iflag)\
  {\
    NotImplemented(iflag);\
    return;\
  }\
  static void ApplyT(ST alpha, void const* P, void const* X, ST beta, void* Y, int* iflag)\
  {\
    NotImplemented(iflag);\
    return;\
  }\
  static void ApplyShifted(ST alpha, const void* P, ST const * sigma,\
          void const* X, ST beta,  void* Y, int* iflag)\
  {\
    NotImplemented(iflag);\
    return;\
  }\
};
#endif
