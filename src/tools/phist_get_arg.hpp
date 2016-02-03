#ifndef PHIST_GETARG_HPP
#define PHIST_GETARG_HPP

#ifdef __cplusplus

#include <iostream>

// command-line parsing tool for C++ main programs
// if the given position is smaller than argc, parse
// the argument argv[pos], otherwise use the default value
// given.
template<typename T>
T get_arg(int argc, char** argv, int pos, T& default_val)
{
  T result;
  if (pos>=argc)
  {
    result=default_val;
  }
  else
  {
    std::stringstream ss(argv[pos]);;
    ss>>result;
  }
return result;
}

// macro intended for use in CXX main program.
// output argument should have a default value
#define GET_ARG(_a_,_p_,_valid_) _a_=get_arg(argc,argv,_p_,_a_); \
{\
std::stringstream ss; ss << #_a_ " = "<<_a_; \
PHIST_SOUT(PHIST_VERBOSE,"%s\n",ss.str().c_str()); \
if (!(_valid_)) \
{\
  PHIST_SOUT(PHIST_ERROR,"parameter %d had an illegal value!\n" \
                         "Typically, phist drivers will print a usage message when run without arguments.",p);\
  return PHIST_INVALID_INPUT;\
}\
}

#else
//TODO: provide something for C main programs,
// like get_darg, get_iarg etc.
#error "GET_ARG only works in C++ code"
#endif
#endif
