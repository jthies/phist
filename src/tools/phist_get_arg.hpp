#ifndef PHIST_GETARG_HPP
#define PHIST_GETARG_HPP

#ifdef __cplusplus
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
    ss>>T;
  }
return result;
}

// macro intended for use in CXX main program.
// output argument should have a default value
#define GET_ARG(_a_,_p_) _a_=get_arg(argc,argv,_p_,_a_); \
std::cout << #_a_ " = "<<_a_<<std::endl;

#else
//TODO: provide something for C main programs,
// like get_darg, get_iarg etc.
#endif
#endif
