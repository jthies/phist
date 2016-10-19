#include <string>
#include <iostream>
#include <sstream>
#include <stdarg.h> 
#include <cstdlib>

#include <mpi.h>

#include "include/cpHelperFuncs.hpp"

template std::string numberToString<int>(int number);
template std::string numberToString<double>(double number);
template std::string numberToString<float>(float number);

template int stringToNumber<int>(const std::string& text);
template double stringToNumber<double>(const std::string& text);
template float stringToNumber<float>(const std::string& text);


template <typename T>
std::string numberToString ( T number )
{
	std::ostringstream ss;
	ss << number;
	return ss.str();
}

template <typename T>
T stringToNumber ( const std::string &text )
{
  std::istringstream ss(text);
  T result;
  return ss >> result ? result : 0;
}

inline int mpiCommWorldRank(){
  int myRank=-1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);  
  return myRank;
}
inline int mpiCommWorldNumProcs(){
  int numProcs=-1;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);  
  return numProcs;
}

void craftErr(const char *fmt, ...)
{
  va_list argp;
  fprintf(stderr, "%d: CRAFT_ERROR: ", mpiCommWorldRank());
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);
  fprintf(stderr, "\n");
}

void craftDbg(int level, const char *fmt, ...)
{
  va_list argp;
  static int craftDebug, craftDebugProc;
  static int t=true;
  if(t==true){
    getEnvVal(craftDebug    , "CRAFT_DEBUG");
    getEnvVal(craftDebugProc, "CRAFT_DEBUG_PROC");
    t = false;
  }
  if (level == 0 || (craftDebug > 0 && craftDebug >= level)) {
    if(mpiCommWorldRank() == craftDebugProc){
      fprintf(stdout, "%d: CRAFT_DEBUG :", mpiCommWorldRank());
      va_start(argp, fmt);
      vfprintf(stdout, fmt, argp);
      va_end(argp);
      fprintf(stdout, "\n");
    }
  }
}

void craftAbort(int rc, const char *fmt, ...)
{
  va_list argp;
  fprintf(stderr, "%d: CRAFT_ABORT: ", mpiCommWorldRank());
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);
  fprintf(stderr, "\n");

  MPI_Abort(MPI_COMM_WORLD, 0);
}


void getEnvVal(int &var, const char * str){
  if(const char* env_p = std::getenv(str)){
    var = stringToNumber<int>(env_p);
    std::cout << str << " val is : " << var << '\n';
  }else{
    std::cout << "Env. val of " << str << " is not set.\n";
  }
  return;
}


