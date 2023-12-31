#include <string>
#include <iostream>
#include <sstream>
#include <stdarg.h> 
#include <cstdlib>

#include <mpi.h>

#include "include/cpHelperFuncs.hpp"
#include "include/timing.h"

char * craftLogFile = "craftLog";

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
  if ((craftDebug > 0 )) {
    if(mpiCommWorldRank() == craftDebugProc){
      fprintf(stderr, "%d: CRAFT_ERROR: ", mpiCommWorldRank());
      va_start(argp, fmt);
      vfprintf(stderr, fmt, argp);
      va_end(argp);
      fprintf(stderr, "\n");
    }
  }
}

void craftDbg(int level, const char *fmt, ...)
{
  va_list argp;
  if ((craftDebug > 0 && craftDebug >= level)) {
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

int craftLog(MPI_Comm const * comm, const char *fmt, ...)
{
  int myrank=-1;
  MPI_Comm_rank(*comm, &myrank);
  if(myrank == craftDebugProc){
    if (( craftDebug > 0 )) {
      FILE * fstrL;
      fstrL = fopen(craftLogFile, "a");
      // ===== WRITE RESCUE LIST ===== // 
      va_list argp;
      va_start(argp, fmt);
      vfprintf(fstrL, fmt, argp);
      va_end(argp);
      fclose(fstrL); 
    }
  }
  return EXIT_SUCCESS;
}

void craftTime(const std::string str)
{ 
//  if(craftTiming){
  if(1){
    double t=0.0;
    get_walltime_ (&t);
    va_list argp;
    MPI_Comm parent;
	  MPI_Comm_get_parent( &parent );
    if(mpiCommWorldRank() == craftDebugProc ){ 
      printf("craftTime:: %s: %f\n", str.c_str(), t);
    }
  }
}

void craftTime(const std::string str, const MPI_Comm * const comm)
{ 
//  if(craftTiming){
  if(1){
    int myrank = -1;
    MPI_Comm_rank(*comm, &myrank);
    double t=0.0;
    get_walltime_ (&t);
    va_list argp;
    MPI_Comm parent;
	  MPI_Comm_get_parent( &parent );
    if(myrank == craftDebugProc ){      // only processes of original MPI_COMM_WORLD can do craft Time. This is to avoid having multiple entries after spawning
      printf("craftTime:: %s: %f\n", str.c_str(), t);
    }
  }
}



void getEnvVal(int &var, const char * str){
  int myrank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if(myrank==0) 
  {
    if(const char* env_p = std::getenv(str)){
      var = stringToNumber<int>(env_p);
      std::cout << str << " val is : " << var << '\n';
    }else{
      std::cout << "Env. val of " << str << " is not set.\n";
    }
  }
  return;
}

void checkDirectoryName(std::string* const s)
{
  std::string::iterator it;
  for (it = s->begin() ; it < s->end() ; ++it){
    switch(*it){
      case '/':case '\\':case ':':case '?':case '"':case '<':case '>':case '|':
      {
		    craftDbg(0, "Invalid character in directroy name: %s", s );
        craftAbort(0, "Invalid character in directroy name: %s", s );
      }
    }
  }
}
void checkPathName(std::string* const s)
{
  std::string::iterator it;
  for (it = s->begin() ; it < s->end() ; ++it){
    switch(*it){
      case '\\':case ':':case '?':case '"':case '<':case '>':case '|':
      {
		    craftDbg(0, "Invalid character in directroy name: %s", s );
        craftAbort(0, "Invalid character in directroy name: %s", s );
      }
    }
  }
}

std::string exec(const char* cmd){ 
  char buffer[128];
  std::string result = "";
  FILE* pipe = popen(cmd, "r");
//  if (!pipe) 
//				throw std::runtime_error("popen() failed!");
  try {
    while (!feof(pipe)) {
       if (fgets(buffer, 128, pipe) != NULL)
         result = buffer;
  	}
  } 
	catch (...) {
    pclose(pipe);
  	throw;
	}
	pclose(pipe);
	result.erase(result.end()-1);
	return result;
}

/*#define MAX_PATH_LEN 256

#define IN
#define OUT

typedef int Bool_T;

static
int GetCharsNumInPath( IN const char * path )
{
int cnt;

int i;
int n = strlen( path );

for( i = 0; i < n; i++ )
{
if( ( path[i] < 0x80 ) || ( path[i] > 0xbf ) )
cnt++;
}

return cnt;
}

Bool_T PathValid( IN const char * path )
{
  if( path != NULL )
  {
    if( GetCharsNumInPath( IN path ) <= MAX_PATH_LEN )
    {
      int i;
      int n = strlen( path );

      for( i = 0; i < n; i++ )
      {
        switch( path[i] )
        {
        // ? " / < > * |
        // these characters can not be used in file or folder names
        //
          case '?':
          case '\"':
          case '/':
          case '<':
          case '>':
          case '*':
          case '|':
              return false;
          
          // Can meet only between a local disk letter and full path
          // for example D:\folder\file.txt
          //
          case ':':
          {
              if( i != 1 )
              {
                return false;
              }
              else{
              break;
              }
          }
// Space and point can not be the last character of a file or folder names
//
          case ' ':
          case '.':
          {
            if( ( i + 1 == n ) || ( path[i+1] == PATH_SEPERATOR_CHAR ) )
            {
               return false;
            }
            else{
              break;
            }
          }
// two backslashes can not go straight
//
          case PATH_SEPERATOR_CHAR:
          {
            if( i > 0 && path[i - 1] == PATH_SEPERATOR_CHAR )
            {
              return false;
            }else{
              break;
            }
          }
        }
      }
    return true;
    }
    else{ // if( GetCharsNumInPath( IN path ) <= MAX_PATH_LEN )
      LOG_ERROR( "PathValid FAILURE --> path is too long" );
      return false;
    }
  }else{ // if( path != NULL )
// wrong argument
//
  return false;
  }
}  
*/
