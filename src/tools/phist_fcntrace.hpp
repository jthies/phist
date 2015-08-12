#ifndef PHIST_FCN_TRACE_HPP
#define PHIST_FCN_TRACE_HPP

#ifndef DOXYGEN
#include "phist_config.h"

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include <string>
#include <vector>

#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif

#endif /* DOXYGEN */

// this is an object that prints a message when created and when destroyed.
// Its sole purpose is to implement a macro PHIST_ENTER_FCN (phist_macros.h) that
// will generate output like this:
// ENTER xyz
// ...
// LEAVE xyz
// for some functino xyz. The macro is defined empty of PHIST_OUTLEV<PHIST_TRACE
class phist_FcnTrace
{
  public:
    phist_FcnTrace(const char* fcn)
    {
#pragma omp master
      {
#ifdef PHIST_TIMINGS_FULL_TRACE
        fcnTrace_.push_back(fcn);
        fcn_ = getFullTraceString();
#else
        fcn_ = fcn;
#endif
      }

#if (PHIST_OUTLEV>=PHIST_TRACE)
      PHIST_OUT(0,"PHIST ENTER %s\n",fcn_.c_str());
#endif
#ifdef PHIST_HAVE_LIKWID
#pragma omp parallel
      {
        LIKWID_MARKER_START(fcn_.c_str());
      }
#endif
    }


    ~phist_FcnTrace()
    { 
#ifdef PHIST_HAVE_LIKWID
#pragma omp parallel
      {
        LIKWID_MARKER_STOP(fcn_.c_str());
      }
#endif

#if (PHIST_OUTLEV>=PHIST_TRACE)
      PHIST_OUT(0,"PHIST LEAVE %s\n",fcn_.c_str()); 
#endif

#ifdef PHIST_TIMINGS_FULL_TRACE
#pragma omp master
      {
        fcnTrace_.pop_back();
      }
#endif
    }


    const std::string& fcn() const {return fcn_;}

  private:
    std::string fcn_;


#ifdef PHIST_TIMINGS_FULL_TRACE
    static std::vector<const char*> fcnTrace_;


    // function that constructs a string from the fcnTrace_ stack
    std::string getFullTraceString() const
    {
      std::string s = "";
      for(std::vector<const char*>::const_iterator it = fcnTrace_.begin(); it != fcnTrace_.end(); it++)
      {
        s += " > ";
        s += *it;
      }
      return s;
    }
#endif
};



// Small helper class for PHIST_ENTER_KERNEL_FCN (phist_macros.h),
// which returns a modified name when nested
class phist_CheckKernelFcnNesting
{
  public:
    phist_CheckKernelFcnNesting(const char* fcn)
    {
      if( !nestedKernelCall_ )
        fcn_ = fcn;

      if( !fcn_.empty() )
        nestedKernelCall_ = true;
    }

    ~phist_CheckKernelFcnNesting()
    {
      if( !fcn_.empty() )
        nestedKernelCall_ = false;
    }

    const char* str() const {return fcn_.c_str();}

  private:
#ifdef PHIST_HAVE_CXX11_THREADLOCAL
    thread_local static bool nestedKernelCall_;
#else
    static bool nestedKernelCall_;
#endif
    std::string fcn_;

};

#endif

