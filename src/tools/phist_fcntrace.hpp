#ifndef PHIST_FCN_TRACE_HPP
#define PHIST_FCN_TRACE_HPP

#ifndef NO_INCLUDES_IN_HEADERS
#include <string>
#ifdef PHIST_HAVE_LIKWID
#include <likwid.h>
#endif
#endif

// this is an object that prints a message when created and when destroyed.
// Its sole purpose is to implement a macro ENTER_FCN (phist_macros.h) that
// will generate output like this:
// ENTER xyz
// ...
// LEAVE xyz
// for some functino xyz. The macro is defined empty of PHIST_OUTLEV<PHIST_TRACE
class FcnTracer
  {
  public:
  
  FcnTracer(const char* fcn) : fcn_(fcn)
    {
#if (PHIST_OUTLEV>=PHIST_TRACE)
    PHIST_OUT(0,"PHIST ENTER %s\n",fcn_.c_str());
#endif
    LIKWID_MARKER_START(fcn_.c_str());
    }

  ~FcnTracer()
    { 
    LIKWID_MARKER_STOP(fcn_.c_str());
#if (PHIST_OUTLEV>=PHIST_TRACE)
    PHIST_OUT(0,"PHIST LEAVE %s\n",fcn_.c_str()); 
#endif
    }

  std::string fcn_;
  };

#endif

