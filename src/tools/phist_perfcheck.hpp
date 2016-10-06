#ifndef PHIST_PERFCHECK_HPP
#define PHIST_PERFCHECK_HPP

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#ifdef PHIST_KERNEL_LIB_GHOST
#include <ghost.h>
#endif

#ifndef PHIST_PERFCHECK

#define PHIST_PERFCHECK_VERIFY(functionName, benchFormula, flops)
#define PHIST_PERFCHECK_BENCHMARK(benchName, benchFunction)
#define PHIST_PERFCHECK_SUMMARIZE(verbosity)

#else

#include "phist_timemonitor.hpp"
#include <map>
#include <sstream>
#include <string>


/*! Defines a performance check, e.g. measure time until scope is left and compare it to the expected time.
 * Results (and especially large deviations) can be reported later.
 * \param functionName name of the function, ideally with required block sizes appended
 * \param benchFormula formula to calculate expected time
 */
#define PHIST_PERFCHECK_VERIFY(functionName,n1,n2,n3,n4,n5,n6,n7, benchFormula,flops) \
  phist_PerfCheck::PerfCheckTimer YouCanOnlyHaveOnePerfCheckInOneScope(functionName,#n1,n1,#n2,n2,#n3,n3,#n4,n4,#n5,n5,#n6,n6,#n7,n7, #benchFormula, benchFormula,flops);


/*! Defines a new benchmark for the performance check
 * \param benchName the name of the benchmark
 * \param benchFunction the benchmark function, should be "void fun(double* benchVal, int*ierr)"
 */
#define PHIST_PERFCHECK_BENCHMARK(benchName, benchFunction) \
namespace { \
double benchName(double factor) { \
  /* calculate benchmark value only if not done yet */ \
  double benchResult = phist_PerfCheck::benchmarks[__FUNCTION__]; \
  if( benchResult == 0 ) { \
    int ierr = 0; \
    double mean_bw; \
    benchFunction(&mean_bw,&benchResult,&ierr); \
    phist_PerfCheck::benchmarks[__FUNCTION__] = benchResult; \
  } \
  return factor / benchResult; \
} \
}

/*! Print gathered data from performance checks
 */
#define PHIST_PERFCHECK_SUMMARIZE(verbosity) phist_PerfCheck::PerfCheckTimer::summarize(verbosity);


//! namespace for performance check stuff
namespace phist_PerfCheck
{
  using phist_TimeMonitor::Timer;

  //! Map for benchmark results
  typedef std::map<std::string,double> BenchmarkMap;

  //! Lookup table for performance benchmark data
  static BenchmarkMap benchmarks;

  //! Performance checking timer
  class PerfCheckTimer : public Timer
  {
    public:
      PerfCheckTimer(const char* name, const char* sn1, double n1, 
                                       const char* sn2, double n2, 
                                       const char* sn3, double n3, 
                                       const char* sn4, double n4, 
                                       const char* sn5, double n5, 
                                       const char* sn6, double n6, 
                                       const char* sn7, double n7, 
                                       const char* formula, double expectedTime, double flops) :
        Timer(constructName(name,sn1,n1,sn2,n2,sn3,n3,sn4,n4,sn5,n5,sn6,n6,sn7,n7,formula).c_str())
      {
        double peak=PHIST_PEAK_FLOPS_CPU;
#ifdef PHIST_KERNEL_LIB_GHOST
        ghost_type gtype;
        ghost_type_get(&gtype);
        if (gtype==GHOST_TYPE_CUDA) peak=PHIST_PEAK_FLOPS_GPU;
#endif
        double t_peak = flops/peak;
        expectedResults_[name_].update(std::max(expectedTime,t_peak));
      }

      // generate and print statistics
      static void summarize(int verbosity);

    protected:
      static Timer::TimeDataMap expectedResults_;

      static std::string constructName(const char* name, 
                                       const char* sn1, double n1, 
                                       const char* sn2, double n2, 
                                       const char* sn3, double n3, 
                                       const char* sn4, double n4, 
                                       const char* sn5, double n5, 
                                       const char* sn6, double n6, 
                                       const char* sn7, double n7, 
                                       const char* formula)
      {
        std::ostringstream oss;

        oss << name << "(";
        if( sn1 != std::string("0") )
          oss << sn1 << "=" << n1;
        if( sn2 != std::string("0") )
          oss << "," << sn2 << "=" << n2;
        if( sn3 != std::string("0") )
          oss << "," << sn3 << "=" << n3;
        if( sn4 != std::string("0") )
          oss << "," << sn4 << "=" << n4;
        if( sn5 != std::string("0") )
          oss << "," << sn5 << "=" << n5;
        if( sn6 != std::string("0") )
          oss << "," << sn6 << "=" << n6;
        if( sn7 != std::string("0") )
          oss << "," << sn7 << "=" << n7;
        oss << ") " << formula;

        return oss.str();
      }
  };
}

#endif /* PHIST_PERFCHECK */

#endif /* PHIST_PERFCHECK_HPP */
