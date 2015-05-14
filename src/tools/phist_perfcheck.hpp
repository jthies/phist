#ifndef PHIST_PERFCHECK_HPP
#define PHIST_PERFCHECK_HPP

#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif


#ifndef PHIST_PERFCHECK

#define PHIST_PERFCHECK_VERIFY(functionName, benchFormula)
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
#define PHIST_PERFCHECK_VERIFY(functionName,n1,n2,n3, benchFormula) \
  phist_PerfCheck::PerfCheckTimer YouCanOnlyHaveOnePerfCheckInOneScope(functionName,n1,n2,n3, #benchFormula, benchFormula);


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
    benchFunction(&benchResult,&ierr); \
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
  BenchmarkMap benchmarks;

  //! Performance checking timer
  class PerfCheckTimer : public Timer
  {
    public:
      PerfCheckTimer(const char* name, double n1, double n2, double n3, const char* formula, double expectedTime) :
        Timer(constructName(name,n1,n2,n3,formula).c_str())
      {
        expectedResults_[name_].update(expectedTime);
      }

      // generate and print statistics
      static void summarize(int verbosity);

    protected:
      static Timer::TimeDataMap expectedResults_;

      static std::string constructName(const char* name, double n1, double n2, double n3, const char* formula)
      {
        std::ostringstream oss;
        oss << name << "(" << n1 << "," << n2 << "," << n3 << ") " << formula;
        return oss.str();
      }
  };
}

#endif /* PHIST_PERFCHECK */

#endif /* PHIST_PERFCHECK_HPP */
