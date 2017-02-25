#ifndef PHIST_TIMEMONITOR_HPP
#define PHIST_TIMEMONITOR_HPP

#include "phist_config.h"
#ifdef PHIST_HAVE_OPENMP
#include <omp.h>
#endif

extern "C" const char* phist_retcode2str(int code);

#ifdef PHIST_USE_TEUCHOS_TIMEMONITOR

// use the Trilinos implementation of a TimeMonitor

// the Teuchos time monitor is a bit fancier when it comes
// to non bulk-synchronous execution models like master/slave
# include "Teuchos_TimeMonitor.hpp"
# include "Teuchos_DefaultComm.hpp"

# define PHIST_CXX_TIMER(s) static Teuchos::RCP<Teuchos::Time> TimeFrom_PHIST_CXX_TIMER; \
                            if( TimeFrom_PHIST_CXX_TIMER.is_null() ) \
                                TimeFrom_PHIST_CXX_TIMER = Teuchos::TimeMonitor::getNewTimer(s); \
                            Teuchos::TimeMonitor TimeMonFrom_PHIST_CXX_TIMER(*TimeFrom_PHIST_CXX_TIMER);
# if TRILINOS_MAJOR_MINOR_VERSION < 110800
#  define PHIST_CXX_TIMER_SUMMARIZE \
                                    Teuchos::TimeMonitor::summarize(Teuchos::DefaultComm<int>::getComm().ptr(), \
                                    std::cout,true,true,false, Teuchos::Union,"");
# else
#  define PHIST_CXX_TIMER_SUMMARIZE \
                                    Teuchos::TimeMonitor::summarize(Teuchos::DefaultComm<int>::getComm().ptr(), \
                                    std::cout,true,true,false, Teuchos::Union,"",true);
# endif

#else

// use our own implementation of a TimeMonitor
//
# define PHIST_CXX_TIMER(s) phist_TimeMonitor::Timer TimerFrom_PHIST_CXX_TIMER(s);
# define PHIST_CXX_TIMER_SUMMARIZE phist_TimeMonitor::Timer::summarize();

#endif /* PHIST_USE_TEUCHOS_TIMEMONITOR */

# include <map>
# include <string>


// encapsulates functions and classes needed for function level time monitoring
namespace phist_TimeMonitor
{
  // used to time a function
  class Timer
  {
    public:
      // start timer
      Timer(const char* s);

      // stop timer
      ~Timer();

      // calculates the results and prints them
      static void summarize(void);

    protected:
      std::string name_;
      double wtime_;
      // simple struct that stores data
      struct TimeData
      {
        unsigned long numberOfCalls;
        double minTime;
        double maxTime;
#ifdef PHIST_TIMINGS_FULL_TRACE
        double selfTime;
#endif
        double totalTime;

        TimeData() :
          numberOfCalls(0),
          minTime(1.e100),
          maxTime(0.),
#ifdef PHIST_TIMINGS_FULL_TRACE
          selfTime(0.),
#endif
          totalTime(0.)
        {}

        void update(double time)
        {
          numberOfCalls++;
          minTime = std::min(minTime, time);
          maxTime = std::max(maxTime, time);
          totalTime += time;
#ifdef PHIST_TIMINGS_FULL_TRACE
          selfTime += time;
#endif
        }

#ifdef PHIST_TIMINGS_FULL_TRACE
        void updateChild(double time)
        {
          selfTime -= time;
        }
#endif
      };


      typedef std::map<std::string,TimeData> TimeDataMap;
      static TimeDataMap timingResults_;


      static double get_wtime()
      {
#if defined(PHIST_HAVE_OPENMP)
        return omp_get_wtime();
#elif defined(PHIST_HAVE_MPI)
        return MPI_Wtime();
#else
#error "phist_TimeMonitor needs OpenMP or MPI"
#endif
      }

      // helper function to determine if mpi is initialized
      static bool wtime_available()
      {
#if defined(PHIST_HAVE_OPENMP)
        return true;
#elif defined(PHIST_HAVE_MPI)
        int mpi_ini, mpi_fini;
        MPI_Initialized(&mpi_ini);
        if( !mpi_ini )
          return false;
        MPI_Finalized(&mpi_fini);
        if( mpi_fini )
          return false;
        return true;
#else
        return false;
#endif
      }

#ifdef PHIST_TIMINGS_FULL_TRACE
      static std::string getNameOfParent(const std::string& child)
      {
        size_t len = child.find_last_of(" > ");
        return child.substr(0, len-2);
      }
#endif
  };


}

#endif /* PHIST_TIMEMONITOR_HPP */
