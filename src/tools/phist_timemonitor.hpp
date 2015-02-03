#ifndef PHIST_TIMEMONITOR_HPP
#define PHIST_TIMEMONITOR_HPP

#include "phist_config.h"

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#ifdef PHIST_HAVE_OPENMP
#include <omp.h>
#endif


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
      Timer(const char* s)
      {
#pragma omp barrier
#pragma omp master
        if( wtime_available() )
        {
          name = s;
          wtime = get_wtime();
        }
      }

      // stop timer
      ~Timer()
      {
#pragma omp barrier
#pragma omp master
        if( !name.empty() && wtime_available() )
        {
          wtime = get_wtime() - wtime;
          _timingResults[name].update(wtime);
        }
      }

      // calculates the results and prints them
      static void summarize(void);

    private:
      std::string name;
      double wtime;
      // simple struct that stores data
      struct TimeData
      {
        unsigned long numberOfCalls;
        double minTime;
        double maxTime;
        double totalTime;

        TimeData() :
          numberOfCalls(0),
          minTime(1.e100),
          maxTime(0.),
          totalTime(0.)
        {}

        void update(double time)
        {
          numberOfCalls++;
          minTime = std::min(minTime, time);
          maxTime = std::max(maxTime, time);
          totalTime += time;
        }
      };


      typedef std::map<std::string,TimeData> TimeDataMap;
      static TimeDataMap _timingResults;


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
#ifdef PHIST_HAVE_MPI
        int mpi_ini;
        MPI_Initialized(&mpi_ini);
        if( !mpi_ini )
          return false;
        MPI_Finalized(&mpi_ini);
        if( mpi_ini )
          return false;
#endif

        return true;
      }

  };


}
#endif /* PHIST_USE_TEUCHOS_TIMEMONITOR */

#endif /* PHIST_TIMEMONITOR_HPP */
