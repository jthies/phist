#ifndef PHIST_TIMEMONITOR_HPP
#define PHIST_TIMEMONITOR_HPP

#include "phist_config.h"

/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#ifndef PHIST_HAVE_MPI
#error "phist_TimeMonitor does only work with MPI_Wtime(), but no MPI found"
#endif

#include <mpi.h>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>


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
        if( !mpi_wtime_available() )
          return;
#pragma omp master
{
        fcnName = s;
        wtime = MPI_Wtime();
}
      }

      // stop timer
      ~Timer()
      {
#pragma omp barrier
        if( fcnName.empty() || !mpi_wtime_available() ) // not ready yet
          return;
#pragma omp master
{
        wtime = MPI_Wtime() - wtime;
        //PHIST_OUT(PHIST_DEBUG, "Measured %10.4e wtime for function %s\n", wtime, fcnName.c_str());
        _timingResults[fcnName].update(wtime);
}
      }

    private:
      std::string fcnName;
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


      // functor for sorting by keys in another vector
      struct SortClass
      {
        SortClass(const std::vector<double>& keys) : _keys(keys) {}
        bool operator() (int i, int j) {return _keys[i] > _keys[j];}
        const std::vector<double>& _keys;
      };

      // helper function to determine if mpi is initialized
      static bool mpi_wtime_available()
      {
        int mpi_ini;
        MPI_Initialized(&mpi_ini);
        if( !mpi_ini )
          return false;
        MPI_Finalized(&mpi_ini);
        if( mpi_ini )
          return false;

        return true;
      }
    public:

      // calculates the results and prints them
      static void summarize(void)
      {
        int nTimers = _timingResults.size();

        // consider only timers from the first process!
        MPI_Bcast(&nTimers, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if( nTimers == 0 )
          return;

        // get timer names from the first process
        std::string fcnNameList;
        for(TimeDataMap::const_iterator it = _timingResults.begin(); it != _timingResults.end(); it++)
          fcnNameList.append( it->first + '\n' );
        int strLen = fcnNameList.length();
        MPI_Bcast(&strLen, 1, MPI_INT, 0, MPI_COMM_WORLD);
        char* strBuf = new char[strLen];
        fcnNameList.copy(strBuf, strLen);
        strBuf[strLen-1] = '\0';
        MPI_Bcast(strBuf, strLen, MPI_CHAR, 0, MPI_COMM_WORLD);
        std::istringstream iss(strBuf);

        std::vector<std::string> fcnName(nTimers);
        std::vector<unsigned long> numberOfCalls(nTimers);
        std::vector<double> minTime(nTimers);
        std::vector<double> maxTime(nTimers);
        std::vector<double> maxTotalTime(nTimers);
        std::vector<double> sumTotalTime(nTimers);

        std::vector<double> tmp(nTimers);

        // convert / copy everything into vectors
        for(int i = 0; i < nTimers; i++)
        {
          iss >> fcnName.at(i);
          numberOfCalls.at(i) = _timingResults[fcnName.at(i)].numberOfCalls;
          minTime.at(i) = _timingResults[fcnName.at(i)].minTime;
          maxTime.at(i) = _timingResults[fcnName.at(i)].maxTime;
          maxTotalTime.at(i) = _timingResults[fcnName.at(i)].totalTime;
        }

        // reductions over mpi processes
        tmp = minTime;
        MPI_Reduce(&tmp[0],&minTime[0],nTimers,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
        tmp = maxTime;
        MPI_Reduce(&tmp[0],&maxTime[0],nTimers,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
        tmp = maxTotalTime;
        MPI_Reduce(&tmp[0],&sumTotalTime[0],nTimers,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        MPI_Reduce(&tmp[0],&maxTotalTime[0],nTimers,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

        // sort by total time
        std::vector<int> sortedIndex(nTimers);
        for(int i = 0; i < nTimers; i++)
          sortedIndex.at(i) = i;
        SortClass sortFunc(maxTotalTime);
        std::sort(sortedIndex.begin(), sortedIndex.end(), sortFunc);

        // print result on proc 0
        PHIST_SOUT(PHIST_INFO, "======================================== TIMING RESULTS ============================================\n");
        PHIST_SOUT(PHIST_INFO, "%40s  %10s  %10s  %10s  %10s  %10s\n", "function", "mtot.time", "count", "max.time", "avg.time", "min.time");
        int nprocs;
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        for(int i_ = 0; i_ < nTimers; i_++)
        {
          int i = sortedIndex.at(i_);
          if( i_ % 10 == 0 )
          {
            PHIST_SOUT(PHIST_INFO, "----------------------------------------------------------------------------------------------------\n");
          }
          PHIST_SOUT(PHIST_INFO, "%40s  %10.3e  %10lu  %10.3e  %10.3e  %10.3e\n", fcnName.at(i).c_str(), maxTotalTime.at(i), numberOfCalls.at(i), maxTime.at(i), sumTotalTime.at(i)/numberOfCalls.at(i)/nprocs, minTime.at(i));
        }
        PHIST_SOUT(PHIST_INFO, "====================================================================================================\n");
      }

  };


}

#endif /* PHIST_TIMEMONITOR_HPP */
