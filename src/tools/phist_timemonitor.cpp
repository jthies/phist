/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#ifdef PHIST_KERNEL_LIB_GHOST
//#include <ghost/config.h>
#include <ghost.h>
# ifdef GHOST_HAVE_CUDA
# include <ghost/cu_util.h>
# endif
#endif
#include "phist_timemonitor.hpp"
#include "phist_macros.h"

#include <algorithm>
#include <sstream>


namespace phist_TimeMonitor
{

      Timer::Timer(const char *s)
      {
        if( wtime_available() )
        {
          name_ = s;
          wtime_ = get_wtime();
        }
      }

      Timer::~Timer()
      {
        if( !name_.empty() && wtime_available() )
        {
#ifdef PHIST_KERNEL_LIB_GHOST
# ifdef GHOST_HAVE_CUDA
  ghost_type gtype;
  ghost_type_get(&gtype);
  if (gtype==GHOST_TYPE_CUDA)
  {
    ghost_cu_barrier();
  }
# endif
#endif
          wtime_ = get_wtime() - wtime_;
#pragma omp critical (phist_timemonitor)
          {
            timingResults_[name_].update(wtime_);
#ifdef PHIST_TIMINGS_FULL_TRACE
            std::string parent = getNameOfParent(name_);
            if( !parent.empty() )
              timingResults_[parent].updateChild(wtime_);
#endif
          }
        }
      }

  Timer::TimeDataMap Timer::timingResults_;

  // functor for sorting by keys in another vector
  struct SortClass
  {
    SortClass(const std::vector<double>& keys) : _keys(keys) {}
    bool operator() (int i, int j) {return _keys[i] > _keys[j];}
    const std::vector<double>& _keys;
  };


  // calculates the results and prints them
  void Timer::summarize(void)
  {
    int ierr = 0;
    int nTimers = timingResults_.size();

    // consider only timers from the first process!
#ifdef PHIST_HAVE_MPI
    PHIST_CHK_IERR(ierr = MPI_Bcast(&nTimers, 1, MPI_INT, 0, MPI_COMM_WORLD), ierr);
#endif

    if( nTimers == 0 )
      return;

    // get timer names from the first process
    std::string fcnNameList;
    for(TimeDataMap::const_iterator it = timingResults_.begin(); it != timingResults_.end(); it++)
      fcnNameList.append( it->first + '\n' );
    int strLen = fcnNameList.length();
#ifdef PHIST_HAVE_MPI
    PHIST_CHK_IERR(ierr = MPI_Bcast(&strLen, 1, MPI_INT, 0, MPI_COMM_WORLD), ierr);
#endif
    char* strBuf = new char[strLen];
    fcnNameList.copy(strBuf, strLen);
    strBuf[strLen-1] = '\0';
#ifdef PHIST_HAVE_MPI
    PHIST_CHK_IERR(ierr = MPI_Bcast(strBuf, strLen, MPI_CHAR, 0, MPI_COMM_WORLD), ierr);
#endif
    std::istringstream iss(strBuf);

    std::vector<std::string> fcnName(nTimers);
    std::vector<unsigned long> numberOfCalls(nTimers);
    std::vector<double> minTime(nTimers);
    std::vector<double> maxTime(nTimers);
    std::vector<double> maxTotalTime(nTimers);
#ifdef PHIST_TIMINGS_FULL_TRACE
    std::vector<double> maxSelfTime(nTimers);
#endif
    std::vector<double> sumTotalTime(nTimers);

    std::vector<double> tmp(nTimers);

    // convert / copy everything into vectors
    for(int i = 0; i < nTimers; i++)
    {
      std::getline(iss, fcnName.at(i));
      numberOfCalls.at(i) = timingResults_[fcnName.at(i)].numberOfCalls;
      minTime.at(i) = timingResults_[fcnName.at(i)].minTime;
      maxTime.at(i) = timingResults_[fcnName.at(i)].maxTime;
      maxTotalTime.at(i) = timingResults_[fcnName.at(i)].totalTime;
#ifdef PHIST_TIMINGS_FULL_TRACE
      maxSelfTime.at(i) = timingResults_[fcnName.at(i)].selfTime;
#endif
    }
    delete[] strBuf;

#ifdef PHIST_HAVE_MPI
    // reductions over mpi processes
    tmp = minTime;
    PHIST_CHK_IERR(ierr = MPI_Reduce(&tmp[0],&minTime[0],nTimers,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD), ierr);
    tmp = maxTime;
    PHIST_CHK_IERR(ierr = MPI_Reduce(&tmp[0],&maxTime[0],nTimers,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD), ierr);
    tmp = maxTotalTime;
    PHIST_CHK_IERR(ierr = MPI_Reduce(&tmp[0],&sumTotalTime[0],nTimers,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD), ierr);
    PHIST_CHK_IERR(ierr = MPI_Reduce(&tmp[0],&maxTotalTime[0],nTimers,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD), ierr);
#ifdef PHIST_TIMINGS_FULL_TRACE
    tmp = maxSelfTime;
    MPI_Reduce(&tmp[0],&maxSelfTime[0],nTimers,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
#endif
#endif

    // sort by total time
    std::vector<int> sortedIndex(nTimers);
    for(int i = 0; i < nTimers; i++)
      sortedIndex.at(i) = i;
    SortClass sortFunc(maxTotalTime);
    std::sort(sortedIndex.begin(), sortedIndex.end(), sortFunc);

    // make all fcnName strings the same length (-> nicer output)
    int maxNameLen = 35;
    for(int i = 0; i < nTimers; i++)
      maxNameLen = std::max(maxNameLen, (int)fcnName.at(i).length());
    maxNameLen = maxNameLen + 5;
    for(int i = 0; i < nTimers; i++)
      fcnName.at(i).resize(maxNameLen,' ');
    // print result on proc 0
    std::string function = "function";
    function.resize(maxNameLen, ' ');
    PHIST_SOUT(PHIST_INFO, "======================================== TIMING RESULTS ============================================\n");
#ifdef PHIST_TIMINGS_FULL_TRACE
    PHIST_SOUT(PHIST_INFO, "%s  %10s  %10s  %10s  %10s  %10s  %10s\n", function.c_str(), "mtot.time", "mself.time", "count", "max.time", "avg.time", "min.time");
#else
    PHIST_SOUT(PHIST_INFO, "%s  %10s  %10s  %10s  %10s  %10s\n", function.c_str(), "mtot.time", "count", "max.time", "avg.time", "min.time");
#endif
    int nprocs = 1;
#ifdef PHIST_HAVE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif
    for(int i_ = 0; i_ < nTimers; i_++)
    {
      int i = sortedIndex.at(i_);
      if( i_ % 10 == 0 )
      {
        PHIST_SOUT(PHIST_INFO, "----------------------------------------------------------------------------------------------------\n");
      }
#ifdef PHIST_TIMINGS_FULL_TRACE
      PHIST_SOUT(PHIST_INFO, "%s  %10.3e  %10.3e  %10lu  %10.3e  %10.3e  %10.3e\n", fcnName.at(i).c_str(), maxTotalTime.at(i), maxSelfTime.at(i), numberOfCalls.at(i), maxTime.at(i), sumTotalTime.at(i)/numberOfCalls.at(i)/nprocs, minTime.at(i));
#else
      PHIST_SOUT(PHIST_INFO, "%s  %10.3e  %10lu  %10.3e  %10.3e  %10.3e\n", fcnName.at(i).c_str(), maxTotalTime.at(i), numberOfCalls.at(i), maxTime.at(i), sumTotalTime.at(i)/numberOfCalls.at(i)/nprocs, minTime.at(i));
#endif
    }
    PHIST_SOUT(PHIST_INFO, "====================================================================================================\n");
  }
}
