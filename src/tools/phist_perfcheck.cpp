/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (Jonas.Thies@DLR.de)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif

#include "phist_perfcheck.hpp"
#include "phist_tools.h"

#include <algorithm>
#include <sstream>
#include <cmath>
#include <cstdarg>
#include <cstdio>

#ifdef PHIST_PERFCHECK
# ifndef PHIST_USE_TEUCHOS_TIMEMONITOR
namespace phist_PerfCheck
{
  Timer::TimeDataMap PerfCheckTimer::expectedResults_;

  // functor for sorting by keys in another vector
  struct SortClass
  {
    SortClass(const std::vector<double>& keys) : _keys(keys) {}
    bool operator() (int i, int j) {return _keys[i] > _keys[j];}
    const std::vector<double>& _keys;
  };


    double PerfCheckTimer::total_GByte_transferred_in_kernels=0.0;
    double PerfCheckTimer::total_Gflops_performed_in_kernels=0.0;
    double PerfCheckTimer::total_time_spent_in_kernels=0.0;

  // calculates the results and prints them
  void PerfCheckTimer::summarize(int verbosity)
  {
    if (PHIST_OUTLEV<verbosity) return;
    int ierr = 0;
    int nTimers = timingResults_.size();
    int nExpectedTimers = expectedResults_.size();

#ifdef PHIST_HAVE_MPI
    MPI_Comm mpi_comm=phist_get_default_comm();
    // consider only timers from the first process!
    PHIST_CHK_MPIERR(ierr = MPI_Bcast(&nTimers, 1, MPI_INT, 0, mpi_comm), ierr);
    PHIST_CHK_MPIERR(ierr = MPI_Bcast(&nExpectedTimers, 1, MPI_INT, 0, mpi_comm), ierr);
#endif

    PHIST_CHK_IERR(ierr = (nTimers != nExpectedTimers) ? -1 : 0, ierr);

    if( nTimers == 0 ) return;
    
    int rank=0,nproc=1;
#ifdef PHIST_HAVE_MPI
    MPI_Comm_rank(mpi_comm,&rank);
    MPI_Comm_size(mpi_comm,&nproc);
#endif

    FILE* ofile= (rank==0)? stdout: NULL;
        
#ifdef PHIST_PERFCHECK_SEPARATE_OUTPUT
    int ndigits=std::ceil(std::log(nproc)/std::log(10));
    char fname_fmt[32], fname[32];

    snprintf(fname_fmt,32,"phist-perfcheck-P%%%dd.txt",ndigits);
    snprintf(fname,32,fname_fmt,rank);
    PHIST_SOUT(PHIST_INFO,"perfcheck results are written to searate files of the form '%s'\n",fname);
    ofile=fopen(fname,"w+");
#endif
    // get timer names from the first process
    std::string fcnNameList;
    for(TimeDataMap::const_iterator it = timingResults_.begin(); it != timingResults_.end(); it++)
      fcnNameList.append( it->first + '\n' );
    int strLen = fcnNameList.length();
#ifdef PHIST_HAVE_MPI
    PHIST_CHK_MPIERR(ierr = MPI_Bcast(&strLen, 1, MPI_INT, 0, mpi_comm), ierr);
#endif
    char* strBuf = new char[strLen];
    fcnNameList.copy(strBuf, strLen);
    strBuf[strLen-1] = '\0';
#ifdef PHIST_HAVE_MPI
    PHIST_CHK_MPIERR(ierr = MPI_Bcast(strBuf, strLen, MPI_CHAR, 0, mpi_comm), ierr);
#endif
    std::istringstream iss(strBuf);

    std::vector<std::string> fcnName(nTimers);
    std::vector<unsigned long> numberOfCalls(nTimers);
    std::vector<double> minTime(nTimers), minExpected(nTimers);
    std::vector<double> maxTime(nTimers), maxExpected(nTimers);
    std::vector<double> maxTotalTime(nTimers), maxTotalExpected(nTimers);
    std::vector<double> sumTotalTime(nTimers), sumTotalExpected(nTimers);

    std::vector<double> tmp(nTimers);

    // convert / copy everything into vectors
    for(int i = 0; i < nTimers; i++)
    {
      std::getline(iss, fcnName.at(i));
      numberOfCalls.at(i) = timingResults_[fcnName.at(i)].numberOfCalls;
      minTime.at(i) = timingResults_[fcnName.at(i)].minTime;
      minExpected.at(i) = expectedResults_[fcnName.at(i)].minTime;
      maxTime.at(i) = timingResults_[fcnName.at(i)].maxTime;
      maxExpected.at(i) = expectedResults_[fcnName.at(i)].maxTime;
      maxTotalTime.at(i) = timingResults_[fcnName.at(i)].totalTime;
      maxTotalExpected.at(i) = expectedResults_[fcnName.at(i)].totalTime;
    }
    delete[] strBuf;

#ifdef PHIST_HAVE_MPI
    // reductions over mpi processes
    tmp = minTime;
    PHIST_CHK_MPIERR(ierr = MPI_Reduce(&tmp[0],&minTime[0],nTimers,MPI_DOUBLE,MPI_MIN,0,mpi_comm), ierr);
    tmp = minExpected;
    PHIST_CHK_MPIERR(ierr = MPI_Reduce(&tmp[0],&minExpected[0],nTimers,MPI_DOUBLE,MPI_MIN,0,mpi_comm), ierr);

    tmp = maxTime;
    PHIST_CHK_MPIERR(ierr = MPI_Reduce(&tmp[0],&maxTime[0],nTimers,MPI_DOUBLE,MPI_MAX,0,mpi_comm), ierr);
    tmp = maxExpected;
    PHIST_CHK_MPIERR(ierr = MPI_Reduce(&tmp[0],&maxExpected[0],nTimers,MPI_DOUBLE,MPI_MAX,0,mpi_comm), ierr);

    tmp = maxTotalTime;
    PHIST_CHK_MPIERR(ierr = MPI_Reduce(&tmp[0],&sumTotalTime[0],nTimers,MPI_DOUBLE,MPI_SUM,0,mpi_comm), ierr);
    PHIST_CHK_MPIERR(ierr = MPI_Reduce(&tmp[0],&maxTotalTime[0],nTimers,MPI_DOUBLE,MPI_MAX,0,mpi_comm), ierr);
    tmp = maxTotalExpected;
    PHIST_CHK_MPIERR(ierr = MPI_Reduce(&tmp[0],&sumTotalExpected[0],nTimers,MPI_DOUBLE,MPI_SUM,0,mpi_comm), ierr);
    PHIST_CHK_MPIERR(ierr = MPI_Reduce(&tmp[0],&maxTotalExpected[0],nTimers,MPI_DOUBLE,MPI_MAX,0,mpi_comm), ierr);
#endif

    std::vector<double> sortBy(nTimers);
    for(int i = 0; i < nTimers; i++)
    {
      // sort by difference to expectation
      sortBy[i] = std::abs(maxTotalTime[i]-maxTotalExpected[i]);
      // sort by total time consumed by kernel
      //sortBy[i]=maxTotalTime[i];
    }
    std::vector<int> sortedIndex(nTimers);
    for(int i = 0; i < nTimers; i++)
    {
      sortedIndex.at(i) = i;
    }
    SortClass sortFunc(sortBy);
    std::sort(sortedIndex.begin(), sortedIndex.end(), sortFunc);

    // make all fcnName strings the same length (-> nicer output) (and split fcnName and formula)
    int maxNameLen = 55;
    std::vector<std::string> fcnFormula(nTimers);
    for(int i = 0; i < nTimers; i++)
    {
      int j = fcnName.at(i).find(' ');
      fcnFormula.at(i) = fcnName.at(i).substr(j+1,std::string::npos);
      fcnName.at(i) = fcnName.at(i).substr(0,j);
      maxNameLen = std::max(maxNameLen, (int)fcnName.at(i).length());
      maxNameLen = std::max(maxNameLen, (int)fcnFormula.at(i).length());
    }
    maxNameLen = maxNameLen + 5;
    for(int i = 0; i < nTimers; i++)
    {
      fcnName.at(i).resize(maxNameLen,' ');
    }

    double global_GBytes=total_GByte_transferred_in_kernels;
    double global_Gflops=total_Gflops_performed_in_kernels;
    double global_time=total_time_spent_in_kernels;

#ifdef PHIST_HAVE_MPI
# ifdef PHIST_PERFCHECK_SEPARATE_OUTPUT
    PHIST_CHK_MPIERR(ierr=MPI_Allreduce(&total_GByte_transferred_in_kernels,
                                            &global_GBytes, 1, MPI_DOUBLE,MPI_SUM,mpi_comm),ierr);

    PHIST_CHK_MPIERR(ierr=MPI_Allreduce(&total_Gflops_performed_in_kernels,
                                            &global_Gflops, 1, MPI_DOUBLE,MPI_SUM,mpi_comm),ierr);
# endif
    PHIST_CHK_MPIERR(ierr=MPI_Allreduce(&total_time_spent_in_kernels,
                                            &global_time, 1, MPI_DOUBLE,MPI_MAX,mpi_comm),ierr);
#endif
    
    // print result on proc 0 or on all procs to a file.
    // We could also implement e.g. that only the ranks on node 0 print the results
    // in this way
    if (ofile==NULL) return;
    
    std::string function = "function(dim) / (formula)";
    function.resize(maxNameLen, ' ');
    fprintf(ofile, "================================================== PERFORMANCE CHECK RESULTS =====================================================\n");
    fprintf(ofile, "%s  %10s  %10s  %10s  %10s  %10s\n", function.c_str(), "total time", "%roofline", "count", "max.%roofline", "min.%roofline");
    int nprocs = 1;
#ifdef PHIST_HAVE_MPI
    MPI_Comm_size(mpi_comm, &nprocs);
#endif
    double sumMaxTotalExpected = 0., sumMaxTotalTime = 0.;
    for(int i_ = 0; i_ < nTimers; i_++)
    {
      int i = sortedIndex.at(i_);
      if( i_ % 10 == 0 )
      {
        fprintf(ofile, "----------------------------------------------------------------------------------------------------------------------------------\n");
      }
      fprintf(ofile, "%s  %10.3e  %10.3g  %10lu  %10.3g  %10.3g\n", fcnName.at(i).c_str(), maxTotalTime.at(i), 100*maxTotalExpected.at(i)/maxTotalTime.at(i), numberOfCalls.at(i),
          100*minExpected.at(i)/minTime.at(i), 100*maxExpected.at(i)/maxTime.at(i));
      fprintf(ofile, " %s\n", fcnFormula.at(i).c_str());
      sumMaxTotalExpected += maxTotalExpected.at(i);
      sumMaxTotalTime += maxTotalTime.at(i);
    }
    fprintf(ofile, "==================================================================================================================================\n");
    std::string strTotal = "total";
    strTotal.resize(maxNameLen, ' ');
    fprintf(ofile, "%s  %10.3e  %10.3g\n", strTotal.c_str(), sumMaxTotalExpected, 100*sumMaxTotalExpected/sumMaxTotalTime);
    fprintf(ofile, "==================================================================================================================================\n");

    fprintf(ofile, "\nPERFORMANCE SUMMARY FOR PERFCHECKED KERNELS:\n");
    fprintf(ofile, "----------------------------------------------------------------------------------------------------------------------------------\n");
#ifdef PHIST_PERFCHECK_SEPARATE_OUTPUT
// note: we do not sum up all the memory transfers yet
//    fprintf(ofile, "estimated GB transferred to/from memory  on this process:\t%8.4e\n", 
//        total_GByte_transferred_in_kernels);
//    fprintf(ofile, "estimated bandwidth [GB/s] achieved      on this process:\t%8.4e\n", 
//        total_GByte_transferred_in_kernels / total_time_spent_in_kernels);
    fprintf(ofile, "estimated performance [Gflop/s] achieved on this process:\t%8.4e\n", total_Gflops_performed_in_kernels  / total_time_spent_in_kernels);
#endif

//    fprintf(ofile, "estimated GB transferred to/from memory  in total:\t%8.4e\n", global_GBytes);
//    fprintf(ofile, "estimated bandwidth [GB/s] achieved      in total:\t%8.4e\n", global_GBytes/global_time);
    fprintf(ofile, "Estimated performance achieved in measured kernels: %6.2e GFlop/s\n"
                   "Total time spent               in measured kernels: %6.2e seconds\n"
                   "(conservative estimate, some flops in fused kernels are not yet counted)\n",
        global_Gflops/global_time, global_time);
    fprintf(ofile, "----------------------------------------------------------------------------------------------------------------------------------\n");

    if (ofile!=stdout) fclose(ofile);
  }
}
# endif  
#endif
