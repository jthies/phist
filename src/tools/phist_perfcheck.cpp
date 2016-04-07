#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#else
#error "this file currently requires MPI"
#endif

#include "phist_perfcheck.hpp"
#include "phist_macros.h"

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


  // calculates the results and prints them
  void PerfCheckTimer::summarize(int verbosity)
  {
    if (PHIST_OUTLEV<verbosity) return;
    int ierr = 0;
    int nTimers = timingResults_.size();
    int nExpectedTimers = expectedResults_.size();

    // consider only timers from the first process!
    PHIST_CHK_IERR(ierr = MPI_Bcast(&nTimers, 1, MPI_INT, 0, MPI_COMM_WORLD), ierr);
    PHIST_CHK_IERR(ierr = MPI_Bcast(&nExpectedTimers, 1, MPI_INT, 0, MPI_COMM_WORLD), ierr);

    PHIST_CHK_IERR(ierr = (nTimers != nExpectedTimers) ? -1 : 0, ierr);

    if( nTimers == 0 ) return;
    
    int rank=0,nproc=1;
#ifdef PHIST_HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
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
    PHIST_CHK_IERR(ierr = MPI_Bcast(&strLen, 1, MPI_INT, 0, MPI_COMM_WORLD), ierr);
    char* strBuf = new char[strLen];
    fcnNameList.copy(strBuf, strLen);
    strBuf[strLen-1] = '\0';
    PHIST_CHK_IERR(ierr = MPI_Bcast(strBuf, strLen, MPI_CHAR, 0, MPI_COMM_WORLD), ierr);
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

    // reductions over mpi processes
    tmp = minTime;
    PHIST_CHK_IERR(ierr = MPI_Reduce(&tmp[0],&minTime[0],nTimers,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD), ierr);
    tmp = minExpected;
    PHIST_CHK_IERR(ierr = MPI_Reduce(&tmp[0],&minExpected[0],nTimers,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD), ierr);

    tmp = maxTime;
    PHIST_CHK_IERR(ierr = MPI_Reduce(&tmp[0],&maxTime[0],nTimers,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD), ierr);
    tmp = maxExpected;
    PHIST_CHK_IERR(ierr = MPI_Reduce(&tmp[0],&maxExpected[0],nTimers,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD), ierr);

    tmp = maxTotalTime;
    PHIST_CHK_IERR(ierr = MPI_Reduce(&tmp[0],&sumTotalTime[0],nTimers,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD), ierr);
    PHIST_CHK_IERR(ierr = MPI_Reduce(&tmp[0],&maxTotalTime[0],nTimers,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD), ierr);
    tmp = maxTotalExpected;
    PHIST_CHK_IERR(ierr = MPI_Reduce(&tmp[0],&sumTotalExpected[0],nTimers,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD), ierr);
    PHIST_CHK_IERR(ierr = MPI_Reduce(&tmp[0],&maxTotalExpected[0],nTimers,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD), ierr);

    // sort by difference to expectation
    std::vector<double> maxTotalDiff(nTimers);
    for(int i = 0; i < nTimers; i++)
      maxTotalDiff[i] = std::abs(maxTotalTime[i]-maxTotalExpected[i]);
    std::vector<int> sortedIndex(nTimers);
    for(int i = 0; i < nTimers; i++)
      sortedIndex.at(i) = i;
    SortClass sortFunc(maxTotalDiff);
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
    // print result on proc 0 or on all procs to a file.
    // We could also implement e.g. that only the ranks on node 0 print the results
    // in this way
    if (ofile==NULL) return;
    std::string function = "function(dim) / (formula)";
    function.resize(maxNameLen, ' ');
    fprintf(ofile, "================================================== PERFORMANCE CHECK RESULTS =====================================================\n");
    fprintf(ofile, "%s  %10s  %10s  %10s  %10s  %10s\n", function.c_str(), "mtot.exp", "%peak-perf", "count", "max.%peak", "min.%peak");
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    double sumMaxTotalExpected = 0., sumMaxTotalTime = 0.;
    for(int i_ = 0; i_ < nTimers; i_++)
    {
      int i = sortedIndex.at(i_);
      if( i_ % 10 == 0 )
      {
        fprintf(ofile, "----------------------------------------------------------------------------------------------------------------------------------\n");
      }
      fprintf(ofile, "%s  %10.3e  %10.3g  %10lu  %10.3g  %10.3g\n", fcnName.at(i).c_str(), maxTotalExpected.at(i), 100*maxTotalExpected.at(i)/maxTotalTime.at(i), numberOfCalls.at(i),
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
    if (ofile!=stdout) fclose(ofile);
  }
}
# endif  
#endif
