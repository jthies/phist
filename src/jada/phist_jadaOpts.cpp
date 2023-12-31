/*******************************************************************************************/
/* This file is part of the PHIST software available at https://bitbucket.org/essex/phist/ */
/* You may redistribute it and/or modify it under the terms of the BSD-style licence       */
/* included in this software.                                                              */
/*                                                                                         */
/* Contact: Jonas Thies (j.thies@tudelft.nl)                                               */
/*                                                                                         */
/*******************************************************************************************/
#include "phist_config.h"
#include "phist_tools.h"
#include "phist_jadaOpts.h"
#include <iostream>
#include <fstream>
#include <sstream>

#include <cstdio>
#include <cstring>

extern "C" void phist_jadaOpts_setDefaults(phist_jadaOpts *opts)
{
  opts->symmetry=phist_GENERAL;
  opts->numEigs=6; 
  opts->which=phist_LM; 
  opts->how=phist_STANDARD;

  opts->maxIters=300;
  opts->blockSize=1;
  opts->lookAhead=-1;
  opts->minBas=10;
  opts->maxBas=20;
  opts->convTol=1.0e-12;
  opts->relConvTol=0;

  opts->v0=NULL;
  opts->arno=1;
  opts->initialShift_r=0.0;
  opts->initialShift_i=0.0;
  opts->initialShiftIters=0;

  opts->innerSolvType=phist_GMRES;
  opts->innerSolvMaxBas=-1;
  opts->innerSolvMaxIters=20;
  opts->innerSolvBlockSize=-1;
  opts->innerSolvStopAfterFirstConverged=0;
  opts->innerSolvRobust=1;
  opts->innerSolvMaxProjectionSpace=-1;
  opts->innerSolvBaseTol=0.1;
  
  opts->preconOp=NULL;
  opts->preconType=phist_NO_PRECON;
  strcpy(opts->preconOpts,"none\0");
  opts->preconSkewProject=1;
  opts->preconUpdate=0;

  opts->customSolver=NULL;
  opts->customSolver_run1=NULL;
  opts->customSolver_run=NULL;

}

// helper function, finds and entry in a string and reads its value (which must sit next to it 
// separated only by whitespace
template<typename entryType>
void set_value(std::string key, entryType& entry, std::string file)
{
  size_t found = file.find(key);
  if (found!=std::string::npos)
  {
    std::istringstream iss(file.substr(found+key.length()));
    iss >> entry;
#if PHIST_OUTLEV>=PHIST_DEBUG
    std::ostringstream oss;
    oss << entry;
    PHIST_SOUT(PHIST_DEBUG, "Found entry %s=%s\n",key.c_str(),oss.str().c_str());
#endif
  }
  else
  {
    PHIST_SOUT(PHIST_DEBUG, "Did not find entry %s\n",key.c_str());
  }
}

extern "C" void phist_jadaOpts_fromFile(phist_jadaOpts* opts, const char* filename, int* iflag)
{
  PHIST_ENTER_FCN(__FUNCTION__);

  // this converts the input file to a string
  std::ifstream ifs;

  try 
  {
    ifs.open(filename);
  } catch (...) {*iflag=PHIST_CAUGHT_EXCEPTION; return;}
  
  std::string file((std::istreambuf_iterator<char>(ifs)),
                 std::istreambuf_iterator<char>());

  set_value("symmetry",opts->symmetry,file);
  set_value("numEigs",opts->numEigs,file);
  set_value("which",opts->which,file); 
  set_value("how",opts->how,file);

  set_value("maxIters",opts->maxIters,file);
  set_value("blockSize",opts->blockSize,file);
  set_value("lookAhead",opts->lookAhead,file);
  set_value("minBas",opts->minBas,file);
  set_value("maxBas",opts->maxBas,file);
  set_value("convTol",opts->convTol,file);
  set_value("relConvTol",opts->relConvTol,file);

  set_value("arno",opts->arno,file);
  set_value("initialShift_r",opts->initialShift_r,file);
  set_value("initialShift_i",opts->initialShift_i,file);
  set_value("initialShiftIters",opts->initialShiftIters,file);

  set_value("innerSolvType",opts->innerSolvType,file);
  set_value("innerSolvMaxBas",opts->innerSolvMaxBas,file);
  set_value("innerSolvMaxIters",opts->innerSolvMaxIters,file);
  set_value("innerSolvBlockSize",opts->innerSolvBlockSize,file);
  set_value("innerSolvBaseTol",opts->innerSolvBaseTol,file);
  
  set_value("preconType",opts->preconType,file);
  set_value("preconOpts",opts->preconOpts,file);
  set_value("preconSkewProject",opts->preconSkewProject,file);
  set_value("preconUpdate",opts->preconUpdate,file);

  set_value("innerSolvStopAfterFirstConverged",opts->innerSolvStopAfterFirstConverged,file);
  set_value("innerSolvRobust",opts->innerSolvRobust,file);
  set_value("innerSolvMaxProjectionSpace",opts->innerSolvMaxProjectionSpace,file);
}

// print jadaOpts to a file or stream. The result can be used as input for subsequent runs
extern "C" void phist_jadaOpts_toFile(phist_jadaOpts const *opts, FILE* stream)
{
  if (stream==NULL) return;
  fprintf(stream,"symmetry\t%s\n",matSym2str(opts->symmetry));
  fprintf(stream,"numEigs\t%d\n",opts->numEigs);
  fprintf(stream,"which\t%s\n",eigSort2str(opts->which)); 
  fprintf(stream,"how\t%s\n",eigExtr2str(opts->how));

  fprintf(stream,"maxIters\t%d\n",opts->maxIters);
  fprintf(stream,"blockSize\t%d\n",opts->blockSize);
  fprintf(stream,"lookAhead\t%d\n",opts->lookAhead);
  fprintf(stream,"minBas\t%d\n",opts->minBas);
  fprintf(stream,"maxBas\t%d\n",opts->maxBas);
  fprintf(stream,"convTol\t%4.2e\n",opts->convTol);
  fprintf(stream,"relConvTol\t%d\n",opts->relConvTol);

  fprintf(stream,"arno\t%d\n",opts->arno);
  fprintf(stream,"initialShift_r\t%e\n",opts->initialShift_r);
  fprintf(stream,"initialShift_i\t%e\n",opts->initialShift_i);
  fprintf(stream,"initialShiftIters\t%d\n",opts->initialShiftIters);

  fprintf(stream,"innerSolvType\t%s\n",linSolv2str(opts->innerSolvType));
  fprintf(stream,"innerSolvMaxBas\t%d\n",opts->innerSolvMaxBas);
  fprintf(stream,"innerSolvMaxIters\t%d\n",opts->innerSolvMaxIters);
  fprintf(stream,"innerSolvBlockSize\t%d\n",opts->innerSolvBlockSize);
  fprintf(stream,"innerSolvBaseTol\t%e\n",opts->innerSolvBaseTol);
  fprintf(stream,"preconType\t%s\n",precon2str(opts->preconType));
  fprintf(stream,"preconOpts\t%s\n",opts->preconOpts==NULL?"<NULL>":opts->preconOpts);
  fprintf(stream,"preconSkewProject\t%d\n",opts->preconSkewProject);
  fprintf(stream,"preconUpdate\t%d\n",opts->preconUpdate);
  fprintf(stream,"innerSolvStopAfterFirstConverged\t%d\n",opts->innerSolvStopAfterFirstConverged);
  fprintf(stream,"innerSolvRobust\t%d\n",opts->innerSolvRobust);
  fprintf(stream,"innerSolvMaxProjectionSpace\t%d\n",opts->innerSolvMaxProjectionSpace);
}

