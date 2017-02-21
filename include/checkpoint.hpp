#ifndef __CHECKPOINT_HPP__
#define __CHECKPOINT_HPP__

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <complex>
#include <sys/stat.h>
#include <utility>

#include "cpOptions.h"
#include "dataType.h"
#include "cpHelperFuncs.hpp"
#include "craftConf.h"
// ====== CP TYPES HEADERS ===== //  
#include "cpEnum.h"
#include "cpBase.hpp"

#ifdef SCR
extern "C"{
	#include <scr.h>
}
#endif


#define SSTR( x ) static_cast< std::ostringstream & >( \
								        ( std::ostringstream() << std::dec << x ) ).str()





class Checkpoint
{
protected:
	MPI_Comm cpMpiComm;
	std::string cpPath;
	std::string cpBasePath;	
	std::string cpPathVersion;
	std::string cpName;
	std::string cpVersionPrefix;
	bool cpUseSCR;
	bool cpCommitted;
	size_t cpVersion; 	
	size_t numBufferCps;
  bool restartStatus;

	int writeCpMetaData(const std::string filename);
	int mkCpDir(std::string path);
	int readCpMetaData();
  int SCRreadCpMetaData();
	int deleteBackupCp();
	int PFSwrite();
	int SCRwrite();
	int PFSread();
	int SCRread();
  int getCpVersion(); 
 
  typedef std::map<const std::string,CpBase *> cp_const_map;
  cp_const_map objects;
  
public:
	Checkpoint(const std::string name_, const MPI_Comm cpMpiComm);
	~Checkpoint();  
	void disableSCR();
	void commit();
  bool needRestart();
  int addToMap(std::string label, CpBase * p);
  MPI_Comm getCpComm();

  int read();
 	int write();
	int update();
  
  template <typename ...Params>
  void add(Params&&... params)
  {
    addCpType(this, std::forward<Params>(params)...);
  }

};

#include "addedCpTypes.hpp"

// implementation of addToMap() for anything that is CpBase,
// anything else will give an error message
int Checkpoint::addToMap(std::string label, CpBase * p)
{
  if( craftEnabled == false){
		craftDbg(0, "Checkpoint::addToMap(): CRAFT_ENABLE: %d", craftEnabled);
    return EXIT_FAILURE;
  }
  if (cpCommitted == true ){
			craftDbg(0, "This checkpoint is already committed. No data can be added to checkpoint after commit() call of a checkpoint");
    return EXIT_FAILURE;
  }
  objects[label] = p;
}




#endif
