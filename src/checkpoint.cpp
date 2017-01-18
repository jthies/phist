#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <cstdlib>

#include "include/checkpoint.hpp"

Checkpoint::Checkpoint( const std::string name_,
                        const MPI_Comm cpMpiComm_=MPI_COMM_WORLD
                       )
{
#ifndef AFT
  getEnvParam();
#endif  
  static std::set<std::string> nameKeys;    // to determine the uniqueness of different checkpoint names.
  std::pair<std::set<std::string>::iterator,bool> ret;
  ret         =  nameKeys.insert( name_);
  if(ret.second==false){
		craftDbg(3, "checkpoint with name '%s' already exists.", name_.c_str());
  } 
	cpName 				= name_;
  checkDirectoryName(&cpName);                        // checks if the given name makes a valid directory name.
	cpMpiComm 	= cpMpiComm_;	
	cpBasePath 	= craftCpPath;
	cpPath 			= cpBasePath + "/" + "Checkpoint-" + cpName;
	craftDbg(3, "cpPath: %s", cpPath.c_str());
	cpVersionPrefix = "v-";
	cpCommitted = false;
	cpVersion 	= 0;
	numBufferCps= 2;
  restartStatus = false; 

  cpUseSCR    = craftUseSCR ;

	if(cpMpiComm == MPI_COMM_WORLD) {
		craftDbg(3, "cpMpiComm = MPI_COMM_WORLD");
  }
}	

Checkpoint::~Checkpoint(){
#ifdef SCR
	if(cpUseSCR){
    static bool firstRun=true;
    if(firstRun){
		  craftDbg(3, "====== SCR_Finalize() ====== ");
		  SCR_Finalize();
      firstRun = false;
    }
	}
#endif
}

void Checkpoint::disableSCR(){
#ifdef SCR
  craftDbg(1, "SCR is disabled for '%s' checkpoint. Use SCR to reduce the impact of Checkpoint/restart", cpName.c_str());
  cpUseSCR = false;
#endif
}

void Checkpoint::commit(){
  if (cpUseSCR){
#ifdef SCR
    craftDbg(3, "====== SCR_Init() cpUseSCR is true %d ====== ", cpUseSCR);
    SCR_Init(&cpMpiComm);     // SCR has to be initialized only once. It is important in case of multilevelCP
#endif
  }
  else{
    mkCpDir(cpPath);	
	  MPI_Barrier(cpMpiComm);			// makes sure that cp-directory is created before any other process tries to write in it.
  }
	cpCommitted = true;	
}

int Checkpoint::mkCpDir(std::string path){
	std::string cmd = "mkdir -p " + path;
	system (cmd.c_str());
  sync();
	return EXIT_SUCCESS;
}

int Checkpoint::writeCpMetaData(const std::string filename){			// writes a tmp file first and then replaces it after complete write
  std::string filenameTmp;
	filenameTmp = filename + ".tmp"; 
	std::ofstream fstr;
	fstr.open ((filenameTmp).c_str(), std::ios::out );	
	if(fstr.is_open()){
		fstr << cpVersion << std::endl;
		fstr.close();
	}
   else{
		std::cerr << "Can't open file " << filenameTmp << std::endl;			
		return EXIT_FAILURE;
	 }
	std::string cmd = "mv " + filenameTmp + " " + filename;
	system (cmd.c_str());
  return EXIT_SUCCESS;
}

int Checkpoint::readCpMetaData(){
	std::string filename, line;
	filename = cpPath + "/" + "metadata.ckpt"; 
	std::ifstream fstr;
	fstr.open ((filename).c_str(), std::ios::in );	
	craftDbg(3, "Metadata filename: %s", (filename).c_str());
	if(fstr.is_open()){
		while ( getline (fstr,line) )
		{
			cpVersion = atoi(line.c_str());
		}
		fstr.close();
		craftDbg(1, "cpVersion to read: %d", cpVersion);
    return EXIT_SUCCESS;
	}
  else{
	  std::cerr << "Can't open file " << filename << std::endl;			
	  return EXIT_FAILURE;
  }
}

int Checkpoint::SCRreadCpMetaData(){
#ifdef SCR
  // === writeMetaData file === // 
	int myrank_ = -1;	
	MPI_Comm_rank(cpMpiComm, &myrank_);
	std::string myrank = numberToString(myrank_);
	std::string filename = cpPath + "/" + "metadata" + myrank + ".ckpt"; 
  char * fNameScrPath = new char[256]; 
	char * tmpFilename  = new char[256];
	strcpy(tmpFilename, filename.c_str()); 
	SCR_Route_file(tmpFilename, fNameScrPath);
	filename = fNameScrPath;
	std::ifstream fstr;
	fstr.open ((filename).c_str(), std::ios::in );	
	craftDbg(3, "Metadata filename: %s", (filename).c_str());
  std::string line;
	if(fstr.is_open()){
		while ( getline (fstr,line) )
		{
			cpVersion = atoi(line.c_str());
		}
		fstr.close();
		craftDbg(1, "cpVersion to read: %d", cpVersion);
	  return EXIT_SUCCESS;
	}
  else{
	  std::cerr << "Can't open file " << filename << std::endl;			
	  return EXIT_FAILURE;
  }
#endif
}

int Checkpoint::deleteBackupCp(){
	int myrank = -1;	
	MPI_Comm_rank(cpMpiComm, &myrank);
	if(myrank == 0){
		std::string toRmDir = cpPath + "/" + cpVersionPrefix + SSTR(cpVersion-numBufferCps-1);
		std::string cmd = "rm -r " + toRmDir;	
		struct stat sb;
		if (stat(toRmDir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
		{
			system ( cmd.c_str());
 		}
	}
	return EXIT_SUCCESS;
}


int Checkpoint::update(){
  if( craftEnabled == 0){
		craftDbg(3, "Checkpoint::update(): CRAFT_ENABLE: %d", craftEnabled);
    return EXIT_FAILURE;
  }
  craftDbg(1, "updating Checkpoint...");
	for(cp_const_map:: iterator it = objects.begin(); it != objects.end(); ++it)
	{
		it->second->update();
	}	
	return EXIT_SUCCESS;
}

int Checkpoint::write()				// TODO: make two version of write. 1) PFS 2) SCR. and do MPI/IO for PFS 
{	
  if( !craftEnabled ){
		craftDbg(3, "Checkpoint::write(): CRAFT_ENABLE: %d", craftEnabled);
    return EXIT_FAILURE;
  }				
	craftDbg(1, "writing Checkpoint...");
  int ret=EXIT_FAILURE;
  if( cpUseSCR ) {
#ifdef SCR
		ret = SCRwrite();
#endif
  }
  else{
		ret = PFSwrite();
  }
	return EXIT_SUCCESS;
}

int Checkpoint::SCRwrite(){
#ifdef SCR
	SCR_Start_checkpoint();
	int myrank_ = -1;
	MPI_Comm_rank(cpMpiComm, &myrank_);
	std::string myrank = numberToString(myrank_);
	std::string filename;
	cpPathVersion = cpPath + "/" + cpVersionPrefix + SSTR(cpVersion);

  for (cp_const_map::iterator it=objects.begin(); it!=objects.end(); it++)
 	{
		filename = cpPathVersion + "/" + (it->first) + myrank + ".ckpt"; 
		char * fNameScrPath = new char[256]; 
		char * tmpFilename  = new char[256];
	  strcpy(tmpFilename, filename.c_str()); 
		SCR_Route_file(tmpFilename, fNameScrPath);
		filename = fNameScrPath;
		craftDbg(2, "SCRwrite: writing file: %s", filename.c_str());
    if(EXIT_SUCCESS != it->second->write(&filename)){
      craftErr("SCRwrite not successful @ %s:%d",
              __FILE__, __LINE__
      );
      return EXIT_FAILURE;
    }
	}
  // === writeMetaData file === // 
	filename = cpPath + "/" + "metadata" + myrank + ".ckpt"; 
  char * fNameScrPath = new char[256]; 
	char * tmpFilename  = new char[256];
	strcpy(tmpFilename, filename.c_str()); 
	SCR_Route_file(tmpFilename, fNameScrPath);
	filename = fNameScrPath;
	craftDbg(2, "SCRwrite: writing file: %s", filename.c_str());
	writeCpMetaData(filename);
	int valid = 1;
	SCR_Complete_checkpoint(valid);     // valid flag should be 1 if every CP is successfully written. TODO: check the flag from each write call.
	++cpVersion;
#endif
	return EXIT_SUCCESS;
}


int Checkpoint::PFSwrite(){
	int myrank_ = -1;
	MPI_Comm_rank(cpMpiComm, &myrank_);
	std::string myrank = numberToString(myrank_);
	std::string filename;
	cpPathVersion = cpPath + "/" + cpVersionPrefix + SSTR(cpVersion);
	mkCpDir(cpPathVersion);	
	MPI_Barrier(cpMpiComm);			// makes sure that cpVersion directory is created before any process tries to write in it.
  for (cp_const_map::iterator it=objects.begin(); it!=objects.end(); it++)
 	{
#ifdef MPIIO																										// Parallel PFS IO
		filename = cpPathVersion + "/" + (it->first) + ".ckpt"; 
#else																														// Serial PFS IO
		filename = cpPathVersion + "/" + (it->first) + myrank + ".ckpt"; 
#endif
    craftDbg(2, "PFSwrite(): writing file: %s", filename.c_str());
		if(EXIT_SUCCESS != it->second->write(&filename)){
        craftErr("PFSwrite not successful @ %s:%d",
                __FILE__, __LINE__
        );
				return EXIT_FAILURE;
    }
	}
	MPI_Barrier(cpMpiComm);			// TODO: do MPI_Gather here 
  // === writeMetaData file === // 
	if(myrank_ == 0){
		filename = cpPath + "/" + "metadata.ckpt"; 
	  writeCpMetaData(filename);
	}
	++cpVersion;
	deleteBackupCp();
	return EXIT_SUCCESS;
}


int Checkpoint::read()				// TODO: make two version of read. 1) PFS 2) SCR
{
  if( craftEnabled == 0){
		craftDbg(3, "Checkpoint::read(): CRAFT_ENABLE: %d", craftEnabled);
    return EXIT_FAILURE;
  }
	craftDbg(1, "reading Checkpoint...");
  int ret=EXIT_FAILURE;
  if( cpUseSCR ) {
#ifdef SCR
		ret = SCRread();
#endif
  }
  else{
		ret = PFSread();
  }
  return ret;
}

int Checkpoint::PFSread(){
  int ret=EXIT_FAILURE;
	ret = readCpMetaData();									// cpVersion is updated here
  if (ret == EXIT_FAILURE){  return EXIT_FAILURE; }
	int myrank_ = -1;
	MPI_Comm_rank(cpMpiComm, &myrank_);
	std::string myrank = numberToString(myrank_);
	std::string filename;
	cpPathVersion = cpPath + "/" + cpVersionPrefix + SSTR(cpVersion); 
  for (cp_const_map::iterator it=objects.begin(); it!=objects.end(); it++)
 	{			
#ifdef MPIIO																										// Parallel PFS IO 
		filename = cpPathVersion + "/" + (it->first) + ".ckpt";				
#else																														// Serial PFS IO 
		filename = cpPathVersion + "/" + (it->first) + myrank + ".ckpt"; 
#endif
    craftDbg(2, "PFSread(): reading file: %s", filename.c_str());
		if(EXIT_SUCCESS != it->second->read(&filename))	{
      craftErr("PFSread not successful @ %s:%d",
              __FILE__, __LINE__
      );
			return EXIT_FAILURE;
    }
  }
	MPI_Barrier(cpMpiComm);
	++cpVersion;
	return EXIT_SUCCESS;
}

int Checkpoint::SCRread(){
#ifdef SCR
	int myrank_ = -1;
	MPI_Comm_rank(cpMpiComm, &myrank_);
	std::string myrank = numberToString(myrank_);
	std::string filename;
	cpPathVersion = cpPath + "/" + cpVersionPrefix + SSTR(cpVersion); 
  for (cp_const_map::iterator it=objects.begin(); it!=objects.end(); it++)
 	{			
		filename = cpPathVersion + "/" + (it->first) + myrank + ".ckpt"; 
		char * fNameScrPath = new char[256]; 
		char * tmpFilename = new char[256];
		strcpy(tmpFilename, filename.c_str()); 
		SCR_Route_file(tmpFilename, fNameScrPath);
		filename = fNameScrPath;
		craftDbg(2, "SCRread: reading file %s", filename.c_str());
		if(EXIT_SUCCESS != it->second->read(&filename))	{
      craftErr("SCRread not successful @ %s:%d",
              __FILE__, __LINE__
      );
			return EXIT_FAILURE;
    }
  }

	MPI_Barrier(cpMpiComm);
	++cpVersion;
#endif

	return EXIT_SUCCESS;
}

bool Checkpoint::needRestart(){
  if( craftReadCpOnRestart == 0 ){    // in case, that checkpoints exist but user wants to restart from the beginning. This flag is set by ENV variables.
    craftDbg(3, "Environment variable CRAFT_READ_CP_ON_RESTART is disabled(0) => Not restarting from last checkpoint");
    return false;   
  }

  int ret=-1;
  if(cpVersion!=0) { 
    craftDbg(3, "%s::cpVersion != 0 cpVersion= %d", cpName.c_str(), cpVersion);
    return false;     // cpVersion is set to 0 at initialization. This means that if the inner L2 loop rotates. The checkpoint will not be read further.
  }

  if(cpVersion==0)   // means either its the 1st-run or restarted-run. Now we differentiate between them
  {     //  In case of 1st-run, there should be no CP, else there should be CP

// === SCR case === // 
#ifdef SCR
    if(cpUseSCR){
        ret = SCRreadCpMetaData();
        if(ret == EXIT_FAILURE){  
          craftDbg(3, "No SCRcpMetaData file found, thus its 1st-run.");
          return false;
        }
        else{    
          craftDbg(3, "SCRcpMetaData file is found,thus Checkpoint exists.");
          return true;        // SCR only has one level of CP. So if it is requesting a restart, app needs restart
        }
    }
#endif

// === PFS case === // 
    ret = readCpMetaData();
    if(ret == EXIT_FAILURE){  
      craftDbg(3, "No PFScpMetaData file found, thus its 1st-run.");
      return false;
    }
    else {
      craftDbg(3, "PFScpMetaData file is found,thus Checkpoint exists.");
      return true;
    }
  }
}

