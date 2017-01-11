#include "include/craftConf.h"
#include "include/cpHelperFuncs.hpp"

int craftDebug            = CRAFT_DEBUG;
int craftDebugLog         = CRAFT_DEBUG_LOG;
int craftDebugProc        = CRAFT_DEBUG_PROC;
int craftReadCpOnRestart  = CRAFT_READ_CP_ON_RESTART;
std::string craftCpPath   = CRAFT_CP_PATH;
int craftEnabled          = CRAFT_ENABLE;
int craftUseSCR           = CRAFT_USE_SCR;
std::string craftCommSpawnPolicy    = CRAFT_COMM_SPAWN_POLICY;
std::string craftCommRecoveryPolicy = CRAFT_COMM_RECOVERY_POLICY;
std::string pbsNodeFile    = PBS_NODEFILE;

void getEnvParam(){
  char* envIn = new char[256];
 
  envIn = getenv ("CRAFT_DEBUG");
  if (envIn!=NULL){
    std::string envInStr;
    envInStr = envIn; 
    craftDebug= stringToNumber<int>(envInStr);
    craftDbg(4, "CRAFT_DEBUG: %d", craftDebug);
  }

  envIn = getenv ("CRAFT_DEBUG_PROC");
  if (envIn!=NULL){
    std::string envInStr;
    envInStr = envIn; 
    craftDebugProc = stringToNumber<int>(envInStr);
    craftDbg(4, "CRAFT_DEBUG_PROC: %d", craftDebugProc);
  }

  envIn = getenv ("CRAFT_ENABLE");
  if (envIn!=NULL){
    std::string envInStr;
    envInStr = envIn; 
    craftEnabled = stringToNumber<int>(envInStr);
    craftDbg(4, "CRAFT_ENABLE: %d", craftEnabled);
  }

  envIn = getenv ("CRAFT_READ_CP_ON_RESTART");
  if (envIn!=NULL){
    std::string envInStr;
    envInStr = envIn; 
    craftReadCpOnRestart = stringToNumber<int>(envInStr);
    craftDbg(4, "CRAFT_READ_CP_ON_RESTART: %d", craftReadCpOnRestart);
  }

  envIn = getenv ("CRAFT_CP_PATH");
  if (envIn!=NULL){
    craftCpPath = envIn; 
    craftDbg(4, "CRAFT_CP_PATH: %s", (craftCpPath).c_str());
  }

  envIn = getenv ("CRAFT_USE_SCR");
  if (envIn!=NULL){
    std::string envInStr;
    envInStr = envIn; 
    craftUseSCR = stringToNumber<int>(envInStr);
    craftDbg(4, "craftUseSCR : %d", craftUseSCR);
    if(craftUseSCR){
#ifndef SCR
      craftDbg(1, "craft is not compiled with SCR, thus env. variable CRAFT_USE_SCR can not be enabled(1)");
      craftUseSCR = 0;
#endif
    }
  }
  envIn = getenv ("CRAFT_COMM_RECOVERY_POLICY");
  if (envIn!=NULL){
    craftCommRecoveryPolicy = envIn; 
    craftDbg(4, "CRAFT_COMM_RECOVERY_POLICY: %s", (craftCommRecoveryPolicy).c_str());
  }
  envIn = getenv ("CRAFT_COMM_SPAWN_POLICY");
  if (envIn!=NULL){
    craftCommSpawnPolicy= envIn; 
    craftDbg(4, "CRAFT_COMM_SPAWN_POLICY: %s", (craftCommSpawnPolicy).c_str());
  }
  envIn = getenv ("PBS_NODEFILE");
  if (envIn!=NULL){
    pbsNodeFile = envIn; 
    craftDbg(4, "pbsNodeFile: %s", (pbsNodeFile).c_str());
  }
}


