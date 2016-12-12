#include "include/craftConf.h"
#include "include/cpHelperFuncs.hpp"

int craftCpReadOnRestart  = CRAFT_CP_READ_ON_RESTART;
std::string craftCpPath   = CRAFT_CP_PATH;
int craftEnabled          = CRAFT_ENABLE;
int craftUseSCR           = CRAFT_USE_SCR;

void getEnvParam(){
  char* envIn = new char[256];

  envIn = getenv ("CRAFT_ENABLE");
  if (envIn!=NULL){
    std::string envInStr;
    envInStr = envIn; 
    craftEnabled = stringToNumber<int>(envInStr);
    craftDbg(4, "CRAFT_ENABLE: %d", craftEnabled);
  }

  envIn = getenv ("CRAFT_CP_READ_ON_RESTART");
  if (envIn!=NULL){
    std::string envInStr;
    envInStr = envIn; 
    craftCpReadOnRestart = stringToNumber<int>(envInStr);
    craftDbg(4, "CRAFT_CP_READ_ON_RESTART: %d", craftCpReadOnRestart);
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
  
}


