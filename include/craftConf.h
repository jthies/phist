#ifndef __CRAFTCONF_H__
#define __CRAFTCONF_H__

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

// 1 = enable , 0 = disable
// 1 = true   , 0 = false


#ifndef CRAFT_ENABLE
#define CRAFT_ENABLE (1)
#endif

#ifndef CRAFT_DEBUG
#define CRAFT_DEBUG (0)
#endif

#ifndef CRAFT_DEBUG_LOG
#define CRAFT_DEBUG_LOG (1)
#endif


#ifndef CRAFT_DEBUG_PROC
#define CRAFT_DEBUG_PROC (0)
#endif

#ifndef CRAFT_READ_CP_ON_RESTART 
#define CRAFT_READ_CP_ON_RESTART (1)
#endif

#ifndef CRAFT_CP_PATH
#define CRAFT_CP_PATH "./"
#endif

#ifndef CRAFT_USE_SCR
#ifdef SCR         // if SCR is define, then by default its status is ON (1)
#define CRAFT_USE_SCR (1)
#else
#define CRAFT_USE_SCR (0)
#endif
#endif

#ifndef CRAFT_COMM_RECOVERY_POLICY 
#define CRAFT_COMM_RECOVERY_POLICY "NON-SHRINKING"
#endif

#ifndef CRAFT_COMM_SPAWN_POLICY
#define CRAFT_COMM_SPAWN_POLICY   "NO-REUSE"
#endif

#ifndef PBS_NODEFILE
#define PBS_NODEFILE "PBS_NODEFILE"
#endif

void getEnvParam();

extern int craftDebug;
extern int craftDebugLog;
extern int craftDebugProc;
extern int craftEnabled;
extern int craftReadCpOnRestart;
extern std::string craftCpPath;
extern int craftUseSCR;
extern std::string craftCommRecoveryPolicy;
extern std::string craftCommSpawnPolicy;
extern std::string pbsNodeFile;

extern char * craftLogFile;

#endif
