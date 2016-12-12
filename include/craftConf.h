#ifndef __CRAFTCONF__
#define __CRAFTCONF__

#include <stdio.h>
#include <stdlib.h>
#include <string>

// 1 = enable , 2 = disable

#ifndef CRAFT_CP_READ_ON_RESTART 
#define CRAFT_CP_READ_ON_RESTART (1)
#endif

#ifndef CRAFT_CP_PATH
#define CRAFT_CP_PATH "./"
#endif

#ifndef CRAFT_ENABLE
#define CRAFT_ENABLE (1)
#endif

#ifndef CRAFT_USE_SCR

#ifdef SCR         // if SCR is define, then by default its status is ON (1)
#define CRAFT_USE_SCR (1)
#else
#define CRAFT_USE_SCR (0)
#endif

#endif


void getEnvParam();

extern int craftEnabled;
extern int craftCpReadOnRestart;
extern std::string craftCpPath;
extern int craftUseSCR;

#endif
