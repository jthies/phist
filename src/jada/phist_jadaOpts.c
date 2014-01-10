#include "phist_jadaOpts.h"
#include <stdlib.h>

void phist_jadaOpts_setDefaults(phist_jadaOpts_t *opts)
{
opts->numEigs=6; 
opts->which=LM; 

opts->maxIters=300; 
opts->blockSize=1; 
opts->minBas=10;
opts->maxBas=20;

opts->v0=NULL; 
opts->arno=1;
opts->initialShift=0.0; 

}
