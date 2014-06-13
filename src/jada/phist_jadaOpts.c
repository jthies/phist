#include "phist_config.h"
/* needs to be included before system headers for some intel compilers+mpi */
#ifdef PHIST_HAVE_MPI
#include <mpi.h>
#endif
#include "phist_jadaOpts.h"
#include <stdlib.h>

void phist_jadaOpts_setDefaults(phist_jadaOpts_t *opts)
{
opts->numEigs=6; 
opts->which=LM; 

opts->innerSolvType=GMRES;

opts->maxIters=300; 
opts->blockSize=1; 
opts->minBas=10;
opts->maxBas=20;
opts->convTol=1.0e-12;

opts->v0=NULL; 
opts->arno=1;
opts->initialShift=0.0; 

}
