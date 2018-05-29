#!/bin/bash
export MPIEXEC=orterun
export MPIFLAGS="-np 1 --bind-to none"

# with GHOST
#export ENV="env OMP_SCHEDULE=static"

#with Tpetra+CUDA
#export ENV="env OMP_PROC_BIND=false env CUDA_MANAGED_FORCE_DEVICE_ALLOC=1 env CUDA_LAUNCH_BLOCKING=1"

#with Tpetra (CPU)
#export MPIFLAGS="-np 4 --bind-to numa"
#export ENV="env OMP_PROC_BIND=false"


