#!/bin/bash

# create a temporary build environment
spack env create ci-PhistEnv SpackEnv-gcc-9.2.0-openmpi-4.0.2.yaml
spack env activate ci-PhistEnv

# install all dependencies if not yet available
spack concretize
spack install

