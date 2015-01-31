#!/bin/bash
export OMP_NUM_THREADS=10
# ./phist_Dsubspacejada spinSZ26 1 20 SR 1.0e-8 300 1 30 50 1 10 0 0 1 1 &> primme_comparison_bs1.out
# ./phist_Dsubspacejada spinSZ26 1 20 SR 1.0e-8 300 4 30 50 4 10 0 0 1 1 &> primme_comparison_bs4.out

./phist_Dsubspacejada /elxfs/essex/hpc023/melven/spin26rcm.bin 1 20 SR 1.0e-8 300 4 28 60 4 8 0 0 1 1
