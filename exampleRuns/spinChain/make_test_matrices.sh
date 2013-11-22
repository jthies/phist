#!/bin/bash
./SpinChainCRS 1 1 1 0 0 20 spin1.crs
./SpinChainCRS 1 1 0 0 0 20 spin2.crs
./SpinChainCRS 0 0 1 0 0 20 spin3.crs
./SpinChainCRS 0 0 0 1 0 20 spin4.crs
./SpinChainCRS 1 1 1 1 1 20 spin5.crs

binCRS2mm.py spin1.crs > spin1.mm
binCRS2mm.py spin2.crs > spin2.mm
binCRS2mm.py spin3.crs > spin3.mm
binCRS2mm.py spin4.crs > spin4.mm
binCRS2mm.py spin5.crs > spin5.mm

