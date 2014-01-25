#!/usr/bin/env python
import string
import math

# open the file
with open( 'spinChain_fortran.pbs.in' ) as infile:
  # read it
  intext = string.Template( infile.read() )

  for nn in [1, 2, 4, 8, 16, 32, 64]:
    for L in [22, 24, 26, 28, 30, 32, 34]:
      if math.log(nn,2) < (L-22)/2:
        continue
      with open( 'spinChain{0}_fortran_nn{1}.pbs'.format(L,nn), 'w+' ) as outfile:
        outfile.write( intext.safe_substitute( dict(CHAIN_SIZE=L, NNODES=nn) ) )

