#!/bin/bash
# as we can't write an empty matrix in HB format, we fix it like this:
find . -name *spzero*.*ua -exec sed -i s/1.00/0.00/ {} ;
