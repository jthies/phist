#!/bin/bash
# as we can't write an empty matrix in HB format, we fix it like this:
find . -name "*spzero*.*ua" -exec sed -i s/1.00/0.00/ {} \;

# put 'general' into the speye and spzero MatrixMarket files as we can't usually read symmetric matrices by half
find . -name "*spzero*.mm" -or -name "*speye*.mm" -exec sed -i s/symmetric/general/ {} \;

