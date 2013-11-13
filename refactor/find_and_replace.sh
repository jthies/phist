#!/bin/bash
REPLACE=_IS_COMPLEX_
REPL_BY=IS_COMPLEX

SED_STRING=s/${REPLACE}/${REPL_BY}/g
echo $SED_STRING > tmp.sed

find . -name "*.c" -exec sed -i -f tmp.sed {} \;
find . -name "*.c.in" -exec sed -i -f tmp.sed {} \;
find . -name "*.h" -exec sed -i -f tmp.sed {} \;
find . -name "*.cpp" -exec sed -i -f tmp.sed {} \;
find . -name "*.hpp" -exec sed -i -f tmp.sed {} \;

rm tmp.sed
