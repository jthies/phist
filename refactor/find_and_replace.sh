#!/bin/bash
REPLACE=IS_COMPLEX
REPL_BY=PHIST_COMPLEX

SED_STRING=s/${REPLACE}/${REPL_BY}/g
echo $SED_STRING > tmp.sed

find . -name "*.c" -exec sed -i -f tmp.sed {} \;
find . -name "*.c.in" -exec sed -i -f tmp.sed {} \;
find . -name "*.h" -exec sed -i -f tmp.sed {} \;
find . -name "*.cpp" -exec sed -i -f tmp.sed {} \;
find . -name "*.hpp" -exec sed -i -f tmp.sed {} \;

rm tmp.sed
