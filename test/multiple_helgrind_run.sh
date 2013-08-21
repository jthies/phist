#!/bin/bash
#Shell script for running helgrind multiple times

#usage: ./multiple_helgrind_run.sh SCRIPT_PATH BUILD_PATH
#Helgrind output in LOG_FILE

#returns 0 if no error was found

if test $# -ne 2
then
   echo "wrong usage: ./"$basename" SCRIPT_PATH BUILD_PATH"
   exit 3
fi

SCRIPT_PATH=$1
BUILD_PATH=$2
LOG_FILE=$SCRIPT_PATH"/helgrind.log"
REPEAT=3

if test -f $LOG_FILE; then
      echo "deleting "$LOG_FILE
      rm -rf $LOG_FILE
fi  
for((var=0;var<$REPEAT;var++))
do
   echo "running valgrind #"$var
   valgrind --tool=helgrind --suppressions=$SCRIPT_PATH/helgrind_for_main_task_model.supp $BUILD_PATH/main_task_model 2>> $LOG_FILE	
   if test $? -ne 0
   then
      echo "error while running valgrind" >&2
      exit 2
   fi
done

#parsing helgrind.log
matches_race=$(grep -c "Possible data race during" $LOG_FILE)
matches_lock=$(grep -c "Lock at" $LOG_FILE)
if test $matches_race -gt 0
then
    echo "Found possible data races" >&2
    if test $matches_lock -gt 0
    then
       echo "Found possible dead locks" >&2
    fi
    echo "See "$LOG_FILE" for more information" >&2
    exit 1
fi
echo "no errors found"
exit 0
