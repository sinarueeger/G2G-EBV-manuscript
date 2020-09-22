#!/bin/bash

STRING=$1 ## takes first argument
DIR_SCRATCH=$2 ## takes first argument
DEBUGGING=$3 ## debugging?

FILES=$(ls $DIR_SCRATCH/lmm_$TOOL_*$STRING*.mlma | grep -v 'subset')

## use GNU parallel
##  | parallel --pipe
## https://stackoverflow.com/questions/20308443/parallel-processing-in-awk

for FILENAME in $FILES
do

  NEWFILENAME="$FILENAME-subset.txt"
  
  if [ ! -f $NEWFILENAME ]; then
    echo $NEWFILENAME
   
    if [ "$DEBUGGING" == "TRUE" ]; then
     # awk '{ if ($9 < 1) { print } }' $FILENAME > $NEWFILENAME &
      cp $FILENAME $NEWFILENAME
      echo debugging
    fi

    if [ "$DEBUGGING" == "FALSE" ]; then
      awk '{ if ($9 < 5e-8) { print } }' $FILENAME > $NEWFILENAME &
      echo nodebugging
    fi
    
  fi
  
  
done
