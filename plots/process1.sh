#!/bin/bash
RUN_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

FOLDER=/storage/WRFOUT/Spain6_1


while [ ! -f $RUN_DIR/STOP ]
do
   echo "Processing The following files:"
   ls ${FOLDER}/wrfout_d01*
   for file in `ls ${FOLDER}/wrfout_d01*`
   do
      sleep 10   # wait 10 seconds in case the files are being written
      file1=`echo $file | sed 's/d01/d02/'`
      echo $file
      time (python3 post_process.py $file & (sleep $[ ( $RANDOM % 5 ) + 1 ] && python3 post_process.py $file1))
      # mv $file ${FOLDER}/processed/
   done
   sleep 1m
done
