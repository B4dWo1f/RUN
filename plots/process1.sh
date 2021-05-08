#!/bin/bash
RUN_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

FOLDER=/storage/WRFOUT/Spain6_1


while [ ! -f $RUN_DIR/STOP ]
do
   echo "Processing The following files:"
   for file in `ls ${FOLDER}/wrfout_d01*`
   do
      sleep 1.5m   # wait 10 seconds in case the files are being written
      file1=`echo $file | sed 's/d01/d02/'`
      ls $file
      ls $file1
      time (python3 post_process.py $file & python3 post_process.py $file1)
      # mv $file ${FOLDER}/processed/
   done
   sleep 1.5m
done
