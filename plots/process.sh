#!/bin/bash

FOLDER=/storage/WRFOUT/Spain6_1

for file in `ls ${FOLDER}/wrfout_d01*`
do
   file1=`echo $file | sed 's/d01/d02/'`
   echo $file
   time (python3 post_process.py $file & (sleep $[ ( $RANDOM % 5 ) + 1 ] && python3 post_process.py $file1))
   # mv $file ${FOLDER}/processed/
done
