#!/bin/bash

FOLDER=/storage/WRFOUT/Spain6_1

mkdir -p ${FOLDER}/processed
for file in `ls ${FOLDER}/wrfout_*`
do
   echo $file
   time python3 post_process.py $file
   # mv $file ${FOLDER}/processed/
done
