#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

#
# This script provides generate_config.py the arguments to generate a
# config file to run run_rasp.sh
# Arguments:
#   $1 = domain
#   $2 = day
#   $3 = hours
#
# Usage:
#   launch_rasp.sh <domain> <+day> <h0,h1>
# Example:
#   launch_rasp.sh Spain6_1 1 8,20
# this command would run WRF for tomorrow since 9:00 till 20:00
#

ID=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1)

FOL="/tmp/METEO_$ID"   # temporary folder to run

killall wrf.exe 2> /dev/null
rm $HOME/errors.txt
touch $HOME/errors.txt

mkdir -p ${FOL}
(
cd ${FOL}
mkdir -p RUN WRF/run WPS dataGFS
#runtime
ln -s $HOME/METEO/RUN/* /tmp/METEO_$ID/RUN/
ln -s $HOME/METEO/WPS/* /tmp/METEO_$ID/WPS/
ln -s $HOME/METEO/WRF/run/* /tmp/METEO_$ID/WRF/run/
ln -s $HOME/METEO/WRF/* /tmp/METEO_$ID/WRF/
)

(
cd ${FOL}/RUN/
echo > run_rasp.err
echo > run_rasp.log

./generate_config.py ${FOL} $1 $2 $3

time ./run_rasp.sh ${FOL} $2 $3
)

rm -r ${FOL}   # delete folder (XXX no debugging!!!)

# Store CFL error files
if [ -s $HOME/errors.txt ]
then
   echo "`date` there were CFL errors"
   cp $HOME/errors.txt $HOME/errors_`date '+%Y%m%d_%H%M'`.txt
else
   echo "`date` NO CFL errors"
fi

