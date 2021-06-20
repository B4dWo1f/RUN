#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

#
# This script provides generate_config.py the arguments to generate a
# config gile to run run_rasp.sh
# Arguments:
#   $1 = domain
#   $2 = day
#   $3 = hours
#
# Usage:
#   launch_rasp.sh <domain> <+day> <h0,h1>
# Example:
#   launch_rasp.sh Spain6_1 1 8,9
#

ID=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1)

FOL="/tmp/METEO_$ID"   # temporary folder to run


mkdir -p ${FOL}
(
cd ${FOL}
mkdir -p RUN WRF/run WPS dataGFS runtime
ln -s /home/aeolus/METEO/RUN/* /tmp/METEO_$ID/RUN/
ln -s /home/aeolus/METEO/WRF/* /tmp/METEO_$ID/WRF/
ln -s /home/aeolus/METEO/WRF/run/* /tmp/METEO_$ID/WRF/run/
ln -s /home/aeolus/METEO/WPS/* /tmp/METEO_$ID/WPS/
)

(
cd ${FOL}/RUN/
echo > run_rasp.err
echo > run_rasp.log

./generate_config.py ${FOL} $1 $2 $3

time ./run_rasp.sh ${FOL} $2 $3
)

rm -r ${FOL}   # delete folder (XXX no debugging!!!)
