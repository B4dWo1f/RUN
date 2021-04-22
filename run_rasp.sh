#!/bin/bash
RUN_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# STDOUT to log
LOG_FILE="${RUN_DIR}/run_rasp.log"
echo > "${LOG_FILE}"
# Close STDOUT file descriptor
exec 1<&-
# Open STDOUT as $LOG_FILE file for read and write.
exec 1<>"${LOG_FILE}"

# STDERR to err
ERR_FILE="${RUN_DIR}/run_rasp.err"
echo > "${ERR_FILE}"
# Close STDERR FD
exec 2<&-
# Redirect STDERR to STDOUT
exec 2<>"${ERR_FILE}"


# changes in namelist.input
# run_days, run_hours...
# start_year,month,day,hour...
# end_year,month,day,hour...
# interval_seconds  3600
# history  60 60
# frames_per_outfile 1, 1
# num_metgrid_levels       = 34,


# changes in namelist.wps
# start_date
# end_date
# interval_seconds


source rasp_env.sh 


DOMAIN=`./get_domain.py | head -n 1 | tail -n 1`
GFSdata=`./get_domain.py | head -n 2 | tail -n 1`
OUTdata=`./get_domain.py | head -n 3 | tail -n 1`
Ncores=`./get_domain.py | head -n 4 | tail -n 1`

echo -e "Starting RUN"
echo -e "Domain         : ${DOMAIN}"
echo -e "GFS Data Folder: ${GFSdata} "
echo -e "Ncores: ${Ncores}\n"

(
echo "Cleaning up previous runs"
cd "${FOLDER_INSTALL}"
#### Clean previous runs
rm WPS/namelist.wps WPS/namelist.input
rm WPS/FILE* WPS/GRIBFILE.AA* WPS/*.log WPS/log.*
rm WRF/run/namelist.input
rm WRF/run/rsl.* WRF/run/wrfout* WRF/run/met_em*
rm "${GFSdata}"/*
rm "${RUN_DIR}"/namelist.wps "${RUN_DIR}"/namelist.input
rm "${DOMAIN}"/met_em* "${DOMAIN}"/geo_em.*
)

(
cd "${RUN_DIR}"
echo "Setting up the inputs for RUN"
#### Prepare namelists
rm namelist.*
python3 inputer.py
ln -s "${RUN_DIR}/namelist.wps" "${FOLDER_INSTALL}/WPS/"
ln -s "${RUN_DIR}/namelist.input" "${FOLDER_INSTALL}/WPS/"
ln -s "${RUN_DIR}/namelist.input" "${FOLDER_INSTALL}/WRF/run/"
echo "WPS/WRF Input files:"
ls ${FOLDER_INSTALL}/WPS/namelist.*
ls ${FOLDER_INSTALL}/WRF/run/namelist.*

#### Download GFS data
echo "Downloading GFS data"
time python3 download_gfs_data.py
if [ $? -eq 0 ]; then
   echo "GFS data downloaded to ${GFSdata}:"
   ls "${GFSdata}"
else
   1>&2 echo "Error downloading GFS data"
   exit 1
fi
)
# Check Input and GFS section
if [ $? -eq 0 ]; then
   echo_success 'Input & GFS'
else
   echo_fail "Error in Input or GFS data"
   exit 1
fi


(
#### WPS
cd "$FOLDER_INSTALL/WPS"
echo "Running geogrid"
time ./geogrid.exe >& log.geogrid
grep "Successful completion of geogrid" log.geogrid
if [ $? -eq 0 ]; then
   echo "Geogrid was successful"
   echo "Geogrid geo_met* files:"
   ls $DOMAIN/geo_em*
   ./link_grib.csh ../dataGFS/
else
   1>&2 echo "Error running Geogrid"
   exit 1
fi

echo "Running ungrib"
ln -sf ungrib/Variable_Tables/Vtable.GFS Vtable
time ./ungrib.exe
echo "Ungrib 'FILES*' files:"
ls FILE*
# XXX ADD checking point

echo "Running metgrid"
time ./metgrid.exe >& log.metgrid && tail log.metgrid
grep -i "Successful completion of program metgrid.exe" metgrid.log
if [ $? -eq 0 ]; then
   echo "Metgrid was successful"
   echo 'Metgrid met_em* files:'
   ls $DOMAIN/met_em*
else
   1>&2 echo "Error running Metgrid"
   exit 1
fi
)
# Check WPS codes
if [ $? -eq 0 ]; then
   echo_sucess "WPS run Ok."
else
   echo_faile "Error during WPS steps"
   exit 1
fi


(
#### WRF
echo "Going for WRF"
cd $FOLDER_INSTALL/WRF/run
ln -sf $DOMAIN/met_em* .
echo "met* files are present:"
ls met_em*
echo -e "\nStarting real.exe"
time mpirun -np 1 ./real.exe
tail -n 1 rsl.error.0000 | grep -w SUCCESS
if [ $? -eq 0 ]; then
   echo "REAL worked!!"
else
   1>&2 echo "Error running real.exe"
fi

# WRF
echo -e "\nStarting wrf.exe"
time mpirun -np $Ncores ./wrf.exe
tail -n 1 rsl.error.0000 | grep -w SUCCESS
if [ $? -eq 0 ]; then
   echo "WRF worked!!"
else
   1>&2 echo "Error running wrf.exe"
fi

mkdir -p "${OUTdata}"
rm wrfoutReady*
# mv wrfout_* "${OUTdata}"
)
# Check WRF
if [ $? -eq 0 ]; then
   echo_success "WRF run Ok."
else
   echo_fail "Error during WRF steps"
   exit 1
fi
