#!/bin/bash

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

DOMAIN=`grep domain_folder $FOLDER_INSTALL/RASP/config.ini  | cut -d '=' -f 2 | sed 's/ //g'`
Ncores=32

cd $FOLDER_INSTALL

rm WPS/namelist.wps WPS/namelist.input
rm WPS/FILE* WPS/GRIBFILE.AA* WPS/*.log WPS/log.*
rm WRF/run/namelist.input 
rm WRF/run/rsl.* WRF/run/wrfout* WRF/run/met_em*
rm dataGFS/*
rm $FOLDER_INSTALL/RASP/namelist.wps $FOLDER_INSTALL/RASP/namelist.input
rm $DOMAIN/met_em* $DOMAIN/geo_em.*



# ln -s $DOMAIN/namelist.wps $FOLDER_INSTALL/WPS/
# ln -s $DOMAIN/namelist.input $FOLDER_INSTALL/WPS/
# ln -s $DOMAIN/namelist.input $FOLDER_INSTALL/WRF/run/

#vim $DOMAIN/namelist.input
#vim $DOMAIN/namelist.wps

cd $FOLDER_INSTALL/RASP
rm namelist.*
python3 inputer.py
ln -s $FOLDER_INSTALL/RASP/namelist.wps $FOLDER_INSTALL/WPS/
ln -s $FOLDER_INSTALL/RASP/namelist.input $FOLDER_INSTALL/WPS/
ln -s $FOLDER_INSTALL/RASP/namelist.input $FOLDER_INSTALL/WRF/run/

time python3 get_gfs_data.py
ls $FOLDER_INSTALL/dataGFS

cd $FOLDER_INSTALL/WPS

time ./geogrid.exe >& log.geogrid && tail log.geogrid
echo "Geogrid geo_met* files:"
ls $DOMAIN/geo_em*
./link_grib.csh ../dataGFS/

ln -sf ungrib/Variable_Tables/Vtable.GFS Vtable
time ./ungrib.exe
echo "Ungrib 'FILES*' files:"
ls FILE*

time ./metgrid.exe >& log.metgrid && tail log.metgrid
echo 'Metgrid met_em* files:'
ls $DOMAIN/met_em*

cd ../WRF/run
ln -sf $DOMAIN/met_em* .
ls met_em*
time mpirun -np 1 ./real.exe
tail -n 1 rsl.error.0000 | grep -w SUCCESS
if [ $? -ne 0 ]
then
   echo "Error running real.exe"
else
   echo "REAL worked!!"
fi


time mpirun -np $Ncores ./wrf.exe
tail -n 1 rsl.error.0000 | grep -w SUCCESS
if [ $? -ne 0 ]
then
   echo "Error running wrf.exe"
else
   echo "WRF worked!!"
fi

mv wrfout* $FOLDER_INSTALL/RASP/mwe_process/data/
