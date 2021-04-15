#!/bin/bash

FOLDER_INSTALL=$HOME/METEO
Ncores=32

export DIR=$FOLDER_INSTALL/Build_WRF/LIBRARIES
export CC=gcc
export CXX=g++
export FC=gfortran
export FCFLAGS=-m64
export F77=gfortran
export FFLAGS=-m64
export JASPERLIB=$DIR/grib2/lib
export JASPERINC=$DIR/grib2/include
export LDFLAGS=-L$DIR/grib2/lib
export CPPFLAGS=-I$DIR/grib2/include
export PATH=$DIR/netcdf/bin:$PATH    #XXX  to .bashrc??
export NETCDF=$DIR/netcdf            #XXX  to .bashrc??
export PATH=$DIR/mpich/bin:$PATH    #XXX  to .bashrc??
export LD_LIBRARY_PATH=$DIR/grib2/lib:$LD_LIBRARY_PATH

# echo "expected PATH:"
# echo "/home/aeolus/RASP/Build_WRF/LIBRARIES/mpich/bin:/home/aeolus/RASP/Build_WRF/LIBRARIES/netcdf/bin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games"
# echo
# echo "current PATH:"
# echo $PATH
