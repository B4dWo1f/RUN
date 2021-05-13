# System Environment Tests
## 1 gfortran cpp gcc
which gfortran
which cpp
which gcc

gcc --version 


## 2 Folders
FOLDER_INSTALL=$HOME/METEO
sudo apt-get install gfortran cpp gcc
mkdir $FOLDER_INSTALL/Build_WRF $FOLDER_INSTALL/TESTS

## 3 Tests
cd $FOLDER_INSTALL/TESTS
wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/Fortran_C_tests.tar
tar -xf Fortran_C_tests.tar

gfortran TEST_1_fortran_only_fixed.f
./a.out
if [ $? -ne 0 ] then
 echo "Error in TEST_1_fortran_only_fixed.f"
fi

gfortran TEST_2_fortran_only_free.f90
./a.out

gcc TEST_3_c_only.c
./a.out

gcc -c -m64 TEST_4_fortran+c_c.c
gfortran -c -m64 TEST_4_fortran+c_f.f90
gfortran -m64 TEST_4_fortran+c_f.o TEST_4_fortran+c_c.o
./a.out

## 4
./TEST_csh.csh
## Fails
sudo apt-get install csh
./TEST_csh.csh
./TEST_perl.pl
./TEST_sh.sh

## 5
which ar
which head
which sed
which awk
which hostname
which sleep
which cat
which ln
which sort
which cd
which ls
which tar
which cp
which make
which touch
which cut
which mkdir
which tr
which expr
which mv
which uname
which file
which nm
which wc
which grep
which printf
which which
which gzip
which rm
which m4 # missing
sudo apt-get install m4
which m4

# 2 Building the libraries
cd ../Build_WRF
mkdir LIBRARIES
echo "Download libraries to $FOLDER_INSTALL/Build_WRF/LIBRARIES"
echo "https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compilation_tutorial.php#STEP2"
echo "Press any key to continue"
wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/mpich-3.0.4.tar.gz
wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/netcdf-4.1.3.tar.gz
wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/jasper-1.900.1.tar.gz
wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/libpng-1.2.50.tar.gz
wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/zlib-1.2.7.tar.gz



read a
## Export variables
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

## NetCDF
tar xzvf netcdf-4.1.3.tar.gz
cd netcdf-4.1.3
./configure --prefix=$DIR/netcdf --disable-dap --disable-netcdf-4 --disable-shared
make
make install
export PATH=$DIR/netcdf/bin:$PATH    #XXX  to .bashrc??
export NETCDF=$DIR/netcdf            #XXX  to .bashrc??
cd .. 

## MPICH
tar xzvf mpich-3.0.4.tar.gz     #or just .tar if no .gz present
cd mpich-3.0.4
./configure --prefix=$DIR/mpich
make
make install
export PATH=$DIR/mpich/bin:$PATH    #XXX  to .bashrc??
cd .. 

## zlib
tar xzvf zlib-1.2.7.tar.gz     #or just .tar if no .gz present
cd zlib-1.2.7
./configure --prefix=$DIR/grib2
make
make install
cd .. 

## libpng
tar xzvf libpng-1.2.50.tar.gz     #or just .tar if no .gz present
cd libpng-1.2.50
./configure --prefix=$DIR/grib2
make
make install
cd .. 

## Jasper
tar xzvf jasper-1.900.1.tar.gz     #or just .tar if no .gz present
cd jasper-1.900.1
./configure --prefix=$DIR/grib2
make
make install
cd .. 

# XXX At this point PATH should be:
# export PATH=/home/aeolus/RASP/Build_WRF/LIBRARIES/mpich/bin:/home/aeolus/RASP/Build_WRF/LIBRARIES/netcdf/bin:/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games


echo "Download tests to $FOLDER_INSTALL/TESTS"
echo "https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compilation_tutorial.php#STEP3"
echo "Press any key to continue"
read a

cd $FOLDER_INSTALL/TESTS
wget https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compile_tutorial/tar_files/Fortran_C_NETCDF_MPI_tests.tar


## Test 1 Fortran + C + NetCDF
tar -xf Fortran_C_NETCDF_MPI_tests.tar
cp ${NETCDF}/include/netcdf.inc .
gfortran -c 01_fortran+c+netcdf_f.f
gcc -c 01_fortran+c+netcdf_c.c
gfortran 01_fortran+c+netcdf_f.o 01_fortran+c+netcdf_c.o -L${NETCDF}/lib -lnetcdff -lnetcdf
./a.out
## Test 2 Fortran + C + NetCDF + MPI
cp ${NETCDF}/include/netcdf.inc .
mpif90 -c 02_fortran+c+netcdf+mpi_f.f
mpicc -c 02_fortran+c+netcdf+mpi_c.c
mpif90 02_fortran+c+netcdf+mpi_f.o 02_fortran+c+netcdf+mpi_c.o -L${NETCDF}/lib -lnetcdff -lnetcdf
mpirun ./a.out

# Building WRF
cd $FOLDER_INSTALL
git clone https://github.com/wrf-model/WRF
cd WRF
echo "when asked select 34. (dmpar) for GNU (gfortran/gcc)"
echo "and select nesting: 1=basic"
./configure
./compile em_real >& log.compile
ls -ls main/*.exe

# Building WPS
cd $FOLDER_INSTALL
git clone https://github.com/wrf-model/WPS 
cd WPS
./clean
echo "when asked select option 1.  Linux x86_64, gfortran (serial)"
echo "it is recommended to use serial regardless of the selection in WRF"
./configure
./compile >& log.compile
ls -ls *.exe
ls -ls geogrid/src/geogrid.exe
ls -ls metgrid/src/metgrid.exe
ls -ls ungrib/src/ungrib.exe


sudo apt-get install python3-rasterio python3-cartopy python3-netcdf4 python3-xarray
sudo pip3 install wrf-python
