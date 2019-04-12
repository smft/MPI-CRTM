#!/usr/bin/env bash
cd libsrc
mpiifort -O3 -xAVX -c -traceback -g module_calculate_tbb.f90 -I${CRTM}/include
mpiifort -O3 -xAVX -c -traceback -g main.f90 -I${NETCDF}/include -I${HDF5}/include -I${ZLIB}/include -I${CRTM}/include
cd ../bin
mpiifort -xAVX -traceback -g ../libsrc/module_calculate_tbb.o ../libsrc/main.o -o CRTM.exe -L${NETCDF}/lib -lnetcdff -lnetcdf -L${HDF5}/lib -lhdf5_hl -lhdf5 -L${ZLIB}/zlib -lz -L${CRTM}/lib -lcrtm






