#!/bin/bash

#Load modules
module --force purge
module load modules/1.59-20220201
module load intel-oneapi-compilers intel-oneapi-mkl intel-mpi/2017.4.196 fftw/3.3.10-mpi hdf5/1.12.1-mpi
export MPICH_CXX=icpc

# Run configuration script and save to logfile
./configure.py --prob=my_setup_constant --cxx=icpc -fft -mpi -hdf5 -h5double

# Compile
make clean
make -j 16
