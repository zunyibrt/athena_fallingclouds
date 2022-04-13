#!/bin/bash

python3 configure.py --prob=my_setup_constant -mpi -hdf5 -h5double --include=$TACC_HDF5_INC --lib_path=$TACC_HDF5_LIB --cxx=icc-phi
