#!/bin/sh

MATLABHOME=/home/guillaume/MATLAB/R2022a

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MATLABHOME/bin/glnxa64
make -f Makefile_own_computer.mk clean
make -f Makefile_own_computer.mk MATLABHOME=$MATLABHOME
export OMP_NUM_THREADS=$1
./UVLM_2_Blade_Rotor_Axial $1
