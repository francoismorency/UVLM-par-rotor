#!/bin/bash
#SBATCH --account=def-morency
#SBATCH --time=10:00:00
#SBATCH --constraint=skylake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48  
#SBATCH --mem=2000M
#SBATCH --job-name=6
#SBATCH --output=buffer.out
#SBATCH --mail-user=<guillaume.damour.1@ens.etsmtl.ca>
#SBATCH --mail-type=ALL

module load matlab/2022a
MATLABHOME=/cvmfs/restricted.computecanada.ca/easybuild/software/2020/Core/matlab/2022a
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MATLABHOME/bin/glnxa64
make -f Makefile_cedar.mk clean
make -f Makefile_cedar.mk MATLABHOME=$MATLABHOME
export OMP_NUM_THREADS=6
./UVLM_2_Blade_Rotor_Axial 6
