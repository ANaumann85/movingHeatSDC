#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=60   # walltime in minutes
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --partition=haswell
#SBATCH --mem-per-cpu=2048   # memory per CPU core
#SBATCH --output=logfile_ros2_%j
#SBATCH --mail-user=<MAIL>   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -A wir


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=1 #$SLURM_CPUS_ON_NODE
#OUTFILE="outName"
source ./runTests.sh

exit 0
