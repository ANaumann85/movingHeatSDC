#!/bin/bash

#Submit this script with: sbatch thefilename

#SBATCH --time=60   # walltime in minutes
#SBATCH --ntasks=NCPU   # number of processor cores (i.e. tasks)
#SBATCH --partition=haswell
#SBATCH --mem-per-cpu=2048   # memory per CPU core
#SBATCH --output=logfile_pfasst_NCPU_%j
#SBATCH --mail-user=<MAIL>   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -A wir


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=1 #$SLURM_CPUS_ON_NODE
#OUTFILE="outName"
nTest=10
numIterV="1 2 3 4"
BINARY="../FE_heat2d_pfasst"
nStep=2
while test ${nTest} -gt 0; do
  for numIter in $numIterV; do
    srun ${BINARY}  --abs_res_tol 10e-10 --num_elements 10 --tend 20 --num_iters ${numIter}  --num_steps ${nStep} > log_2dMovingRect_${nStep}_${numIter}
    mv heat2dMoving_result_mlsdc.vtu heat2dMoving_result_mlsdc_${nStep}_${numIter}.vtu
  done
  nStep=$(($nStep*2))
  nTest=$(($nTest-1))
done

exit 0
