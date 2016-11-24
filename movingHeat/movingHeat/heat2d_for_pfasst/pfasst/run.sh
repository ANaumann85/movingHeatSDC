#!/bin/bash

nTest=10
numIterV="1 2 3 4"
BINARY="../FE_heat2d_pfasst"
nStep=2
while test ${nTest} -gt 0; do
  for numIter in $numIterV; do
    srun ${BINARY}  --abs_res_tol 10e-10 --num_elements 10 --tend 20 --num_iters ${numIter}  --num_steps ${nStep} > log_2dMovingRect_${nStep}_${numIter}
    mv heat2dMoving_result.vtu heat2dMoving_result_${nStep}_${numIter}.vtu
  done
  nStep=$(($nStep*2))
  nTest=$(($nTest-1))
done

