#!/bin/bash
nTest=14

nStep=1
while test $nTest -gt 0; do  
  ../release/src/heat_pfasst_ros2 $nStep > log_${nStep}
  nStep=$(($nStep*2))
  nTest=$(($nTest-1))
done
