#!/bin/bash

nTest=10
n=5

while test $nTest -ge 0; do
  ../release/src/heat_sdc $n > log_${n}
  nTest=$(($nTest-1))
  n=$(($n*2))
done
