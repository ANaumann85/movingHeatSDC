#!/bin/bash

Mv="2 3"
for M in $Mv; do
nTest=10
n=5
while test $nTest -ge 0; do
  ../../release/src/heat_sdc_standard $n $M > log_${n}_${M}
  nTest=$(($nTest-1))
  n=$(($n*2))
done
done
