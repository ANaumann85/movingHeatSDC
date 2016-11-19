#!/bin/bash

Mv="2 3"
kv="1 2 3 4"
for k in $kv; do
rm -rf kiter_${k}
mkdir kiter_${k}
cd kiter_${k}
ln -s ../../../../sdc_quad_weights/
for M in $Mv; do
nTest=1
n=1
while test $nTest -ge 0; do
  echo "running $n $M $k "
  ../../../../release/src/heat_sdc_standard $n $M $k 3 > log_${n}_${M}
  nTest=$(($nTest-1))
  n=$(($n*2))
done
nTest=6
n=5
while test $nTest -ge 0; do
  echo "running $n $M $k "
  ../../../../release/src/heat_sdc_standard $n $M $k 3 > log_${n}_${M}
  nTest=$(($nTest-1))
  n=$(($n*2))
done
done
cd ..
done
