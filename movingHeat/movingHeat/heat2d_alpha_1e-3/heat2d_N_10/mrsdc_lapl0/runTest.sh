#!/bin/bash

Mv="2 3"
Pv="2 3 6 8"
kv="1 2 3 4"
for k in $kv; do
cd kiter_${k}
for M in $Mv; do
for P in $Pv; do
nTest=1
n=1
while test $nTest -ge 0; do
  ../../release/src/heat_sdc $n $M $P $k > log_${n}_${M}_${P}
  nTest=$(($nTest-1))
  n=$(($n*2))
done
nTest=8
n=5
while test $nTest -ge 0; do
  ../../release/src/heat_sdc $n $M $P $k > log_${n}_${M}_${P}
  nTest=$(($nTest-1))
  n=$(($n*2))
done
done
done
cd ..
done
