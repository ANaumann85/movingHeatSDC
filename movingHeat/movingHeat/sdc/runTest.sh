#!/bin/bash

Mv="2 3"
Pv="2 3 6 8"
for M in $Mv; do
for P in $Pv; do
nTest=10
n=5
while test $nTest -ge 0; do
  ../../release/src/heat_sdc $n $M $P > log_${n}_${M}_${P}
  nTest=$(($nTest-1))
  n=$(($n*2))
done
done
done
