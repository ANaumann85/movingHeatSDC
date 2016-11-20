#!/bin/bash

Mv="2 3"
Pv="2 4 6 8"
iterV="1 2 3 4"
nStep="1 2 4 8"

Mv="2"
Pv="2"
iterV="1"
nStep="1"

for m in $Mv ; do
  for p in $Pv ; do
    rm -rf M_${m}_P_${p}
    mkdir M_${m}_P_${m}
    cd M_${m}_P_${m}
    for iter in $iterV; do
      for n in $nStep; do
        cdir="iter_${iter}_step_${n}"
	rm -rf $cdir
        mkdir $cdir
	cd $cdir
	  ln -s ../../../sdc_quad_weights
	  ../../../release/src/accuracy_model ${m} ${p} ${n} ${iter}
	cd ..
      done
    done
    cd ..
  done
done
