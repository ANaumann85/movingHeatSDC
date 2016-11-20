#!/bin/bash

Mv="2 3"
Pv="2 4 6 8"
iterV="1 2 3 4"
nStep="1 2 4 8"

for m in $Mv ; do
  for p in $Pv ; do
    for iter in $iterV; do
      for n in $nStep; do
        cdir="M_${m}_P_${p}_iter_${iter}_step_${n}"
	echo "cur: ${cdir}"
	rm -rf $cdir
        mkdir $cdir
	cd $cdir
	ln -s ../../sdc_quad_weights
	../../release/src/accuracy_model ${m} ${p} ${n} ${iter}
	cd ..
      done
    done
  done
done
