#!/bin/bash

nCores="2 4 8 16"

base="jobSkript_taurus_base.sh"
for n in $nCores; do
  curDir="pfasst_${n}"
  rm -rf $curDir
  mkdir $curDir
  sed "s:NCPU:${n}:g" $base > ${curDir}/jobSkript_taurus.sh
  cd $curDir
#sbatch jobSkript_taurus.sh
  cd ..
done
