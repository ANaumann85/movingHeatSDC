#!/bin/bash

K="2 3 4 8 32 64"
#nstep="2 4 8"
nstep="2 4 8 16 32 64 128"

for k in $K; do
  for n in $nstep; do
    python runsdc_drake.py $n $k
  done
done
