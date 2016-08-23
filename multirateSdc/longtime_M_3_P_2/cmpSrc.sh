#/bin/bash

K="2 3 4 8 32 64"
nStep="2 4 8 16 32 64 128"

baseN="T_sdc_K_%d_%d_fixed_%d.vtu"
ref="../T_imex_2000_fixed_2000.vtu"
callB="python ../diff.py ${ref} ${baseN}"
for k in $K; do
  for n in $nStep; do
    call=$(printf "${callB}" $k $n $n);
    res=$(${call} | tail -n 1)
    echo "$k $n $res"
  done
done
