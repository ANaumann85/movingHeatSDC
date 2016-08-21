#/bin/bash

K="2 3 4 8" # 32"
nStep="4 8"

baseN="T_noSrc_sdc_K_%d_%d_fixed_%d.vtu"
ref="../T_noSrc_imex_20000_fixed_20000.vtu"
callB="python ../diff.py ${ref} ${baseN}"
for k in $K; do
  for n in $nStep; do
    call=$(printf "${callB}" $k $n $n);
    res=$(${call} | tail -n 1)
    echo "$k $n $res"
  done
done
