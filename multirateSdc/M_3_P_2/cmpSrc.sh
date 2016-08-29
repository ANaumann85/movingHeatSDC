#/bin/bash

K="2 3 4 8 32 64"
nStep="2 4 8 16 32 64 128"

baseN="T_sdc_K_%d_%d_fixed_%d.vtu"
ref="T_ros2_20000_fixed_1.vtu"
#ref="T_sdc_K_32_64_fixed_64.vtu"
callB="python ../diff.py ${ref} ${baseN}"
for k in $K; do
  for n in $nStep; do
    curF=$(printf ${baseN} $k $n $n)
    if test -e  ${curF}; then
    call=$(printf "${callB}" $k $n $n);
    res=$(${call} | head -n 1)
    echo "$k $n $res" | sed 's/maxNorm://'
    fi
  done
done
