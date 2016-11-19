nS="2 4 8"
for s in $nS; do
  rm -rf stages_${s}
  mkdir stages_${s}
  cd stages_${s}
  nTest=2
  nStep=1
  while test $nTest -gt 0; do
    echo "running: $nStep"
    ../../../../release/src/heat_rkc ${s} ${nStep} > log_${nStep}
    nStep=$(($nStep*2))
    nTest=$(($nTest-1))
  done
  nTest=7
  nStep=5
  while test $nTest -gt 0; do
    echo "running: $nStep"
    ../../../../release/src/heat_rkc ${s} ${nStep} > log_${nStep}
    nStep=$(($nStep*2))
    nTest=$(($nTest-1))
  done
  cd ..
done
