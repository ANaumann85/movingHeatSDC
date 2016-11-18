nTest=14
nStep=5
while test $nTest -gt 0; do
  echo "running: $nStep"
  ../release/src/heat_ros2 ${nStep} > log_${nStep}
  nStep=$(($nStep*2))
  nTest=$(($nTest-1))
done
