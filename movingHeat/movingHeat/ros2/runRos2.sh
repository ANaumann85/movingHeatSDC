nTest=10
nStep=40
while test $nTest -gt 0; do
  nStep=$(($nStep*2))
  nTest=$(($nTest-1))
  echo "running: $nStep"
  ../release/src/heat_ros2 ${nStep} > log_${nStep}
done
