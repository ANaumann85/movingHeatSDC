from diff import getDiff
from math import log

te=20.0
nTest=9
nStepsRef=5*(2**nTest)
refFile="heat_sdc_M-3_P-2_nStep-%d-00001.vtu" % nStepsRef
refFile="../ros2/heat_ros2_40960-00001.vtu"

print "#h l2 max"
for n in range(0, nTest):
    nstep=(5*(2**n))
    dFile1="heat_sdc_M-3_P-2_nStep-%d-00001.vtu" % nstep
    res=getDiff(refFile, dFile1)
    h=te/nstep
    print h, res[0], res[1] 
