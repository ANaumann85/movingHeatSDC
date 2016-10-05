from diff import getDiff
from math import log

te=20.0
nTest=13
nStepsRef=5*(2**nTest)
refFile="heat_ros2_%d-00001.vtu" % nStepsRef

print "#h l2 max"
for n in range(0, nTest):
    nstep=(5*(2**n))
    h=te/nstep
    dFile1="heat_ros2_%d-00001.vtu" % nstep
    res=getDiff(refFile, dFile1)
    print h,res[0],res[1]

