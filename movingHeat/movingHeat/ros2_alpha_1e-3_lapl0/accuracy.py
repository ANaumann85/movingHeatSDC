from diff import getDiff
from math import log
from tail import *

te=20.0
nTest=10
nStepsRef=5*(2**nTest)
refFile="../ros2_alpha_1e-3/heat_ros2_40960-00001.vtu" 

print "#h l2 max nstep runtime[ms]"
for n in range(0, nTest):
    nstep=(5*(2**n))
    h=te/nstep
    dFile1="heat_ros2_%d-00001.vtu" % nstep
    res=getDiff(refFile, dFile1)
    log="log_"+str(nstep)
    logf=open(log)
    lastLine=tail(logf, 1)
    logf.close()
    print h,res[0],res[1],nstep,lastLine[0].split(':')[1].split()[0]

