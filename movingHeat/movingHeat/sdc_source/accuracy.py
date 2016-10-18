from diff import getDiff
from math import log
from tail import *

te=20.0
nTest=9
nStepsRef=5*(2**nTest)
refFile="../kiter_5/heat_sdc_standard_M-3_nStep-5120-00001.vtu"
#refFile="../ros2/heat_ros2_40960-00001.vtu"
M=2
acname="accuracy_M_%d.dat" %(M)
f=open(acname,'w')
f.write("#h l2 max nstep runtime[ms]\n")
for n in range(0, nTest):
    nstep=(5*(2**n))
    dFile1="heat_sdc_standard_M-%d_nStep-%d-00001.vtu" % (M,nstep)
    res=getDiff(refFile, dFile1)
    h=te/nstep
    log="log_%d_%d"%(nstep,M)
    logf=open(log)
    lastLine=tail(logf, 1)
    logf.close()
    data=str(h)+" " +str(res[0])+" " +str(res[1]) +" "+str(nstep)+" " +str(lastLine[0].split(':')[1].split()[0])
    f.write(data+"\n")
f.close()
