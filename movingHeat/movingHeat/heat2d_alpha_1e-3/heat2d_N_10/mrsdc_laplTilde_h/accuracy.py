from diff import getDiff
from math import log
from tail import *

te=20.0
nTest=6
steps=[1,2]
for n in range(nTest):
    steps.append(5*2**n)
P_vec=[2,3,6,8]
#refFile="../../sdc/kiter_5/heat_sdc_standard_M-3_nStep-5120-00001.vtu" 
refFile="../../ros2/heat_ros2_40960-00001.vtu"
Mv=[2,3]
for M in Mv:
    for P in P_vec:
        acname="accuracy_M_%d_P_%d.dat" %(M,P)
        f=open(acname,'w')
        f.write("#h l2 max nstep runtime[ms]\n")
        for nstep in steps:
            dFile1="heat_sdc_M-%d_P-%d_nStep-%d-00001.vtu" % (M,P,nstep)
            res=getDiff(refFile, dFile1)
            h=te/nstep
            log="log_%d_%d_%d"%(nstep,M,P)
            logf=open(log)
            lastLine=tail(logf, 1)
            logf.close()
            data=str(h)+" " +str(res[0])+" " +str(res[1]) +" "+str(nstep)+" " +str(lastLine[0].split(':')[1].split()[0])
            f.write(data+"\n")
        f.close()
