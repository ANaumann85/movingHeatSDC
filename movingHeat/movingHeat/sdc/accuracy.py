from diff import getDiff
from math import log

te=20.0
nTest=9
nStepsRef=5*(2**nTest)
P_vec=[2,3,6,8]
refFile="../kiter_5/heat_sdc_M-3_P-8_nStep-%d-00001.vtu" % nStepsRef
#refFile="../ros2/heat_ros2_40960-00001.vtu"

for P in P_vec:
    acname="accuracy_P_%d.dat" %P
    f=open(acname,'w')
    f.write("#h l2 max\n")
    for n in range(0, nTest):
        nstep=(5*(2**n))
        dFile1="heat_sdc_M-3_P-%d_nStep-%d-00001.vtu" % (P,nstep)
        res=getDiff(refFile, dFile1)
        h=te/nstep
        data=str(h)+" " +str(res[0])+" " +str(res[1])
        f.write(data+"\n")
    f.close()
