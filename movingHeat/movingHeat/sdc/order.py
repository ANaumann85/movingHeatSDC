from diff import getDiff
from math import log

def order(e2h, eh):
    return [log(e2h[0]/eh[0])/log(2), log(e2h[1]/eh[1])/log(2)]

nTest=9
nStepsRef=5*(2**nTest)
refFile="heat_sdc_M-3_P-2_nStep-%d-00001.vtu" % nStepsRef
#refFile="../ros2/heat_ros2_40960-00001.vtu"
dFile="heat_sdc_M-3_P-2_nStep-%d-00001.vtu" % (5)
e2h=getDiff(refFile, dFile)
print e2h, "-", "-"
dFile1="heat_sdc_M-3_P-2_nStep-%d-00001.vtu" % (5*2)
res=getDiff(refFile, dFile1)
e2h2=getDiff(dFile, dFile1)
print res, order(e2h, res), "-" #order(e2h2, res2)
e2h=res
dFile=dFile1

for n in range(2, nTest):
    dFile1="heat_sdc_M-3_P-2_nStep-%d-00001.vtu" % (5*(2**n))
    res=getDiff(refFile, dFile1)
    res2=getDiff(dFile, dFile1)
    print res, order(e2h, res), order(e2h2, res2)
    e2h=res
    e2h2=res2
    dFile=dFile1


