from diff import getDiff
from math import log

def order(e2h, eh):
    return [log(e2h[0]/eh[0])/log(2), log(e2h[1]/eh[1])/log(2)]

nTest=16
nSteps=40
nStepsRef=40*(2**nTest)
refFile="heat_ros2_%d-00001.vtu" % nStepsRef
dFile="heat_ros2_%d-00001.vtu" % (40)
e2h=getDiff(refFile, dFile)
print e2h, "-", "-"
dFile1="heat_ros2_%d-00001.vtu" % (40*2)
res=getDiff(refFile, dFile1)
e2h2=getDiff(dFile, dFile1)
print res, order(e2h, res), "-" #order(e2h2, res2)
e2h=res
dFile=dFile1

for n in range(2, nTest):
    dFile1="heat_ros2_%d-00001.vtu" % (40*(2**n))
    res=getDiff(refFile, dFile1)
    res2=getDiff(dFile, dFile1)
    print res, order(e2h, res), order(e2h2, res2)
    e2h=res
    e2h2=res2
    dFile=dFile1


