from diff import getDiff
from math import log

def order(e2h, eh):
    return [log(e2h[0]/eh[0])/log(2), log(e2h[1]/eh[1])/log(2)]

nStep=map(lambda x: 25*(2**x), range(0,7))


baseN="T_ros2_%d_fixed_1.vtu"
refFile=baseN % (2*nStep[-1])

dFile=baseN % (nStep[0])
e2h=getDiff(refFile, dFile)
print e2h, "-", "-"
dFile1=baseN % (nStep[1])
res=getDiff(refFile, dFile1)
e2h2=getDiff(dFile, dFile1)
print res, order(e2h, res), "-" #order(e2h2, res2)
e2h=res
dFile=dFile1
for i in range(2,len(nStep)):
    n=nStep[i]
    dFile1=baseN % (n)
    res=getDiff(refFile, dFile1)
    res2=getDiff(dFile, dFile1)
    print res, order(e2h, res), order(e2h2, res2)
    e2h=res
    e2h2=res2
    dFile=dFile1
