from diff import getDiff
from math import log

def order(e2h, eh):
    return [log(e2h[0]/eh[0])/log(2), log(e2h[1]/eh[1])/log(2)]

K=[2,3,4,8,32,64]
nStep=[2,4,8,16, 32, 64, 128]


baseN="T_sdc_K_%d_%d_fixed_%d.vtu"
refFile="T_ros2_20000_fixed_1.vtu"

for k in K:
    print "-------",k,"-------"
    dFile=baseN % (k, nStep[0], nStep[0])
    e2h=getDiff(refFile, dFile)
    print e2h, "-", "-"
    dFile1=baseN % (k, nStep[1], nStep[1])
    res=getDiff(refFile, dFile1)
    e2h2=getDiff(dFile, dFile1)
    print res, order(e2h, res), "-" #order(e2h2, res2)
    e2h=res
    dFile=dFile1
    for i in range(2,len(nStep)):
        n=nStep[i]
        dFile1=baseN % (k, n, n)
        res=getDiff(refFile, dFile1)
        res2=getDiff(dFile, dFile1)
        print res, order(e2h, res), order(e2h2, res2)
        e2h=res
        e2h2=res2
        dFile=dFile1
    print "------------------"


