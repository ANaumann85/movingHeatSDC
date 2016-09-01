from scipy.io import *
from scipy.sparse import linalg as sLA
from scipy import linalg as dLA
import numpy as np
import matplotlib.pyplot as plt
lapl=mmread('laplace.mm') 
mass=mmread('mass.mm')
nEig=lapl.shape[0]-1
[eigL, vL]=sLA.eigsh(lapl, nEig, mass)

nT=5
te=15
dt=te/nT
allEigLmv=[]
for k in range(0, nT):
    t=k*dt
    lmvName='lmv-%0.2f.mm' % t
    laplMv=mmread(lmvName)
    #[eigLmv, vL]=sLA.eigsh(laplMv, nEig, mass)
    [eigLmv, vL]=dLA.eig(laplMv.todense(),mass.todense()) #,right=False)
    allEigLmv.append(np.sort(eigLmv))

for k in range(1, nT):
    print max(abs(np.divide(allEigLmv[0]-allEigLmv[k],allEigLmv[0])))
x=range(0,nEig)
fig=plt.figure()
plt.plot(eigL,'k-')
for k in range(0, nT):
    plt.plot(allEigLmv[k],'-')
plt.show()
