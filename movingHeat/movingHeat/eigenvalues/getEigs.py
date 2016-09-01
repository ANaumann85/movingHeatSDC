from scipy.io import *
from scipy.sparse.linalg import *
import matplotlib.pyplot as plt
lapl=mmread('laplace.mm') 
mass=mmread('mass.mm')
nEig=lapl.shape[0]-1
[eigL, vL]=eigsh(lapl, nEig, mass)

nT=5
te=15
dt=te/nT
allEigLmv=[]
for k in range(0, nT):
    t=k*dt
    lmvName='lmv-%0.2f.mm' % t
    laplMv=mmread(lmvName)
    [eigLmv, vL]=eigsh(laplMv, nEig, mass)
    allEigLmv.append(eigLmv)

x=range(0,nEig)
for k in range(0, nT):
    plt.plot(x, eigL, x, allEigLmv[k])
plt.show()
