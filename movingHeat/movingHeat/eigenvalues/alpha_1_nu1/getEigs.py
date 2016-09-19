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

print 'difference of eigenvalues at time t to eigenvalues at time t=0'
for k in range(1, nT):
    print max(abs(np.divide(allEigLmv[0]-allEigLmv[k],allEigLmv[0])))
x=range(0,nEig)
fig=plt.figure(figsize=(8,6))
plt.plot(eigL,'k-', label='lapl')
for k in range(0, nT):
    plt.plot(np.real(allEigLmv[k]),'-', label='t=%d' % (k*dt))
#plt.show()
ax1=plt.gca()
box = ax1.get_position()
ax1.set_position([box.x0, box.y0 + box.height * 0.1,
        box.width, box.height * 0.9])
plt.legend(loc='upper center', bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=3)
plt.savefig('eigenvalues.pdf',bbox_inches='tight')

fig=plt.figure(figsize=(8,6))
for k in range(1, nT):
    plt.plot(abs(np.divide(allEigLmv[0]-allEigLmv[k],allEigLmv[0])), label='t=%d' %(k*dt))
plt.legend(loc='upper center', bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=4)
plt.savefig('eigenvalues_reldiff0.pdf',bbox_inches='tight')

fig=plt.figure(figsize=(8,6))
for k in range(1, nT):
    plt.plot(abs(allEigLmv[0]-allEigLmv[k]), label='t=%d' %(k*dt))
plt.legend(loc='upper center', bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=4)
plt.savefig('eigenvalues_absdiff0.pdf',bbox_inches='tight')

fig=plt.figure(figsize=(8,6))
plt.plot(abs(allEigLmv[0]-allEigLmv[1]), label='t=%d' %(1*dt))
plt.legend(loc='upper center', bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=4)
plt.savefig('eigenvalues_absdiff_0-3.pdf',bbox_inches='tight')

fig=plt.figure(figsize=(8,6))
plt.plot((allEigLmv[0]-allEigLmv[1]), label='t=%d' %(1*dt))
plt.legend(loc='upper center', bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=4)
plt.savefig('eigenvalues_diff_0-3.pdf',bbox_inches='tight')
print np.real(allEigLmv[0][len(allEigLmv[0])-10:len(allEigLmv[0])])
