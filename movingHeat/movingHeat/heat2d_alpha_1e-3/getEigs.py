from scipy.io import *
from scipy.sparse import linalg as sLA
from scipy import linalg as dLA
import numpy as np
import matplotlib.pyplot as plt
lapl=mmread('laplace.mm') 
mass=mmread('mass.mm')
nEig=6
eigL=sLA.eigsh(lapl, nEig, mass, return_eigenvectors=False, tol=1.0e-8)

print "eigsLapl:", eigL
mvName='mv-%0.2f.mm' % 0
mv=mmread(mvName)
eigMv=sLA.eigsh(mv, nEig, mass, return_eigenvectors=False, tol=1.0e-8)
print "eigs Coupling:", eigMv

#nT=1
#te=15
#dt=te/nT
#allEigLmv=[]
#allEigMv=[]
#for k in range(0, nT):
#    t=k*dt
#    lmvName='lmv-%0.2f.mm' % t
#    laplMv=mmread(lmvName)
#    [eigLmv, vL]=sLA.eigsh(laplMv, nEig, mass)
#    #[eigLmv, vL]=dLA.eig(laplMv.todense(),mass.todense()) #,right=False)
#    allEigLmv.append(np.sort(eigLmv))
#    mv=mmread(mvName)
#    [eigLmv, vL]=sLA.eigsh(laplMv, nEig, mass)
#    #[eigMv, vL]=dLA.eig(mv.todense(),mass.todense()) #,right=False)
#    allEigMv.append(np.sort(eigMv))

#print 'difference of eigenvalues at time t to eigenvalues at time t=0'
#for k in range(1, nT):
#    print "diff Lapl:",max(abs(np.divide(allEigLmv[0]-allEigLmv[k],allEigLmv[0])))
#    print "diff Mv:",max(abs(allEigMv[0]-allEigMv[k]))
#x=range(0,nEig)
#fig=plt.figure(figsize=(8,6))
#plt.plot(eigL,'k-', label='lapl')
#for k in range(0, nT):
#    plt.plot(np.real(allEigLmv[k]),'-', label='t=%d' % (k*dt))
##plt.show()
#ax1=plt.gca()
#box = ax1.get_position()
#ax1.set_position([box.x0, box.y0 + box.height * 0.1,
#        box.width, box.height * 0.9])
#plt.legend(loc='upper center', bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=3)
#plt.savefig('eigenvalues.pdf',bbox_inches='tight')
#
#fig=plt.figure(figsize=(8,6))
#for k in range(1, nT):
#    plt.plot(abs(np.divide(allEigLmv[0]-allEigLmv[k],allEigLmv[0])), label='t=%d' %(k*dt))
#plt.legend(loc='upper center', bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=4)
#plt.savefig('eigenvalues_reldiff0.pdf',bbox_inches='tight')
#
#fig=plt.figure(figsize=(8,6))
#for k in range(1, nT):
#    plt.plot(abs(allEigLmv[0]-allEigLmv[k]), label='t=%d' %(k*dt))
#plt.legend(loc='upper center', bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=4)
#plt.savefig('eigenvalues_absdiff0.pdf',bbox_inches='tight')
#
#fig=plt.figure(figsize=(8,6))
#plt.plot(abs(allEigLmv[0]-allEigLmv[1]), label='t=%d' %(1*dt))
#plt.legend(loc='upper center', bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=4)
#plt.savefig('eigenvalues_absdiff_0-3.pdf',bbox_inches='tight')
#
#fig=plt.figure(figsize=(8,6))
#plt.plot((allEigLmv[0]-allEigLmv[1]), label='t=%d' %(1*dt))
#plt.legend(loc='upper center', bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=4)
#plt.savefig('eigenvalues_diff_0-3.pdf',bbox_inches='tight')
#print np.real(allEigLmv[0][len(allEigLmv[0])-10:len(allEigLmv[0])])
#
#fig=plt.figure()
#plt.title(r"eigenvalues of moving part at $\alpha=1$")
#for k in range(0, nT):
#    plt.plot(np.real(allEigMv[k]),'-', label='t=%d' % (k*dt))
#plt.legend(loc='upper center', bbox_to_anchor=(0.5,-0.05), fancybox=True, shadow=True, ncol=4)
#plt.savefig('eigenvalues_mv.pdf',bbox_inches='tight')
#print "max eigenvalues lapl mv:",max(abs(eigL)), max(abs(allEigMv[0]))
