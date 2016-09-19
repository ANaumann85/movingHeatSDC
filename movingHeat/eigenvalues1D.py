from scipy.io import *
from scipy.sparse import linalg as sLA
from scipy import linalg as dLA
import numpy as np
import matplotlib.pyplot as plt



#periodic 1d heat with source al\delta(x-at)(v0-u(t,x)), i.e. 
#\dot{u} = \Delta u + al\delta(x-at)(v0-u(t,x))
#approximate \delta(x-at) as point or as triangle with unit area
#stencil point:
#\dot u_{0}=(u_{N}-2u_0+u_1)/h^2 +al\delta(x_0-at)(v_0-u_0)
#\dot u_{j}=(u_{j-1}-2u_{j}+u_{j+1})/h^2+al\delta(x_j-at)(v_0-u_j), j=N-1..1
#\dot u_{N}=(u_{N-1}-2u_N+u_1)/h^2+al\delta(x_N-at)(v_0-u_N)
#stencil unit area:
#\dot u_{0}=(u_{N}-2u_0+u_1)/h^2 +al*tri(x_0-at,h)(v_0-u_0)/h
#\dot u_{j}=(u_{j-1}-2u_{j}+u_{j+1})/h^2+al*tri(x_j-at,h)(v_0-u_j)/h, j=N-1..1
#\dot u_{N}=(u_{N-1}-2u_N+u_1)/h^2+al*tri(x_N-at,h)(v_0-u_N)/h
#utilizing tri(x,h)=1-|x/h| 0<=x<h, 0 else

alpha=1.0 #5e-4
nu=1.0 #1e-3
L=1.0
N=10
h=L/N
a=1
t0=0.0
te=1.0
nT=5
#main laplacian
A=np.diag(np.ones(N-1),-1)-2*np.diag(np.ones(N))+np.diag(np.ones(N-1),1)
A[0,N-1]=1
A[N-1,0]=1
A=A*nu

eigLapl, v=dLA.eigh(A)
eigLapl=np.sort(eigLapl)
dt=te/nT
allEigs=[]
for k in range(nT):
    t=k*dt
    jc=int(a*t/h)
    B=np.copy(A)
    B[jc,jc] +=-alpha/h
    l,v=dLA.eig(B)
    allEigs.append(np.real(np.sort(l)))

diffEigs=np.zeros(nT)
for k in range(1,nT):
    diffEigs[k-1]=max(abs(allEigs[0]-allEigs[k]))
print "maximal absolute difference of eigenvalues at t and t=0", diffEigs
#print A
plt.figure(num=1,figsize=(8,8))
plt.plot(eigLapl, '--', label='lapl')
for k in range(nT):
    plt.plot(allEigs[k], label='t=%0.1f' % (k*dt))
#eigenvalues from model problem, N=2
eigsN2=np.sort([-0.07957747152, -0.9920221875, -4.167152758, -1.309421233, -4.168043595])
eigsN2=eigsN2[::-1]
posN2=range(N-1,N-len(eigsN2)-1,-1)

eigsN3=np.sort([-0.07957747152, -0.996893609, -9.167314174, -4.154102106, -1.280565482, -9.167871777, -4.188182517])
eigsN3=eigsN3[::-1]
posN3=range(N-1,N-len(eigsN3)-1,-1)

eigsN4=np.sort([-0.07957747152, -1.272238367, -4.341755173, -9.173960806, -16.16697531, -0.998328048, -3.993141516, -9.160222213, -16.16661811])
eigsN4=eigsN4[::-1]
posN4=range(N-1,N-len(eigsN4)-1,-1)

plt.plot(posN2, eigsN2, 'x', label='N=2')
plt.plot(posN3, eigsN3, 'x', label='N=3')
plt.plot(posN4, eigsN4, 'x', label='N=4')
plt.legend(loc='upper center', bbox_to_anchor=(0.5,-0.0), fancybox=True, shadow=True, ncol=3)
#plt.legend(loc='best')#, ncol=6)
#plt.show()
plt.savefig('eigenvalues1D.pdf',bbox_inches='tight')

#for k in range(nT):
#    print max(np.imag(allEigs[k]))
#    plt.plot(np.imag(allEigs[k]))
#plt.show()
