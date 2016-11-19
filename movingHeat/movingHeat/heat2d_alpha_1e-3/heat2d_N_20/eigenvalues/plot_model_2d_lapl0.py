import numpy as np
import copy

#from pylab import rcParams
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from subprocess import call

#from http://stackoverflow.com/questions/14733605/python-matplotlib-label-title-wrong-character
#plt.matplotlib.rc('text', usetex=True)

M = 3
P = 8
tstart = 0.0
tend   = 1.0
nue=500
alphae=60
alpha_vec = np.linspace(0.0, alphae, 120)
nu_vec    = np.linspace(0.0, nue, 50)
a = 1.0
K_iter = 2
theta=1.0

te=20
dt=[te, te/2]
for x in range(6):
    dt.append(te/(5*2.0**x))
dt=np.array(dt)
print dt 
fs = 16
#shade the area for the 1d-plot
#fname='../../stability-K-%d-M-%d-P-%d-nue_%d_alphae_5_theta-%g.txt' %(K_iter, M, P, nue, theta)
fname='../../stability-K-%d-M-%d-P-%d-nue_%d_theta-%g.txt' %(K_iter, M, P, nue, theta)
stab=np.loadtxt(fname)
plt.figure(0)
CS2 = plt.contour(nu_vec, alpha_vec, np.absolute(stab), [1.0],  colors='k')
#plt.clabel(CS2, fmt='theta=%g'%theta, manual=[(200,15)])
x2d0=10.98/4*dt #
y2d0=0.079/0.4*dt #alpha=50, lmbdaCouple=1e-4, dt=1
print "p2d0:", x2d0,y2d0
plt.plot(x2d0,y2d0, marker='o',color='r', label='2d')

#x2d0L0=2.8/4*dt #nu=1e-3, lmbda=1.0e-3, dt=1
#y2d0L0=1.7/0.4*dt #alpha=1e-3, lmbdaFast=35, dt=1
#print "p2d0L0:", x2d0L0,y2d0L0
#plt.plot(x2d0L0,y2d0L0, marker='o',color='g', label=r'2d,$\lambda_0$')

plt.title('stability regions of model problem and 2d scaled to model')
plt.xlabel(r'$\nu$', fontsize=fs)
#plt.ylabel(r'$a$', fontsize=fs)
plt.ylabel(r'$\alpha$', fontsize=fs)
plt.legend(ncol=2, loc=0)#bbox_to_anchor=(-0.0, -0.10, 1., -0.015),
plt.ylim((-0.1,10))
plt.xlim((-0.1,100.0))
#plt.show()
filename='stability_model_2dPoints-K-%d-M-%d-P-%d.pdf' %(K_iter, M, P)
plt.savefig(filename)
