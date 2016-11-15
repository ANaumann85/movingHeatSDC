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
alpha_vec = np.linspace(0.0, 30.0, 120)
nu_vec    = np.linspace(0.0, nue, 50)
a = 1.0
K_iter = 2
thetaVec=[0, 0.25, 0.5, 1.0]

fs = 16
for n in range(len(thetaVec)):
    theta=thetaVec[n]
    #shade the area for the 1d-plot
    fname='stability-K-%d-M-%d-P-%d-nue_%d_theta-%g.txt' %(K_iter, M, P, nue, theta)
    stab=np.loadtxt(fname)
    plt.figure(0)
    CS2 = plt.contour(nu_vec, alpha_vec, np.absolute(stab), [1.0],  colors='k')
    plt.clabel(CS2, fmt='theta=%g'%theta, manual=[(200,15)])
#p=CS2.collections[0].get_paths()[1]
#x=p.vertices[:,0]
#y=p.vertices[:,1]
#plt.figure(0)
#lines=plt.plot(x,y,color='r',label='model,P=%d' %P)
#fname2d='../movingHeat/movingHeat/stability_sdc/stabEigs_sdc.mtx'
#stab2d=np.loadtxt(fname2d)
#alpha_vec_2d = stab2d[1:,0]
#nu_vec_2d    = stab2d[0,1:]
#stab2d=np.copy(stab2d[1:,1:])
#plt.figure(1)
#cs2d=plt.contour(nu_vec_2d, alpha_vec_2d, np.absolute(stab2d), [1.0])
##plt.close()
#p2d=cs2d.collections[0].get_paths()[0]
#x2d=p2d.vertices[:,0]
#y2d=p2d.vertices[:,1]
##lines=plt.fill_between(nuv2,y,alv2,where=y>alv2,color='b',label='1d,P=%d' %P)
#plt.figure(0)
#lines=plt.plot(x2d*600,y2d*35/0.4, color='b', label='2d,scaled,P=2')
x3d0=50*10/4 #nu=50, lmbdaStaender=10, dt=1
y3d0=50*1e-4/0.4 #alpha=50, lmbdaCouple=1e-4, dt=1
x3d1=50*26/4 #nu=50, lmbdaStock=26, dt=1
y3d1=50*4e-4/0.4 #alpha=50, lmbdaCouple=1e-4, dt=1
print "p3d0:", x3d0,y3d0
print "p3d1:", x3d1,y3d1
plt.plot(x3d0,y3d0, marker='o',color='r', label='3d,staender,dt=1')
plt.plot(x3d1,y3d1, marker='x',color='r', label='3d,stock,dt=1')

x3d0L0=50*10/4 #nu=50, lmbdaStaender=10, dt=1
y3d0L0=50*0.15/0.4 #alpha=50, lmbdaFast=0.15, dt=1
x3d1L0=50*26/4 #nu=50, lmbdaStock=26, dt=1
y3d1L0=50*0.08/0.4 #alpha=50, lmbdaFast=0.08, dt=1
print "p3d0L0:", x3d0L0,y3d0L0
print "p3d1L0:", x3d1L0,y3d1L0
plt.plot(x3d0L0,y3d0L0, marker='o',color='g', label='3d,staender,dt=1')
plt.plot(x3d1L0,y3d1L0, marker='x',color='g', label='3d,stock,dt=1')

plt.title('stability regions of model problem and 3d scaled to model')
plt.xlabel(r'$\nu$', fontsize=fs)
#plt.ylabel(r'$a$', fontsize=fs)
plt.ylabel(r'$\alpha$', fontsize=fs)
#plt.legend(ncol=2, loc=2)#bbox_to_anchor=(-0.0, -0.10, 1., -0.015),
plt.ylim((-1.0,30))
#plt.show()
filename='stability_model_3dPoints-K-%d-M-%d-P-%d.pdf' %(K_iter, M, P)
plt.savefig(filename)
