import numpy as np
import copy

#from pylab import rcParams
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from subprocess import call

#from http://stackoverflow.com/questions/14733605/python-matplotlib-label-title-wrong-character
plt.matplotlib.rc('text', usetex=True)

M = 3
P = 12
tstart = 0.0
tend   = 1.0
nue=500
alpha_vec = np.linspace(0.0, 35.0, 20)
nu_vec    = np.linspace(0.0, nue, 40)
a = 1.0
K_iter = 2

fs = 16
#shade the area for the 1d-plot
P=12
fname='stability-K-%d-M-%d-P-%d-nue_%d.txt' %(K_iter, M, P, nue)
stab=np.loadtxt(fname)
plt.figure(1)
CS2 = plt.contour(nu_vec, alpha_vec, np.absolute(stab), [1.0]) #,  color=curColor)
p=CS2.collections[0].get_paths()[0]
x=p.vertices[:,0]
y=p.vertices[:,1]
plt.figure(0)
lines=plt.plot(x,y,color='r',label='model,P=%d' %P)
fname1d='stability-1d-K-2-M-3-P-12.txt'
stab1d=np.loadtxt(fname1d)
alpha_vec_1d = np.linspace(0.0, 1.0, 20)
nu_vec_1d    = np.linspace(0.0, 10.0, 20)
plt.figure(1)
cs1d=plt.contour(nu_vec_1d, alpha_vec_1d, np.absolute(stab1d), [1.0])
#plt.close()
p1d=cs1d.collections[0].get_paths()[3]
x1d=p1d.vertices[:,0]
y1d=p1d.vertices[:,1]
#lines=plt.fill_between(nuv2,y,alv2,where=y>alv2,color='b',label='1d,P=%d' %P)
plt.figure(0)
lines=plt.plot(x1d*100,y1d*25, color='b', label='1d,scaled,P=%d' %P)
plt.title('stability regions of model problem and 1d scaled to model')
plt.xlabel(r'$\nu$', fontsize=fs)
#plt.ylabel(r'$a$', fontsize=fs)
plt.ylabel(r'$\alpha$', fontsize=fs)
plt.legend(bbox_to_anchor=(-0.0, -0.10, 1., -0.015),ncol=2, loc=2)
plt.show()
