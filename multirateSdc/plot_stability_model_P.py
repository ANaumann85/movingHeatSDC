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
nue=50
a = 1.0
alpha_vec = np.linspace(0.0, 35.0, 20)
nu_vec    = np.linspace(0.0, nue, 20)
P_vec = range(M,M+11)

K_iter = 2
fig1  = plt.figure(0)
fig2 =plt.figure(1)
cmap=plt.get_cmap('Accent')
for P in P_vec:
  fname='stability-K-%d-M-%d-P-%d-nue_%d.txt' %(K_iter, M, P, nue)
  stab=np.loadtxt(fname)
  scal=(P-M)/(M+11.0)
  curColor=cmap(scal)
  plt.figure(0)
  CS2 = plt.contour(nu_vec, alpha_vec, np.absolute(stab), [1.0]) #,  color=curColor)
  p=CS2.collections[0].get_paths()[0]
  x=p.vertices[:,0]
  y=p.vertices[:,1]
  plt.figure(1)
  lines=plt.plot(x,y,color=curColor,label='P=%d' %P)


fs = 16
plt.xlabel(r'$\nu$', fontsize=fs)
#plt.ylabel(r'$a$', fontsize=fs)
plt.ylabel(r'$\alpha$', fontsize=fs)
plt.legend(bbox_to_anchor=(-0.1, -0.02, 1., -0.102), loc=2,ncol=len(P_vec)/2)
filename = 'stabilityP-K'+str(K_iter)+'-M'+str(M)+'-nue_'+str(nue)+'.pdf'
fig2.savefig(filename, bbox_inches='tight')
call(["pdfcrop", filename, filename])
#plt.show()

