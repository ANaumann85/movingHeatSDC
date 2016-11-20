import numpy as np
import copy

#from pylab import rcParams
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from subprocess import call

alpha_vec=np.linspace(0.0, 30, 60);
nu_vec=np.linspace(0.0, 100, 200);


Mv=[2,3]
Pv=[2,4,6,8]
iterV=[1,2,3,4]
nStepV=[1,2,4,8]
fnameFo='M_%(M)d_P_%(P)d_iter_%(iter)d_step_%(step)d/relSdcMrsdc.dat'
titFo="M=%(M)d, P=%(P)d, iter=%(iter)d, nStep=%(step)d"
levels=[0.1,1.0, 10.0]

#M=Mv[0]
#P=Pv[0]
#for nStep in nStepV:
#    for iter in iterV:
#        data={'M':M, 'P':P, 'iter':iter, 'step':nStep}
#        fname=fnameFo % data
#        sdcMrsdc=np.loadtxt(fname)
#        plt.figure(0)
#        CS=plt.contour(nu_vec, alpha_vec, sdcMrsdc, levels)
#        for p in range(len(levels)):
#            CS.collections[p].set_label(str(levels[p]))
#        plt.legend()
#        plt.xlim((0, 10))
#        plt.title(titFo % data)
#        plt.show()

iter=iterV[2]
nStep=nStepV[2]

for M in Mv:
    for P in Pv:
        data={'M':M, 'P':P, 'iter':iter, 'step':nStep}
        fname=fnameFo % data
        sdcMrsdc=np.loadtxt(fname)
        plt.figure(0)
        CS=plt.contour(nu_vec, alpha_vec, sdcMrsdc, levels)
        for p in range(len(levels)):
            CS.collections[p].set_label(str(levels[p]))
        plt.legend()
        plt.xlim((0, 10))
        plt.title(titFo % data)
        plt.show()
