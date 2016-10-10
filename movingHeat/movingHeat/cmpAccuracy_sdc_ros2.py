#from pylab import rcParams
import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt 
from matplotlib.patches import Polygon
from subprocess import call

import numpy as np


acros2=np.loadtxt('ros2/accuracy.dat')
P_vec=[2,3,6,8]
M_vec=[2,3]
kiter_vec=[1,2]
markers={2:'s', 3:'o'}
colors={2:'k',3:'b',6:'r',8:'g'}
figStep=plt.figure(0)
figTime=plt.figure(1)
for k_iter in kiter_vec:
    for M in M_vec:
        for p in P_vec:
            acsdc=np.loadtxt('sdc/kiter_%d/accuracy_M_%d_P_%d.dat' %(k_iter,M,p))
            plt.figure(0)
            plt.loglog(20.0/acsdc[:,0],acsdc[:,2],label='sdc-M%d_P%d_kiter%d' %(M,p,k_iter),marker=markers[M],color=colors[p])
            plt.figure(1)
            plt.loglog(acsdc[:,4],acsdc[:,2],label='sdc-M%d_P%d_kiter%d' %(M,p,k_iter),marker=markers[M],color=colors[p])
    plt.figure(0)
    plt.loglog(20.0/acros2[:,0],acros2[:,2],label='ros2',color='y')
    plt.xlabel('nStep')
    plt.ylabel('maxErr')
    plt.legend(bbox_to_anchor=(0.5,-0.01),ncol=2)
    plt.figure(1)
    plt.loglog(acros2[:,4],acros2[:,2],label='ros2',color='y')
    plt.xlabel('runtime [ms]')
    plt.ylabel('maxErr')
    plt.legend()
    #plt.show()
    figStep.savefig('cmpAccuracy_sdc_ros2_nstep_kiter_%d.pdf' %(k_iter))
    figTime.savefig('cmpAccuracy_sdc_ros2_time_kiter_%d.pdf' %(k_iter)) 
