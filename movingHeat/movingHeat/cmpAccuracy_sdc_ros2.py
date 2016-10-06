#from pylab import rcParams
import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt 
from matplotlib.patches import Polygon
from subprocess import call

import numpy as np


acros2=np.loadtxt('ros2/accuracy.dat')
P_vec=[2,3,6,8]
fig=plt.figure()
k_iter=2
for p in P_vec:
    acsdc=np.loadtxt('sdc/kiter_%d/accuracy_P_%d.dat' %(k_iter,p))
    plt.loglog(20.0/acsdc[:,0],acsdc[:,2],label='sdc-M3_P%d_kiter%d' %(p,2),marker='o')
    plt.loglog(20.0/acsdc[:,0],acsdc[:,2],label='sdc-M3_P%d_kiter%d' %(p,5),marker='v')
    plt.loglog(20.0/acsdc[:,0],acsdc[:,2],label='sdc-M2_P%d_kiter%d' %(p,2),marker='s')
plt.loglog(20.0/acros2[:,0],acros2[:,2],label='ros2')
plt.xlabel('nStep')
plt.ylabel('maxErr')
plt.legend()
#plt.show()
plt.savefig('cmpAccuracy_sdc_ros2.pdf')

