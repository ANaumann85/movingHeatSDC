from pylab import rcParams
import matplotlib.pyplot as plt 
from matplotlib.patches import Polygon
from subprocess import call

import numpy as np


acsdc=np.loadtxt('sdc/accuracy.dat')
acros2=np.loadtxt('ros2/accuracy.dat')

fig=plt.figure()
plt.loglog(20.0/acsdc[:,0],acsdc[:,2],label='sdc-M3_P2_kiter5')
plt.loglog(20.0/acros2[:,0],acros2[:,2],label='ros2')
plt.xlabel('nStep')
plt.ylabel('maxErr')
plt.legend(loc=2)
#plt.show()
plt.savefig('cmpAccuracy_sdc_ros2.pdf')

