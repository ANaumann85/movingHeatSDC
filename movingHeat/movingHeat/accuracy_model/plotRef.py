import numpy as np
import copy

#from pylab import rcParams
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from subprocess import call

alpha_vec=np.linspace(0.0, 30, 60);
nu_vec=np.linspace(0.0, 100, 200);


fname='relSdcMrsdc.dat'
sdcMrsdc=np.loadtxt(fname)
plt.figure(0)
CS=plt.contour(nu_vec, alpha_vec, sdcMrsdc, [0.1,1.0, 10.0])
plt.clabel(CS)
plt.show()
