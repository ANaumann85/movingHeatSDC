from problem_model import problem_model
from sdc import sdc_step
from sdc_standard import sdc_standard_step
import numpy as np
from scipy.integrate import odeint
import copy

from matplotlib import pyplot as plt
from pylab import rcParams
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator
from subprocess import call
from ros2 import ros2_step

M = 3
P = 5
tstart = 0.0
tend   = 2.5
nsteps = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70]
K_iter = [2, 3, 4]
err    = np.zeros((np.size(K_iter),np.size(nsteps)))
err_ros2 = np.zeros(np.size(nsteps))
err_std = np.zeros((np.size(K_iter),np.size(nsteps)))

order  = np.zeros((np.size(K_iter),np.size(nsteps)))

fs     = 8

a     = 1.0
nu    = 1.1
alpha = 25.0
v0    = 0.25

nsteps=np.loadtxt('convergence_nsteps.txt')
err_ros2=np.loadtxt('convergence_err_ros2.txt')
err=np.loadtxt('convergence_err.txt')
err_std=np.loadtxt('convergence_err_std.txt')
for kk in range(np.size(K_iter)):
  order_p = np.min([K_iter[kk], 2*M-1, 2*P-1])
  for ll in range(np.size(nsteps)):
    order[kk,ll] = err[kk,0]*(float(nsteps[0])/float(nsteps[ll]))**order_p

fig = plt.figure()
plt.loglog(nsteps, err_ros2, 'kd', markersize=fs, label="Ros(2)")

plt.loglog(nsteps, err[0,:], 'bo', markersize=fs, label=("K=%1i" % K_iter[0]))
plt.loglog(nsteps, err_std[0,:], 'b^', markersize=fs, label=("Std-K=%1i" % K_iter[0]))
plt.loglog(nsteps, order[0,:], '-', color='b')

plt.loglog(nsteps, err[1,:], 'ro', markersize=fs, label=("K=%1i" % K_iter[1]))
plt.loglog(nsteps, err_std[1,:], 'r^', markersize=fs, label=("Std-K=%1i" % K_iter[1]))
plt.loglog(nsteps, order[1,:], '-', color='r')

plt.loglog(nsteps, err[2,:], 'go', markersize=fs, label=("K=%1i" % K_iter[2]))
plt.loglog(nsteps, err_std[2,:], 'g^', markersize=fs, label=("Std-K=%1i" % K_iter[2]))
plt.loglog(nsteps, order[2,:], '-', color='g')
plt.xlim([0.95*nsteps[0], 1.05*nsteps[-1]])

plt.legend(loc='lower left', fontsize=fs, prop={'size':fs})

plt.xlabel('nstep')
plt.ylabel('maxErr')
#plt.gca().get_xaxis().get_major_formatter().labelOnlyBase = False
plt.gca().get_xaxis().set_major_formatter(ScalarFormatter())
plt.gca().get_xaxis().set_major_locator(LogLocator(subs=[1.0,2.0,4.0]))
plt.gca().get_xaxis().reset_ticks()
plt.show()
#plt.savefig('convergence_model_M_%d_P_%d_withCorrection.pdf' %(M,P))
