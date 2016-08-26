from problem_model import problem_model
from sdc import sdc_step
import numpy as np
from scipy.integrate import odeint
import copy

from matplotlib import pyplot as plt
from pylab import rcParams
from matplotlib.ticker import ScalarFormatter
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

order  = np.zeros((np.size(K_iter),np.size(nsteps)))

fs     = 8

a     = 1.0
nu    = 1.1
alpha = 25.0
v0    = 0.25
prob  = problem_model(a, nu, alpha, v0)

u_ex = np.zeros(prob.dim)
u0    = np.zeros(prob.dim)
u0[0] = 1.0
sol = odeint(prob.f, u0, [0.0, tend], rtol=1e-12, atol=1e-12)
u_ex = sol[-1]

for kk in range(np.size(K_iter)):

  order_p = np.min([K_iter[kk], 2*M-1, 2*P-1])

  for ll in range(np.size(nsteps)):

    dt    = (tend - tstart)/float(nsteps[ll])
    u0    = np.zeros(prob.dim)
    u0[0] = 1.0
    u_ros = u0
    
    for n in range(nsteps[ll]):
      t_n    = float(n)*dt
      t_np1  = float(n+1)*dt
      sdc    = sdc_step(M, P, t_n, t_np1, prob)
      
      if kk==0:
        ros2   = ros2_step(t_n, t_np1, prob)
        u_ros  = ros2.step(u_ros)
      
      # reset buffers to zero
      u     = np.zeros((M,prob.dim))
      usub  = np.zeros((M,P,prob.dim))
      u_    = np.zeros((M,prob.dim))
      usub_ = np.zeros((M,P,prob.dim))

      sdc.predict(u0, u_, usub_)
      for k in range(K_iter[kk]):
        sdc.sweep(u0, u, usub, u_, usub_)
        u_    = copy.deepcopy(u)
        usub_ = copy.deepcopy(usub)

      u0 = u[M-1]
      
    ###
    err[kk,ll]   = np.linalg.norm(u0 - u_ex, np.inf)/np.linalg.norm(u_ex, np.inf)
    order[kk,ll] = err[kk,0]*(float(nsteps[0])/float(nsteps[ll]))**order_p
    if kk==0:
      err_ros2[ll] = np.linalg.norm(u_ros - u_ex, np.inf)/np.linalg.norm(u_ex, np.inf)

#rcParams['figure.figsize'] = 2.5, 2.5
fig = plt.figure()
plt.loglog(nsteps, err_ros2, 'k^', markersize=fs, label="Ros(2)")

plt.loglog(nsteps, err[0,:], 'bo', markersize=fs, label=("K=%1i" % K_iter[0]))
plt.loglog(nsteps, order[0,:], '-', color='b')
plt.loglog(nsteps, err[1,:], 'rd', markersize=fs, label=("K=%1i" % K_iter[1]))
plt.loglog(nsteps, order[1,:], '-', color='r')
plt.loglog(nsteps, err[2,:], 'gs', markersize=fs, label=("K=%1i" % K_iter[2]))
plt.loglog(nsteps, order[2,:], '-', color='g')
plt.xlim([0.95*nsteps[0], 1.05*nsteps[-1]])
plt.legend(loc='lower left', fontsize=fs, prop={'size':fs})

#plt.gca().get_xaxis().get_major_formatter().labelOnlyBase = False
#plt.gca().get_xaxis().set_major_formatter(ScalarFormatter())
plt.show()
