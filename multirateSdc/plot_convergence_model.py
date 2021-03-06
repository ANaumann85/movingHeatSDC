from problem_model import problem_model
from sdc import sdc_step
from sdc_standard import sdc_standard_step
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
err_std = np.zeros((np.size(K_iter),np.size(nsteps)))

order  = np.zeros((np.size(K_iter),np.size(nsteps)))

fs     = 8

a     = 1.0
nu    = 1.1
alpha = 1.0e-2
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
    u0_ros = copy.deepcopy(u0)
    u0_sdc = copy.deepcopy(u0)
    
    for n in range(nsteps[ll]):
      t_n    = float(n)*dt
      t_np1  = float(n+1)*dt
      sdc    = sdc_step(M, P, t_n, t_np1, prob)
      sdc_std = sdc_standard_step(M, t_n, t_np1, prob)
      
      if kk==0:
        ros2   = ros2_step(t_n, t_np1, prob)
        u0_ros  = ros2.step(u0_ros)
      
      ### Multi-rate SDC
      
      # reset buffers to zero
      u     = np.zeros((M,prob.dim))
      usub  = np.zeros((M,P,prob.dim))
      fu     = np.zeros((M,prob.dim))
      fu_sub  = np.zeros((M,P,prob.dim))
      
      sdc.predict(u0, u, usub, fu, fu_sub)
      for k in range(K_iter[kk]):
        sdc.sweep(u0, u, usub, fu, fu_sub)

      u0 = u[M-1]
      
      ### Single-step IMEX SDC
      u   = np.zeros((M,prob.dim))
      u_  = np.zeros((M,prob.dim))
      sdc_std.predict(u0_sdc, u_)
      for k in range(K_iter[kk]):
        sdc_std.sweep(u0_sdc, u, u_)
        u_  = copy.deepcopy(u)
  
      u0_sdc = u[M-1]
      
    ###
    err[kk,ll]   = np.linalg.norm(u0 - u_ex, np.inf)/np.linalg.norm(u_ex, np.inf)
    order[kk,ll] = err[kk,0]*(float(nsteps[0])/float(nsteps[ll]))**order_p
    err_std[kk,ll] = np.linalg.norm(u0_sdc - u_ex, np.inf)/np.linalg.norm(u_ex, np.inf)
    if kk==0:
      err_ros2[ll] = np.linalg.norm(u0_ros - u_ex, np.inf)/np.linalg.norm(u_ex, np.inf)

#
# Convergence plots
#

#rcParams['figure.figsize'] = 2.5, 2.5
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

#plt.gca().get_xaxis().get_major_formatter().labelOnlyBase = False
#plt.gca().get_xaxis().set_major_formatter(ScalarFormatter())
plt.savefig('convergence_model_M_%d_P_%d.pdf' %(M,P))
#plt.show()

#
# Error constants
#
err_c      = np.zeros((np.size(K_iter),np.size(nsteps)))
err_c_ros2 = np.zeros(np.size(nsteps))
err_c_std  = np.zeros((np.size(K_iter),np.size(nsteps)))
for kk in range(np.size(K_iter)):
  order_p = np.min([K_iter[kk], 2*M-1, 2*P-1])

  for ll in range(np.size(nsteps)):
    ns = nsteps[ll]
    dt = tend/float(ns)
    err_c[kk,ll] = err[kk,ll]/(dt**order_p)
    err_c_std[kk,ll] = err_std[kk,ll]/(dt**order_p)
    if kk==0:
      err_c_ros2[ll] = err_ros2[ll]/(dt**2)

fig = plt.figure()
plt.semilogy(nsteps, err_c_ros2, 'kd-', markersize=fs, label="Ros(2)")

plt.semilogy(nsteps, err_c[0,:], 'bo-', markersize=fs, label=("K=%1i" % K_iter[0]))
plt.semilogy(nsteps, err_c_std[0,:], 'b^-', markersize=fs, label=("Std-K=%1i" % K_iter[0]))

plt.semilogy(nsteps, err_c[1,:], 'ro-', markersize=fs, label=("K=%1i" % K_iter[1]))
plt.semilogy(nsteps, err_c_std[1,:], 'r^-', markersize=fs, label=("Std-K=%1i" % K_iter[1]))

plt.semilogy(nsteps, err_c[2,:], 'go-', markersize=fs, label=("K=%1i" % K_iter[2]))
plt.semilogy(nsteps, err_c_std[2,:], 'g^-', markersize=fs, label=("Std-K=%1i" % K_iter[2]))
plt.xlim([0.95*nsteps[0], 1.05*nsteps[-1]])
plt.legend(loc='lower left', fontsize=fs, prop={'size':fs})

#plt.gca().get_xaxis().get_major_formatter().labelOnlyBase = False
#plt.gca().get_xaxis().set_major_formatter(ScalarFormatter())
plt.show()
