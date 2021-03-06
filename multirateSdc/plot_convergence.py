from problem import problem
from sdc import sdc_step
from ros2 import ros2_step
import numpy as np
import copy

from matplotlib import pyplot as plt
from pylab import rcParams
from matplotlib.ticker import ScalarFormatter
from subprocess import call

M = 3
P = 3
tstart = 0.0
tend   = 2.25
nsteps = [4, 6, 8, 10, 12, 14, 16, 18, 20]
K_iter = [1, 2, 3]
err    = np.zeros((np.size(K_iter),np.size(nsteps)))
order  = np.zeros((np.size(K_iter),np.size(nsteps)))
err_ros = np.zeros(np.size(nsteps))

fs     = 8

lambda_1 = -0.1
lambda_2 = -1.0
prob     = problem(lambda_1, lambda_2)

for kk in range(np.size(K_iter)):

  order_p = np.min([K_iter[kk], 2*M-1, 2*P-1])

  for ll in range(np.size(nsteps)):

    dt      = (tend - tstart)/float(nsteps[ll])
    u0      = 2.0
    u0_ros  = u0
    u_ex    = u0*np.exp(tend*(lambda_1+lambda_2))

    for n in range(nsteps[ll]):
      t_n    = float(n)*dt
      t_np1  = float(n+1)*dt
      sdc    = sdc_step(M, P, t_n, t_np1, prob)
      
      if kk==0:
        ros    = ros2_step(t_n, t_np1, prob)
        u0_ros = ros.step(u0_ros)
    
      # reset buffers to zero
      u      = np.zeros((M,1,1))
      usub   = np.zeros((M,P,1))
      fu     = np.zeros((M,1,1))
      fu_sub = np.zeros((M,P,1))

      sdc.predict(u0, u, usub, fu, fu_sub)
      for k in range(K_iter[kk]):
        sdc.sweep(u0, u, usub, fu, fu_sub)

      u0 = u[M-1]

    ###
    if kk==0:
      err_ros[ll]  = abs(u0_ros - u_ex)/abs(u_ex)
    err[kk,ll]   = abs(u0 - u_ex)/abs(u_ex)
    order[kk,ll] = err[kk,0]*(float(nsteps[0])/float(nsteps[ll]))**order_p

#rcParams['figure.figsize'] = 2.5, 2.5
fig = plt.figure()
plt.loglog(nsteps, err[0,:], 'bo', markersize=fs)
plt.loglog(nsteps, err_ros, 'bx', markersize=fs)
plt.loglog(nsteps, order[0,:], '-', color='b')

plt.loglog(nsteps, err[1,:], 'ro', markersize=fs)
plt.loglog(nsteps, order[1,:], '-', color='r')

plt.loglog(nsteps, err[2,:], 'go', markersize=fs)
plt.loglog(nsteps, order[2,:], '-', color='g')
plt.xlim([0.95*nsteps[0], 1.05*nsteps[-1]])
plt.ylabel('Error')
plt.xlabel('Number of time steps')
#plt.gca().get_xaxis().get_major_formatter().labelOnlyBase = False
#plt.gca().get_xaxis().set_major_formatter(ScalarFormatter())

filename = 'convergence.pdf'
fig.savefig(filename, bbox_inches='tight')
call(["pdfcrop", filename, filename])


fig = plt.figure()
plt.loglog(np.multiply(nsteps, M*K_iter[0]), err[0,:], 'bo', markersize=fs)
plt.loglog(np.multiply(nsteps,2), err_ros, 'bx', markersize=fs)
plt.loglog(np.multiply(nsteps, M*K_iter[1]), err[1,:], 'ro', markersize=fs)
plt.loglog(np.multiply(nsteps, M*K_iter[2]), err[2,:], 'go', markersize=fs)
plt.ylabel('Error')
plt.xlabel('Number of implicit solves')
plt.show()
