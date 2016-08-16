from problem_model import problem_model
from sdc import sdc_step
import numpy as np
import copy

from matplotlib import pyplot as plt
from pylab import rcParams
from matplotlib.ticker import ScalarFormatter
from subprocess import call

def uex(t, y0, a, nu):
  p1 = np.exp(nu*t)
  p2 = np.exp((nu+2*a)*t)
  y = np.zeros(2)
  y[0] = 0.5*(y0[0] + y0[1])*p1 + 0.5*(y0[0] - y0[1])*p2
  y[1] = 0.5*(y0[0] + y0[1])*p1 + 0.5*(y0[1] - y0[0])*p2
  return y

M = 3
P = 3
tstart = 0.0
tend   = 2.25
nsteps = [4, 6, 8, 10, 12, 14, 16, 18, 20]
K_iter = [2, 3, 4]
err    = np.zeros((np.size(K_iter),np.size(nsteps)))
order  = np.zeros((np.size(K_iter),np.size(nsteps)))

fs     = 8

a  = -0.25
nu = -1.1
prob     = problem_model(a, nu)

for kk in range(np.size(K_iter)):

  order_p = np.min([K_iter[kk], 2*M-1, 2*P-1])

  for ll in range(np.size(nsteps)):

    dt    = (tend - tstart)/float(nsteps[ll])
    u0    = [1.0, 0.0]
    u_ex  = uex(tend, u0, a, nu)

    for n in range(nsteps[ll]):
      t_n    = float(n)*dt
      t_np1  = float(n+1)*dt
      sdc    = sdc_step(M, P, t_n, t_np1, prob)

      # reset buffers to zero
      u     = np.zeros((M,2))
      usub  = np.zeros((M,P,2))
      u_    = np.zeros((M,2))
      usub_ = np.zeros((M,P,2))

      sdc.predict(u0, u_, usub_)
      for k in range(K_iter[kk]):
        sdc.sweep(u0, u, usub, u_, usub_)
        u_    = copy.deepcopy(u)
        usub_ = copy.deepcopy(usub)

      u0 = u[M-1]
      
    ###
    err[kk,ll]   = np.linalg.norm(u0 - u_ex, np.inf)/np.linalg.norm(u_ex, np.inf)
    order[kk,ll] = err[kk,0]*(float(nsteps[0])/float(nsteps[ll]))**order_p

#rcParams['figure.figsize'] = 2.5, 2.5
fig = plt.figure()
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
