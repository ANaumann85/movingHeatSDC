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

M = 5
P = 4
tstart = 0.0
tend   = 10.0
nsteps = 25
dt = (tend - tstart)/float(nsteps)

a    = 0.25
nu   = -1.0
prob = problem_model(a, nu)

K_iter = 12
u_    = np.zeros((M,2))
usub_ = np.zeros((M,P,2))
u     = np.zeros((M,2))
usub  = np.zeros((M,P,2))
u0    = [2.0, 1.0]
u_ex  = uex(tend,u0, a, nu)
u_plot = np.zeros((nsteps+1,2))
t_axis = np.zeros(nsteps+1)
t_axis[0] = tstart
u_plot[0] = u0

for n in range(nsteps):
  tstart = float(n)*dt
  tend   = float(n+1)*dt
  sdc    = sdc_step(M, P, tstart, tend, prob)
  
  # reset buffers to zero
  u_    = np.zeros((M,2))
  usub_ = np.zeros((M,P,2))
  u     = np.zeros((M,2))
  usub  = np.zeros((M,P,2))

  sdc.predict(u0, u_, usub_)
  for k in range(K_iter):
    sdc.sweep(u0, u, usub, u_, usub_)
    u_    = copy.deepcopy(u)
    usub_ = copy.deepcopy(usub)
  print ("+++ step: %3i +++ " % n)
  print ("standard residual: %5.3e " % sdc.residual(u0, u))
  print ("embedded residual: %5.3e " % sdc.sub_residual(u0, usub))
  u0 = u[M-1]
  u_plot[n+1,:] = u0
  t_axis[n+1]   = tend

###
print ("error:             %5.3e " % (np.linalg.norm(u0 - u_ex, np.inf)/np.linalg.norm(u_ex, np.inf)))

#rcParams['figure.figsize'] = 2.5, 2.5
fig = plt.figure()
plt.plot(t_axis, u_plot[:,0], 'b')
plt.plot(t_axis, u_plot[:,1], 'r')
plt.show()