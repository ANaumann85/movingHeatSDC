from problem_model import problem_model
from sdc import sdc_step
import numpy as np
import copy

from matplotlib import pyplot as plt
from pylab import rcParams
from matplotlib.ticker import ScalarFormatter
from subprocess import call

def uex(t, y0, a, nu):
  pass

M = 2
P = 5
tstart = 0.0
tend   = 10.0
nsteps = 100
dt = (tend - tstart)/float(nsteps)

a    = 1.0
nu   = -1.0
prob = problem_model(a=a, nu=nu, u0=1.0, v0=0.75)

K_iter = 12

u0    = [1.0, 0.0, 1.0, 0.0]
u_ex  = uex(tend,u0, a, nu)
u_plot = np.zeros((nsteps+1,prob.dim))
u_euler = np.zeros((nsteps+1, prob.dim))
t_axis = np.zeros(nsteps+1)
t_axis[0] = tstart
u_plot[0] = u0
u_euler[0] = u0

for n in range(nsteps):
  tstart = float(n)*dt
  tend   = float(n+1)*dt
  
  u_temp       = u_euler[n] + dt*prob.f2(u_euler[n], tstart)
  u_euler[n+1] = prob.solve_f1(dt, u_temp)
  
  sdc    = sdc_step(M, P, tstart, tend, prob)
  
  # reset buffers to zero
  u_    = np.zeros((M,prob.dim))
  usub_ = np.zeros((M,P,prob.dim))
  u     = np.zeros((M,prob.dim))
  usub  = np.zeros((M,P,prob.dim))

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
#print ("error:             %5.3e " % (np.linalg.norm(u0 - u_ex, np.inf)/np.linalg.norm(u_ex, np.inf)))

#rcParams['figure.figsize'] = 2.5, 2.5
fig = plt.figure()
plt.plot(t_axis, u_plot[:,0], 'b')
plt.plot(t_axis, u_euler[:,0], 'b--')

plt.plot(t_axis, u_plot[:,1], 'r')
plt.plot(t_axis, u_euler[:,1], 'r--')

plt.plot(t_axis, u_plot[:,2], 'g')
plt.plot(t_axis, u_euler[:,2], 'g--')

plt.plot(t_axis, u_plot[:,3], 'k')
plt.plot(t_axis, u_euler[:,3], 'k--')
plt.show()