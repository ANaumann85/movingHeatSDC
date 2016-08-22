from problem_model import problem_model
from problem_full import problem_full
from sdc import sdc_step
import numpy as np
import copy

from matplotlib import pyplot as plt
from pylab import rcParams
from matplotlib.ticker import ScalarFormatter
from subprocess import call

def delta(a, t, x):
  u = 1.0/(2.0*np.pi) + 0.0*x
  for k in range(1,3):
    u+= 1.0/np.pi*( np.cos(float(k)*a*t)*np.cos(float(k)*x) + np.sin(float(k)*a*t)*np.sin(float(k)*x) )
  return u


M      = 2
P      = 5
tstart = 0.0
tend   = 25.0
nsteps = 50
dt     = (tend - tstart)/float(nsteps)

Nx = 50
xaxis = np.linspace(0, 2*np.pi, Nx+1)
xaxis = xaxis[0:Nx]

prob = problem_model(a=1.0, nu=-1.0, alpha=2.0, u0=0.75, v0=1.0)
K_iter = 12

u0         = [2.0*prob.u0, 1.0, 0.0, 1.0, 0.0]
u_plot     = np.zeros((nsteps+1,prob.dim))
u_euler    = np.zeros((nsteps+1, prob.dim))
t_axis     = np.zeros(nsteps+1)
t_axis[0]  = tstart
u_plot[0]  = u0
u_euler[0] = u0

fig = plt.figure()

for n in range(nsteps):
  tstart = float(n)*dt
  tend   = float(n+1)*dt
  
  u_temp       = u_euler[n] + dt*prob.f2(u_euler[n], tstart)
  u_euler[n+1] = prob.solve_f1(dt, u_temp)
  
  #
  # Solve reduced model first
  #
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

  uplot   = prob.get_solution(u0, xaxis)
  usource = delta(prob.a, tend, xaxis)
  fig.clear()
  plt.plot(xaxis, usource, 'k--')
  plt.plot(xaxis, uplot, 'b')
  plt.ylim([-1.0, 5.0])
  plt.title(("t = %3.2f" % tend))
  plt.show(block=False)
  plt.pause(0.04)

  u0 = u[M-1]
  u_plot[n+1,:] = u0
  t_axis[n+1]   = tend

###
#print ("error:             %5.3e " % (np.linalg.norm(u0 - u_ex, np.inf)/np.linalg.norm(u_ex, np.inf)))

#rcParams['figure.figsize'] = 2.5, 2.5
#fig = plt.figure()
#plt.plot(t_axis, u_plot[:,0], 'k')
#plt.plot(t_axis, u_euler[:,0], 'k--')

#plt.plot(t_axis, u_plot[:,1], 'b')
#plt.plot(t_axis, u_euler[:,1], 'b--')

#plt.plot(t_axis, u_plot[:,2], 'r')
#plt.plot(t_axis, u_euler[:,2], 'r--')

#plt.plot(t_axis, u_plot[:,3], 'g')
#plt.plot(t_axis, u_euler[:,3], 'g--')

#plt.plot(t_axis, u_plot[:,4], 'k')
#plt.plot(t_axis, u_euler[:,4], 'k--')
#plt.show()