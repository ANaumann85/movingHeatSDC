from problem_tcoeff import problem_tcoeff
from sdc import sdc_step
from sdc_standard import sdc_standard_step

import numpy as np
import unittest
import copy

from matplotlib import pyplot as plt

# Consider the IVP
# y'(t) = cos(t)*y(t) - nu*y(t); y(0) = c
# with solution
# y(t) = c*exp[sin(t) - nu*t]

nu     = -0.1
a      = 1.0
prob   = problem_tcoeff(nu, a)
M      = 5
P      = 5
tstart = 0.0
tend   = 10.0
nsteps = 10
u0     = 1.0

dt=(tend-tstart)/float(nsteps)
u0_=u0
u0_std = u0
u_mrsdc=[u0]
u_sdc  = [u0]
for n in range(nsteps):
  ts=tstart+n*dt
  te=ts+dt
  sdc = sdc_step(M, P, ts, te, prob)
  sdc_standard = sdc_standard_step(M, ts, te, prob)
  
  ### multirate SDC ###
  u = np.zeros((M,prob.dim))
  usub = np.zeros((M,P,prob.dim))
  fu = np.zeros((M,prob.dim))
  fu_sub = np.zeros((M,P,prob.dim))

# run predictor
  sdc.predict(u0_, u, usub, fu, fu_sub)

  for k in range(15):
    # run standard node sweep...
    u_old = copy.deepcopy(u)
    usub_old = copy.deepcopy(usub)
    sdc.sweep(u0_, u, usub, fu, fu_sub)
  u0_ = u[-1]
  u_mrsdc.append(u0_)

  ### standard SDC ###
  u_std = np.zeros((M, prob.dim))
  u_old_std = np.zeros((M, prob.dim))

  sdc_standard.predict(u0_std, u_std)
  for k in range(15):
    u_old_std = copy.deepcopy(u_std)
    sdc_standard.sweep(u0_std, u_std, u_old_std)
  u0_std = copy.deepcopy(u_std[-1])
  u_sdc.append(u0_std)

uex = u0*np.exp(a*np.sin(tend) + nu*tend)
err =  abs(uex - u[-1])
print ("MR-SDC error: %5.3e" % err)
print ""
err =  abs(uex - u_std[-1])
print ("Standard SDC error: %5.3e" % err)

taxis = np.linspace(0, tend, nsteps+1)
uex_pl = u0*np.exp(a*np.sin(taxis) + nu*taxis)
plt.figure(1)
plt.plot(taxis, u_mrsdc, 'b')
plt.plot(taxis, u_sdc, 'r')
plt.plot(taxis, uex_pl, 'k')
plt.show()
