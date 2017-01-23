from problem_tcoeff import problem_tcoeff
from sdc import sdc_step
from sdc_standard import sdc_standard_step

import numpy as np
import unittest
import copy

nu = -1.0
prob = problem_tcoeff(nu)
M = 5
P = 5
tstart = 0.0
tend   = 0.5
nIter=10
u0 = 1.0

dt=(tend-tstart)/nIter
u0_=u0
u0_std = u0
res=[]
for n in range(nIter):
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
  res.append(u0_[0])

  ### standard SDC ###
  u_std = np.zeros((M, prob.dim))
  u_old_std = np.zeros((M, prob.dim))

  sdc_standard.predict(u0_std, u_std)
  for k in range(15):
    u_old_std = copy.deepcopy(u_std)
    sdc_standard.sweep(u0_std, u_std, u_old_std)
  u0_std = copy.deepcopy(u_std[-1])

uex = u0*(np.exp(np.sin(tend) - nu*tend))
err =  abs(uex - u[-1])
print err
f=open("problem_inh.dat",'w')
print res
for k in res:
  f.write(str(k)+"\n")
f.close()
print ""
print "Standard SDC"
print ""
err =  abs(uex - u_std[-1])
print err
f=open("problem_inh_std.dat",'w')
print res
for k in res:
  f.write(str(k)+"\n")
