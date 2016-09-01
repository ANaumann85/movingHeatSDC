from problem_inh import problem_inh
from sdc import sdc_step
import numpy as np
import unittest
import copy

nu = -1.0
prob = problem_inh(nu)
M = 5
P = 5
tstart = 0.0
tend   = 0.5
nIter=10
u0 = 1.0

dt=(tend-tstart)/nIter
u0_=u0
res=[]
for n in range(nIter):
  ts=tstart+n*dt
  te=ts+dt
  sdc = sdc_step(M, P, ts, te, prob)
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
c1  = u0 + 1.0/(nu**2+1)
uex = c1*np.exp(nu*tend) - (nu*np.sin(tend) + np.cos(tend))/(nu**2+1)
err =  abs(uex - u[-1])
print err
f=open("problem_inh.dat",'w')
print res
for k in res:
  f.write(str(k)+"\n")
f.close()
