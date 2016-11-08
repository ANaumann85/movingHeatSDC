from problem_inh_const import problem_inh_const as problem_inh
from sdc import sdc_step
from sdc_standard import sdc_standard_step
import numpy as np
import unittest
import copy

nu = -1.0
b=2.0
prob = problem_inh(nu, b)
M = 3
P = 8
tstart = 0.0
tend   = 1.0
nIter=2
nKVec=range(1,5)
u0 = -b/nu

dt=(tend-tstart)/nIter
for nK in nKVec:
    u0_=u0
    res=[]
    for n in range(nIter):
      ts=tstart+n*dt
      te=ts+dt
      sdc = sdc_step(M, P, ts, te, prob)
      sdc_stan=sdc_standard_step(M, ts, te, prob)
      u = np.zeros((M,prob.dim))
      usub = np.zeros((M,P,prob.dim))
      fu = np.zeros((M,prob.dim))
      fu_sub = np.zeros((M,P,prob.dim))

    # run predictor
      sdc.predict(u0_, u, usub, fu, fu_sub)

      for k in range(nK):
        # run standard node sweep...
        u_old = copy.deepcopy(u)
        usub_old = copy.deepcopy(usub)
        sdc.sweep(u0_, u, usub, fu, fu_sub)
      u0_ = u[-1]
      res.append(u0_[0])
    uex = u0*np.exp(nu*tend) + b/nu*(np.exp(nu*tend)-1)
    err =  abs(uex - u[-1])
    print "mrsdc:", err,uex, u[-1]

    us0_=u0
    u      = np.zeros((M,1,1))
    u_     = np.zeros((M,1,1))
    fu     = np.zeros((M, 1))

    sdc_stan.predict(us0_, u_)
    for k in range(nK):
        sdc_stan.sweep(us0_, u, u_)
        u_    = copy.deepcopy(u)
    err =  abs(uex - u[-1])
    print "imex:",err,uex, u[-1]
