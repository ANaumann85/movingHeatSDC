from problem_inh_const import problem_inh_const as problem_inh
from sdc import sdc_step
from sdc_standard import sdc_standard_step
import numpy as np
import unittest
import copy

nu = -4.0
b=2.0
prob = problem_inh(nu, b)
M = 3
P = 8
tstart = 0.0
tend   = 1.0
nIter=1
nKVec=range(0,5)
u0 = -b/nu

dt=(tend-tstart)/nIter
for nK in nKVec:
    u0_=u0
    us0_=u0
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

      us      = np.zeros((M,1,1))
      us_     = np.zeros((M,1,1))
      fsu     = np.zeros((M, 1))

    # run predictor
      sdc.predict(u0_, u, usub, fu, fu_sub)

      for k in range(nK):
        # run standard node sweep...
        u_old = copy.deepcopy(u)
        usub_old = copy.deepcopy(usub)
        sdc.sweep(u0_, u, usub, fu, fu_sub)
        #print sdc.residual(u0_, u)
      u0_ = u[-1]

      sdc_stan.predict(us0_, us_)
      for k in range(nK):
        sdc_stan.sweep(us0_, us, us_)
        us_    = copy.deepcopy(us)

    uex = u0*np.exp(nu*tend) + b/nu*(np.exp(nu*tend)-1)
    err =  abs(uex - u[-1])
    print "mrsdc:", err,uex, u[-1]

    err =  abs(uex - us[-1])
    print "imex:",err.flatten(),uex, us[-1]
