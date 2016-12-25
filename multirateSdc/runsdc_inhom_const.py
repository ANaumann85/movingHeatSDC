from problem_inh_const import problem_inh_const as problem_inh
from sdc import sdc_step
from sdc_standard import sdc_standard_step
import numpy as np
import unittest
import copy

nu1 = -4.0
nu2 = -0.1
nu = nu1+nu2
b=10.0
prob = problem_inh(nu1, nu2, b)
M = 2
P = 8
tstart = 0.0
tend   = 1.0
nIter=2
nKVec=range(0,10)
nIterVec=[2**n for n in range(0,5)]
u0 = 0.5 #-b/nu
nK=4
dt=(tend-tstart)/nIter

#for nIter in nIterVec:
#    dt=(tend-tstart)/nIter
for nK in nKVec:
    u0_=u0
    us0_=u0
    res=[]
    for n in range(nIter):
      ts=tstart+n*dt
      te=ts+dt
      sdc = sdc_step(M, P, ts, te, prob,theta=1)
      sdc_stan=sdc_standard_step(M, ts, te, prob)
      u = np.zeros((M,prob.dim))
      usub = np.zeros((M,P,prob.dim))
      fu = np.zeros((M,prob.dim))
      fu_sub = np.zeros((M,P,prob.dim))

      us      = np.zeros((M,1))
      us_     = np.zeros((M,1))

    # run predictor
      sdc.predict(u0_, u, usub, fu, fu_sub)

      for k in range(nK):
        # run standard node sweep...
        #u_old = copy.deepcopy(u)
        #usub_old = copy.deepcopy(usub)
        sdc.sweep(u0_, u, usub, fu, fu_sub)
        #print sdc.residual(u0_, u)
      u0_ = copy.deepcopy(u[-1])

      sdc_stan.predict(us0_, us_)
      us = copy.deepcopy(us_)
      for k in range(nK):
        sdc_stan.sweep(us0_, us, us_)
        us_    = copy.deepcopy(us)
        #print "resid:",sdc_stan.residual(us0_, us)
      us0_=us[-1]

    uex = u0*np.exp(nu*tend) + b/nu*(np.exp(nu*tend)-1)
    errM =  abs(uex - u[-1])
    #print "mrsdc:", nK, errM,uex, u[-1]

    errS =  abs(uex - us[-1])
    #print "imex: ", nK,errS.flatten(),uex, us[-1]
    print "err:", nK, errM.flatten(), errS.flatten()
