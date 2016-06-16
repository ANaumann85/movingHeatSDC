from problem import problem
from sdc import sdc_step
import numpy as np
import unittest

class test_sdc_step(unittest.TestCase):

  def setUp(self, lambda_1=None, lambda_2=None):
    if lambda_1 is None: 
      lambda_1 = -1.0
    if lambda_2 is None:
      lambda_2 = -0.1
    self.prob = problem(lambda_1, lambda_2)
    self.M = 2
    self.P = 2
    tstart = 0.0
    tend   = 1.0
    self.sdc = sdc_step(self.M, self.P, tstart, tend, self.prob)

  def test_predict(self):
    u0 = np.random.rand(1)
    u  = self.sdc.predict(u0)
    symb = 1.0
    for m in range(self.M):
      symb *= 1.0/(1.0 - self.prob.lambda_1*self.sdc.coll.coll.delta_m[m])
      err = abs(u[m] - symb*u0)
      assert err<1e-14, ("Predict on coarse level did not reproduce implicit Euler. Error: %5.3e" % err)

  def test_predict_sub(self):
    self.setUp(lambda_1=0.0)
    u0 = np.random.rand(1)
    # for lambda_1, u should be constant
    u  = self.sdc.predict(u0)   
    err = np.linalg.norm(u - u0*np.ones(self.M), np.inf)
    assert err<1e-14, ("predict for lambda_1 = 0 failed to generate a constant solution. Err: %5.3e" % err)
    assert np.linalg.norm(self.prob.f1(u), np.inf)<1e-14
    usub = self.sdc.sub_predict(u0, u)
    symb = 1.0
    for m in range(self.M):
      for p in range(self.P):
        symb *= (1.0 + self.prob.lambda_2*self.sdc.coll.coll_sub[m].delta_m[p])
        err = abs(usub[m,p] - symb*u0)
        assert err<1e-14, ("sub_predict for lambda_1 = 0 failed to reproduce explicit Euler. Error: 5.3e" % err)

  '''
  Collocation solution has to be invariant under SDC sweeps
  '''
  def test_sweep_coll_invariant(self):
    self.setUp(lambda_2=0.0)
    u0 = np.random.rand(1)
    Q = self.sdc.coll.coll.Qmat
    Q = Q[1:,1:]
    Mcoll = np.eye(self.M) - self.sdc.dt*Q*self.prob.lambda_1
    ucoll = np.linalg.inv(Mcoll).dot(u0*np.ones(self.M))
    self.sdc.update_I_m_mp1(ucoll, np.zeros((self.M,self.P)))
    usweep = self.sdc.sweep(u0, ucoll)
    err = np.linalg.norm(usweep - ucoll, np.inf)
    assert err<1e-14, ("Collocation solution not invariant under coarse level SDC sweep with lambda_2=0. Error: %5.3e" % err)

  def test_subsweep_coll_invariant(self):
    self.setUp(lambda_1=0.0)
    u0 = np.random.rand(1)
    u0_step = u0
    ucoll = np.zeros((self.M, self.P))
    for m in range(self.M):
      Q          = self.sdc.coll.coll_sub[m].Qmat
      dt         = self.sdc.coll.coll_sub[m].tright - self.sdc.coll.coll_sub[m].tleft
      Q          = Q[1:,1:]
      Mcoll      = np.eye(self.P) - dt*Q*self.prob.lambda_2
      ucoll[m,:] = np.linalg.inv(Mcoll).dot(u0_step*np.ones(self.P))
      u0_step    = ucoll[m,-1]
    self.sdc.update_I_p_pp1(np.zeros(self.M), ucoll)
    usweep = self.sdc.sub_sweep(u0, np.zeros(self.M), np.zeros(self.M), ucoll)
    err = np.linalg.norm(usweep[0,:] - ucoll[0,:], np.inf)
    print err
