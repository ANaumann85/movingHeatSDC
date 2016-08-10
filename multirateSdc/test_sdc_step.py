from problem import problem
from sdc import sdc_step
import numpy as np
import unittest
import copy

class test_sdc_step(unittest.TestCase):

  def setUp(self, lambda_1=None, lambda_2=None):
    if lambda_1 is None: 
      lambda_1 = -1.0
    if lambda_2 is None:
      lambda_2 = -0.1
    self.prob = problem(lambda_1, lambda_2)
    self.M = 2
    self.P = 3
    tstart = 0.0
    tend   = 0.5
    self.sdc = sdc_step(self.M, self.P, tstart, tend, self.prob)

  # Predictor for standard nodes is implicit Euler
  def test_predict(self):
    u0 = np.random.rand(1)
    u  = self.sdc.predict(u0)
    symb = 1.0
    for m in range(self.M):
      symb *= 1.0/(1.0 - self.prob.lambda_1*self.sdc.coll.coll.delta_m[m])
      err = abs(u[m] - symb*u0)
      assert err<1e-14, ("Predict on coarse level did not reproduce implicit Euler. Error: %5.3e" % err)

  # Predictor for embedded nodes is explicit Euler
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
  Standard collocation solution has to be invariant under SDC sweep across standard nodes
  '''
  def test_sweep_coll_invariant(self):
    self.setUp(lambda_2=0.0)
    u0 = np.random.rand(1)
    ucoll = self.sdc.get_collocation_solution(u0)
    self.sdc.update_I_m_mp1(ucoll, np.zeros((self.M,self.P)))
    usweep = self.sdc.sweep(u0, ucoll)
    err = np.linalg.norm(usweep - ucoll, np.inf)
    assert err<1e-14, ("Collocation solution not invariant under standard node SDC sweep with lambda_2=0. Error: %5.3e" % err)

  '''
  Standard collocation solution must have zero residual
  '''
  def test_collocation_residual(self):
    self.setUp(lambda_2=0.0)
    u0 = np.random.rand(1)
    ucoll = self.sdc.get_collocation_solution(u0)
    self.sdc.update_I_m_mp1(ucoll, np.zeros((self.M,self.P)))
    res = self.sdc.residual(u0, ucoll)
    assert res<1e-14, ("Residual not zero for collocation solution. Error: %5.3e" % res)

  '''
  Embedded collocation solution must have zero residual
  '''
  def test_sub_collocation_residual_p(self):
    self.setUp(lambda_1 = 0.0, lambda_2=-1.0)
    u0 = 1.0
    for m in range(self.M):
      ucoll     = self.sdc.get_collocation_solution_sub(u0, m)
      usub      = np.zeros((self.M,self.P))
      usub[m,:] = ucoll
      self.sdc.update_I_p_pp1(np.zeros(self.M), usub)
      res = self.sdc.sub_residual_m(u0, ucoll, m)
      assert res<1e-14, ("Individual sub-step collocation solution failed to result in zero sub-step residual. Error: %5.3e" % res)

  '''
  Make sure the collocation solution computed from matrix inversion leads to a zero residual with respect to the integral operators I_m_mp1 and I_p_pp1
  '''
  def test_sub_collocation_residual(self):
    self.setUp(lambda_1 = 0.0, lambda_2=-1.0)
    u0     = np.random.rand(1)
    usub   = np.zeros((self.M,self.P))
    u0_sub = u0
    for m in range(self.M):
      usub[m,:]  = self.sdc.get_collocation_solution_sub(u0_sub, m)
      u0_sub     = usub[m,-1]
    self.sdc.update_I_p_pp1(np.zeros(self.M), usub)
    res = self.sdc.sub_residual(u0, usub)
    assert res<1e-14, ("Solution composed of sub-step collocation solutions failed to produce zero for sub_residual. Error: %5.3e" % res)

  '''
  '''
  def test_subsweep_coll_invariant(self):
    self.setUp(lambda_1=0.0, lambda_2=-1.0)
    u0 = np.random.rand(1)
    u0_step = u0
    ucoll = np.zeros((self.M, self.P))
    for m in range(self.M):
      ucoll[m,:] = self.sdc.get_collocation_solution_sub(u0_step, m)
      u0_step    = ucoll[m,-1]
    self.sdc.update_I_p_pp1(np.zeros(self.M), ucoll)
    usweep = self.sdc.sub_sweep(u0, np.zeros(self.M), np.zeros(self.M), ucoll)
    err = np.linalg.norm(usweep - ucoll, np.inf)
    assert err<1e-14, ("Solution composed of sub-step collocation solutions not invariant under sub-step sweep. Error: %5.3e" % err)

  ''' 
  '''
  def test_converge_to_fixpoint(self):
    self.setUp(lambda_2=0.0)   
    u0    = 1.0
    # run standard node predictor...
    u_    = self.sdc.predict(u0)
    # ... then fill in embedded nodes
    usub_ = self.sdc.sub_predict(u0, u_)
    for k in range(3):
      # update integral operators
      self.sdc.update_I_m_mp1(u_, usub_)
      self.sdc.update_I_p_pp1(u_, usub_)
      # run standard node sweep...
      u    = self.sdc.sweep(u0, u_)
      # ...followed by embedded node sweep
      usub = self.sdc.sub_sweep(u0, u, u_, usub_)
      print ("Standard node update: %5.3e" % np.linalg.norm(u - u_, np.inf))
      print ("Embedded node update: %5.3e" % np.linalg.norm(usub - usub_, np.inf))
      print ("Standard node residual: %5.3e" % self.sdc.residual(u0, u))
      print ("Embedded node residual: %5.3e" % self.sdc.sub_residual(u0, usub))
      print ""
      u = copy.deepcopy(u_)
      usub = copy.deepcopy(usub_)

