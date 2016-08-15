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
    self.setUp(lambda_2=0.0)
    u0   = np.random.rand(1)
    u    = np.zeros((self.M, 1))
    usub = np.zeros((self.M,self.P))
    self.sdc.predict(u0, u, usub)
    symb = 1.0
    for m in range(self.M):
      symb *= 1.0/(1.0 - self.prob.lambda_1*self.sdc.coll.coll.delta_m[m])
      err = np.linalg.norm(u[m,:] - symb*u0, np.inf)
      assert err<1e-14, ("Predict for standard nodes did not reproduce implicit Euler. Error: %5.3e" % err)

  # Predictor for embedded nodes is explicit Euler
  def test_predict_sub(self):
    self.setUp(lambda_1=0.0)
    u0 = np.random.rand(1)
    u    = np.zeros((self.M, 1))
    usub = np.zeros((self.M,self.P,1))
    # for lambda_1, u should be constant
    self.sdc.predict(u0, u, usub)
    symb = 1.0
    for m in range(self.M):
      for p in range(self.P):
        symb *= (1.0 + self.prob.lambda_2*self.sdc.coll.coll_sub[m].delta_m[p])
        err = abs(usub[m,p] - symb*u0)
        assert err<1e-14, ("sub_predict for lambda_1 = 0 failed to reproduce explicit Euler. Error: 5.3e" % err)

  '''
  Standard collocation solution has to be invariant under SDC sweep across standard nodes
  '''
  def test_sweep_coll_standard_invariant(self):
    self.setUp(lambda_2=0.0)
    u0         = np.random.rand(1)
    ucoll_     = self.sdc.get_collocation_solution(u0)
    ucoll_sub_ = np.zeros((self.M,self.P,1))
    
    ucoll     = np.zeros((self.M,1))
    ucoll_sub = np.zeros((self.M,self.P,1))
    
    self.sdc.sweep(u0, ucoll, ucoll_sub, ucoll_, ucoll_sub_)
    err = np.linalg.norm((ucoll - ucoll_).flatten(), np.inf)
    assert err<1e-14, ("Collocation solution not invariant under standard node SDC sweep with lambda_2=0. Error: %5.3e" % err)

  '''
  Standard collocation solution must have zero residual
  '''
  def test_collocation_residual(self):
    self.setUp(lambda_2=0.0)
    u0 = np.random.rand(1)
    ucoll = self.sdc.get_collocation_solution(u0)
    self.sdc.update_I_m_mp1(ucoll, np.zeros((self.M,self.P,1)))
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
  def test_converge_to_fixpoint(self):
    self.setUp(lambda_2=0.0)   
    u0    = 1.0
    u_    = np.zeros((self.M,1))
    u     = np.zeros((self.M,1))
    usub_ = np.zeros((self.M,self.P,1))
    usub  = np.zeros((self.M,self.P,1))
    # run predictor
    self.sdc.predict(u0, u_, usub_)

    for k in range(3):
      # update integral operators
      self.sdc.update_I_m_mp1(u_, usub_)
      self.sdc.update_I_p_pp1(u_, usub_)
      # run standard node sweep...
      self.sdc.sweep(u0, u, usub, u_, usub_)
      print ("Standard node update: %5.3e" % np.linalg.norm( (u - u_).flatten(), np.inf))
      print ("Embedded node update: %5.3e" % np.linalg.norm( (usub - usub_).flatten(), np.inf))
      print ("Standard node residual: %5.3e" % self.sdc.residual(u0, u))
      print ("Embedded node residual: %5.3e" % self.sdc.sub_residual(u0, usub))
      print ""
      u = copy.deepcopy(u_)
      usub = copy.deepcopy(usub_)

