from problem import problem
from problem_model import problem_model
from problem_inh import problem_inh
from sdc_standard import sdc_standard_step
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
    tstart = 0.0
    tend   = 0.5
    self.sdc = sdc_standard_step(self.M, tstart, tend, self.prob)

  # Predictor for lambda_2=0 is implicit Euler
  def test_predict(self):
    self.setUp(lambda_2=0.0)
    u0   = np.random.rand(1)
    u    = np.zeros((self.M, 1))
    self.sdc.predict(u0, u)
    symb = 1.0
    for m in range(self.M):
      symb *= 1.0/(1.0 - self.prob.lambda_1*self.sdc.coll.delta_m[m])
      err = np.linalg.norm(u[m,:] - symb*u0, np.inf)
      assert err<1e-14, ("predict for standard nodes did not reproduce implicit Euler. Error: %5.3e" % err)

  # Predictor for lambda_1=0 is explicit Euler
  def test_predict_sub(self):
    self.setUp(lambda_1=0.0)
    u0 = np.random.rand(1)
    u    = np.zeros((self.M, 1))
    # for lambda_1, u should be constant
    self.sdc.predict(u0, u)
    symb = 1.0
    for m in range(self.M):
      symb *= (1.0 + self.prob.lambda_2*self.sdc.coll.delta_m[m])
      err = abs(u[m,:] - symb*u0)
      assert err<1e-14, ("predict for lambda_1 = 0 failed to reproduce explicit Euler. Error: 5.3e" % err)

  '''
  Standard collocation solution has to be invariant under SDC sweep across standard nodes
  '''
  def test_sweep_coll_standard_invariant(self):
    self.setUp()
    u0         = np.random.rand(1)
    u0 = 1.0
    ucoll_     = self.sdc.get_collocation_solution(u0)
    ucoll     = np.zeros((self.M,1))
    self.sdc.sweep(u0, ucoll, ucoll_)
    err = np.linalg.norm((ucoll - ucoll_).flatten(), np.inf)
    assert err<1e-14, ("Collocation solution not invariant under standard node SDC sweep with lambda_2=0. Error: %5.3e" % err)

  '''
  Standard collocation solution must have zero residual
  '''
  def test_collocation_residual(self):
    self.setUp()
    u0 = np.random.rand(1)
    ucoll = self.sdc.get_collocation_solution(u0)
    self.sdc.update_I_m_mp1(ucoll)
    res = self.sdc.residual(u0, ucoll)
    assert res<1e-14, ("Residual not zero for collocation solution. Error: %5.3e" % res)

  ''' 
  '''
  def test_converge_to_fixpoint_inh(self):
    nu = -1.0
    self.prob = problem_inh(nu)
    self.M = 8
    tstart = 0.0
    tend   = 0.1
    self.sdc = sdc_standard_step(self.M, tstart, tend, self.prob)
  
    u0 = 1.0
    u_ = np.zeros((self.M,self.prob.dim))
    u = np.zeros((self.M,self.prob.dim))


    # run predictor
    self.sdc.predict(u0, u_)

    for k in range(15):
      # run standard node sweep...
      self.sdc.sweep(u0, u, u_)
      update_standard = np.linalg.norm( (u-u_).flatten(), np.inf)
      res_standard    = self.sdc.residual(u0, u)
      u_              = copy.deepcopy(u)

    c1  = u0 + 1.0/(nu**2+1)
    uex = c1*np.exp(nu*tend) - (nu*np.sin(tend) + np.cos(tend))/(nu**2+1)
    err =  abs(uex - u[-1])
    assert err<1e-13, ("Larger than expected error for inhomogenous problem. Error: %5.3e" % err)
    assert update_standard<1e-13, ("Standard update failed to converge to zero. Value: %5.3e" % update_standard)
    assert res_standard<1e-13, ("Standard residual failed to converge to zero. Value: %5.3e" % res_standard)
  
  ''' 
  '''
  def test_converge_to_fixpoint(self):
    self.prob = problem_model(a = -0.1, nu = -1.0, alpha = 1.0, v0 = 1.0)
    self.M = 3
    tstart = 0.0
    tend   = 0.2
    self.sdc = sdc_standard_step(self.M, tstart, tend, self.prob)
    
    u0    = np.reshape([2.0, 1.0, 0.0, 1.0, 0.0], (self.prob.dim,))
    u_    = np.zeros((self.M,self.prob.dim))
    u     = np.zeros((self.M,self.prob.dim))
    # run predictor
    self.sdc.predict(u0, u_)

    for k in range(25):
      # run standard node sweep...
      self.sdc.sweep(u0, u, u_)
      update_standard = np.linalg.norm( (u-u_).flatten(), np.inf)
      res_standard    = self.sdc.residual(u0, u)
      u_    = copy.deepcopy(u)

    assert update_standard<1e-12, ("Standard update failed to converge to zero. Value: %5.3e" % update_standard)
    assert res_standard<1e-12, ("Standard residual failed to converge to zero. Value: %5.3e" % res_standard)
