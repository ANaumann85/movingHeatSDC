from problem import problem
from problem_model import problem_model
from problem_inh import problem_inh
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
    u0 = 1.0
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
    fu, fu_sub = self.sdc.evaluate_f(ucoll, np.zeros((self.M,self.P,1)))
    self.sdc.update_I_m_mp1(fu, fu_sub)
    res = self.sdc.residual(u0, ucoll)
    assert res<1e-14, ("Residual not zero for collocation solution. Error: %5.3e" % res)

  '''
  Embedded collocation solution must have zero residual
  '''
  def test_sub_collocation_residual_p(self):
    self.setUp(lambda_1 = 0.0, lambda_2=-1.0)
    u0 = 1.0
    for m in range(self.M):
      ucoll         = self.sdc.get_collocation_solution_sub(u0, m)
      usub          = np.zeros((self.M,self.P,self.prob.dim))
      usub[m,:,:]   = ucoll
      fu, fu_sub = self.sdc.evaluate_f(np.zeros((self.M,self.prob.dim)), usub)
      self.sdc.update_I_p_pp1(fu, fu_sub)
      res = self.sdc.sub_residual_m(u0, ucoll, m)
      assert res<1e-14, ("Individual sub-step collocation solution failed to result in zero sub-step residual. Error: %5.3e" % res)

  '''
  Make sure the collocation solution computed from matrix inversion leads to a zero residual with respect to the integral operators I_m_mp1 and I_p_pp1
  '''
  def test_sub_collocation_residual(self):
    self.setUp(lambda_1 = 0.0, lambda_2=-1.0)
    u0     = np.random.rand(1)
    usub   = np.zeros((self.M,self.P,self.prob.dim))
    u0_sub = u0
    for m in range(self.M):
      usub[m,:]  = self.sdc.get_collocation_solution_sub(u0_sub, m)
      u0_sub     = usub[m,-1]
    fu, fu_sub = self.sdc.evaluate_f(np.zeros((self.M,self.prob.dim)), usub)
    self.sdc.update_I_p_pp1(fu, fu_sub)
    res = self.sdc.sub_residual(u0, usub)
    assert res<1e-14, ("Solution composed of sub-step collocation solutions failed to produce zero for sub_residual. Error: %5.3e" % res)

  '''
  '''
  def test_collocation_update(self):
    self.setUp(lambda_1 = 0.0, lambda_2=-1.0)
    u0     = np.random.rand(1)
    usub   = np.zeros((self.M,self.P,self.prob.dim))
    u0_sub = u0
    for m in range(self.M):
      usub[m,:]  = self.sdc.get_collocation_solution_sub(u0_sub, m)
      u0_coll    = self.sdc.collocation_update(u0_sub, np.zeros((self.M,self.prob.dim)), usub, m)
      u0_sub     = usub[m,-1]
      err        = np.linalg.norm(u0_coll - u0_sub)
      assert err<1e-14, ("For collocation solution, update formula failed to reproduce last stage. Error: %5.3e" % err)
      
  '''
  '''
  def test_converge_to_fixpoint_inh(self):
    nu = -1.0
    self.prob = problem_inh(nu)
    self.M = 5
    self.P = 5
    tstart = 0.0
    tend   = 0.5
    self.sdc = sdc_step(self.M, self.P, tstart, tend, self.prob)
  
    u0 = 1.0
    u_ = np.zeros((self.M,self.prob.dim))
    u = np.zeros((self.M,self.prob.dim))
    usub_ = np.zeros((self.M,self.P,self.prob.dim))
    usub = np.zeros((self.M,self.P,self.prob.dim))

    # run predictor
    self.sdc.predict(u0, u_, usub_)

    for k in range(15):
      # run standard node sweep...
      self.sdc.sweep(u0, u, usub, u_, usub_)
      update_standard = np.linalg.norm( (u-u_).flatten(), np.inf)
      update_embedded = np.linalg.norm( (usub-usub_).flatten(), np.inf)
      res_standard    = self.sdc.residual(u0, u)
      res_embedded    = self.sdc.sub_residual(u0, usub)
      u_    = copy.deepcopy(u)
      usub_ = copy.deepcopy(usub)

    c1  = u0 + 1.0/(nu**2+1)
    uex = c1*np.exp(nu*tend) - (nu*np.sin(tend) + np.cos(tend))/(nu**2+1)
    err =  abs(uex - u[-1])
    assert err<1e-13, ("Larger than expected error for inhomogenous problem. Error: %5.3e" % err)
    assert update_standard<1e-13, ("Standard update failed to converge to zero. Value: %5.3e" % update_standard)
    assert update_embedded<1e-13, ("Embedded update failed to converge to zero. Value: %5.3e" % update_embedded)
    assert res_standard<1e-13, ("Standard residual failed to converge to zero. Value: %5.3e" % res_standard)
    assert res_embedded<1e-13, ("Embedded residual failed to converge to zero. Value: %5.3e" % res_embedded)
  
  ''' 
  '''
  def test_converge_to_fixpoint(self):
    self.prob = problem_model(a = -0.1, nu = -1.0, alpha = 1.0, v0 = 1.0)
    self.M = 3
    self.P = 4
    tstart = 0.0
    tend   = 0.2
    self.sdc = sdc_step(self.M, self.P, tstart, tend, self.prob)
    
    u0    = np.reshape([2.0, 1.0, 0.0, 1.0, 0.0], (self.prob.dim,))
    u_    = np.zeros((self.M,self.prob.dim))
    u     = np.zeros((self.M,self.prob.dim))
    usub_ = np.zeros((self.M,self.P,self.prob.dim))
    usub  = np.zeros((self.M,self.P,self.prob.dim))
    # run predictor
    self.sdc.predict(u0, u_, usub_)

    for k in range(25):
      # run standard node sweep...
      self.sdc.sweep(u0, u, usub, u_, usub_)
      update_standard = np.linalg.norm( (u-u_).flatten(), np.inf)
      update_embedded = np.linalg.norm( (usub-usub_).flatten(), np.inf)
      res_standard    = self.sdc.residual(u0, u)
      res_embedded    = self.sdc.sub_residual(u0, usub)
      u_    = copy.deepcopy(u)
      usub_ = copy.deepcopy(usub)

    assert update_standard<1e-12, ("Standard update failed to converge to zero. Value: %5.3e" % update_standard)
    assert update_embedded<1e-12, ("Embedded update failed to converge to zero. Value: %5.3e" % update_embedded)
    assert res_standard<1e-12, ("Standard residual failed to converge to zero. Value: %5.3e" % res_standard)
    assert res_embedded<1e-12, ("Embedded residual failed to converge to zero. Value: %5.3e" % res_embedded)

    # Continue test to validate collocation update formula
    u0_sub = u0
    for m in range(self.M):
      u0_coll    = self.sdc.collocation_update(u0_sub, u, usub, m)
      u0_sub     = usub[m,-1]
      err        = np.linalg.norm(u0_coll - u0_sub)
      assert err<1e-12, ("For collocation solution, update formula failed to reproduce last stage. Error: %5.3e" % err)