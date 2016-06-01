from problem import problem
import numpy as np
import unittest

class test_problem(unittest.TestCase):

  def setUp(self):
    M = np.random.randint(8)+2
    P = np.random.randint(3)+2
    [tleft, tright] = np.sort(np.random.rand(2))
    self.prob = problem(M, P, tleft, tright+1.0)

  def test_I_m_mp1_linear(self):
    # evaluate f(u(t)) = lambda*t at nodes
    fu = self.prob.coll.nodes
    fu = np.reshape(fu, (self.prob.M, 1))
    self.prob.update_I_m(fu)
    int_ex = ((self.prob.lambda_fast + self.prob.lambda_slow)/2.0)*(self.prob.coll.nodes[0]**2 - self.prob.coll.tleft**2)
    err    = abs(self.prob.I_m_mp1[0] - int_ex)  
    assert err<1e-10, ("Smat failed to correctly integrate linear function -- error: %5.3e" % err)
    for i in range(1,self.prob.M):
      int_ex = ((self.prob.lambda_fast + self.prob.lambda_slow)/2.0)*(self.prob.coll.nodes[i]**2 - self.prob.coll.nodes[i-1]**2)
      err    = abs(self.prob.I_m_mp1[i] - int_ex)
      assert err<1e-10, ("Smat failed to correctly integrate linear function -- error: %5.3e" % err)

  def test_I_p_pp1_linear(self):
    # evaluate f(u(t)) = lambda*t at sub-nodes
    fu = np.zeros((self.prob.M,self.prob.P))
    for m in range(self.prob.M):
      fu[m,:] = self.prob.coll_fast[m].nodes
    self.prob.update_I_p(fu)
    for m in range(self.prob.M):
      int_ex = ((self.prob.lambda_fast + self.prob.lambda_slow)/2.0)*(self.prob.coll_fast[m].nodes[0]**2 - self.prob.coll_fast[m].tleft**2)
      err    = abs(self.prob.I_p_pp1[m,0] - int_ex)
      assert err<1e-10, ("Smat in coll_fast failed to correctly integrate linear function -- error: %5.3e" % err)
      for p in range(1,self.prob.P):
        int_ex = ((self.prob.lambda_fast + self.prob.lambda_slow)/2.0)*(self.prob.coll_fast[m].nodes[p]**2 - self.prob.coll_fast[m].nodes[p-1]**2)
        err    = abs(self.prob.I_p_pp1[m,p] - int_ex)
        assert err<1e-10, ("Smat in coll_fast failed to correctly integrate linear function -- error: %5.3e" % err)

  def test_get_coll_solution(self):
    err = np.zeros((6,1))
    for M in range(2,8):
      self.prob = problem(M, 2, 0.0, 1.0)
      coll = self.prob.get_coll_solution(1.0)
      uend = self.prob.end_value(coll, 1.0)
      err[M-2,0]  = abs(uend - np.exp(self.prob.lamb))
    assert np.min(err[0:5]/err[1:6])>10, "Error in end value generated from collocation solution seems not to decay geometrically in M"

  def test_get_coll_solution_sub(self):
    M = 4
    tend = 7.6
    err  = np.zeros((6,M))
    for P in range(2,8):
      self.prob  = problem(M, P, 0.0, tend)
      coll = self.prob.get_coll_solution_sub(1.0)
      for m in range(0,M):
        # Compute exact starting value for interval
        u0   = np.exp(self.prob.lamb*self.prob.coll_fast[m].tleft)
        uex  = np.exp(self.prob.lamb*self.prob.coll_fast[m].tright)
        uend = self.prob.end_value(np.reshape(coll[m,:], (P, 1)), u0, coll=self.prob.coll_fast[m])
        err[P-2,m] = abs(uend - uex)
    for m in range(0,M):
      assert np.min(err[0:5,m]/err[1:6,m])>9, "Error in values generated from sub collocation solution seems not to decay geometrically in P"

  def test_def_end_value_sub(self):
    tend = 7.3
    M = 4
    err = np.zeros((6,M))
    for P in range(2,8):
      self.prob  = problem(M, P, 0.0, tend)
      for m in range(0,M):
        # compute exact starting value for sub interval
        coll       = self.prob.get_coll_solution_sub(1.0)
        uend       = self.prob.end_value_sub(coll, 1.0, m)
        err[P-2,m] = abs(uend - np.exp(self.prob.lamb*self.prob.coll_fast[m].tright))
    for m in range(0,M):
      assert np.min(err[0:5,m]/err[1:6,m])>9, "Error in values generated from sub collocation solution seems not to decay geometrically in P"

  def test_get_residual(self):
    coll = self.prob.get_coll_solution(1.0)
    res  = self.prob.get_residual(coll, 1.0)
    assert res<1e-14, "Residual of collocation solution is not zero"

  def test_get_residual_sub(self):
    coll = self.prob.get_coll_solution_sub(1.0)
    res  = self.prob.get_residual_sub(coll, 1.0)
    assert res<1e-14, "Residuals of sub collocation solutions are not zero"

