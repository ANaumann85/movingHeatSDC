from problem import problem
import numpy as np
import unittest

class test_problem(unittest.TestCase):

  def setUp(self):
    M = np.random.randint(8)+2
    P = np.random.randint(3)+2
    P = 2
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

