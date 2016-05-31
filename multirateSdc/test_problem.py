from problem import problem
import numpy as np
import unittest

class test_problem(unittest.TestCase):

  def setUp(self):
    M = np.random.randint(8)+2
    P = np.random.randint(3)+2
    P = 2
    self.prob = problem(M, P)

  def test_I_m_mp1_linear(self):
    # evaluate f(u(t)) = lambda*t at nodes
    fu = self.prob.coll.nodes
    fu = np.reshape(fu, (self.prob.M, 1))
    self.prob.update_I_m(fu)
    for i in range(1,self.prob.M):
      int_ex = ((self.prob.lambda_fast + self.prob.lambda_slow)/2.0)*(self.prob.coll.nodes[i]**2 - self.prob.coll.nodes[i-1]**2)
      assert abs(self.prob.I_m_mp1[i] - int_ex)<1e-10, "Smat failed to correctly integrate linear function"

  @unittest.skip("not yet implemented")
  def test_I_p_pp1_linear(self):
    pass
