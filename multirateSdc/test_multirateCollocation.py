from MultirateCollocation import multirateCollocation
import unittest
import numpy as np

class test_problem(unittest.TestCase):

  def setUp(self):
    self.M = np.random.randint(8)+2
    self.P = np.random.randint(3)+2
    [self.tleft, self.tright] = np.sort(np.random.rand(2))

  def test_caninstantiate(self):
    mc_coll = multirateCollocation(self.M, self.P, self.tleft, self.tright)

  def test_integrate_p_p1_linear(self):
    mc_coll = multirateCollocation(self.M, self.P, self.tleft, self.tright)
    # Make sure method throws if input data does not fit
    with self.assertRaises(TypeError):
      mc_coll.integrate_p_pp1(np.zeros((self.M+1,1)), 0, 0)
    slope = np.random.rand(1)
    intercept = np.random.rand(1)
    fu      = slope*mc_coll.coll.nodes + intercept
    for m in range(self.M):
      for p in range(self.P):
        if p==0:
          ta = mc_coll.coll_sub[m].tleft
        else:
          ta = mc_coll.coll_sub[m].nodes[p-1]
        tb   = mc_coll.coll_sub[m].nodes[p]        
        intval  = mc_coll.integrate_p_pp1(fu, m, p) 
        intex   = 0.5*slope*(tb**2-ta**2) + intercept*(tb - ta)
        err     =  abs(intval - intex)
        assert err<1e-14, ("Function integrate_p_pp1 failed to integrate linear function excatly. Error: %5.3e" % err)

  def test_integrate_p_pp1_sub_linear(self):
    mc_coll   = multirateCollocation(self.M, self.P, self.tleft, self.tright)
    slope     = np.random.rand(1)
    intercept = np.random.rand(1)
    for m in range(self.M):
      fu      = slope*mc_coll.coll_sub[m].nodes + intercept
      for p in range(self.P):
        if p==0:
          ta = mc_coll.coll_sub[m].tleft
        else:
          ta = mc_coll.coll_sub[m].nodes[p-1]
        tb   = mc_coll.coll_sub[m].nodes[p]
        intval = mc_coll.integrate_p_pp1_sub(fu, m, p)
        intex   = 0.5*slope*(tb**2-ta**2) + intercept*(tb - ta)
        err     =  abs(intval - intex)
        assert err<1e-14, ("Function integrate_p_pp1_sub failed to integrate linear function excatly. Error: %5.3e" % err)

  def test_integrate_m_mp1_linear(self):
    mc_coll   = multirateCollocation(self.M, self.P, self.tleft, self.tright)
    slope     = np.random.rand(1)
    intercept = np.random.rand(1)
    for m in range(self.M):
      fu      = slope*mc_coll.coll.nodes + intercept
      if m==0:
        ta = mc_coll.coll.tleft
      else:
        ta = mc_coll.coll.nodes[m-1]
      tb   = mc_coll.coll.nodes[m]
      intval = mc_coll.integrate_m_mp1(fu, m)
      intex   = 0.5*slope*(tb**2-ta**2) + intercept*(tb - ta)
      err     =  abs(intval - intex)
      assert err<1e-14, ("Function integrate_m_mp1 failed to integrate linear function excatly. Error: %5.3e" % err)
