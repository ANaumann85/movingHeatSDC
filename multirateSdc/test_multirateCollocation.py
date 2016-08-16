from MultirateCollocation import multirateCollocation
import unittest
import numpy as np

class test_problem(unittest.TestCase):

  def setUp(self):
    self.M = np.random.randint(8)+2
    self.P = np.random.randint(3)+2
    [self.tleft, self.tright] = np.sort(np.random.rand(2))
    self.tright += 1.0

  def test_caninstantiate(self):
    mc_coll = multirateCollocation(self.M, self.P, self.tleft, self.tright, 1)

  # Sum of embedded sub step lengths equals length of standard step
  def test_timesteps_match(self):
    mc_coll   = multirateCollocation(self.M, self.P, self.tleft, self.tright, 1)
    for m in range(self.M):
      dt_sum_sub = 0.0
      for p in range(self.P):
        dt_sum_sub += mc_coll.coll_sub[m].delta_m[p]
      err = abs(dt_sum_sub - mc_coll.coll.delta_m[m])
      assert err<1e-14, ("Sum of sub-steps did not yield length of step on coarse level. Error: %5.3e" % err)

  # Sum of embedded weights for function at standard nodes matches standard weights
  def test_weights_match(self):
    mc_coll = multirateCollocation(self.M, self.P, self.tleft, self.tright, 1)
    Smat    = mc_coll.coll.Smat
    Smat    = Smat[1:,1:]
    S_mnp   = mc_coll.S_mnp
    for m in range(self.M):
      for j in range(self.M):
        s = Smat[m,j]
        stilde = 0.0
        for p in range(self.P):
          stilde += S_mnp[m,j,p]
        err = abs(s - stilde)
        assert err<1e-14, ("Mismatch between s and s_tilde weights. Error %5.3e" % err)

  # Sum over integrals over embedded sub steps matches integral over standard sub step
  def test_integrates_match(self):
    dim = np.random.randint(1,10)
    mc_coll   = multirateCollocation(self.M, self.P, self.tleft, self.tright, dim)
    slope     = np.random.rand(dim,1)
    intercept = np.random.rand(dim,1)

    fu        = np.zeros((self.M,dim))
    for d in range(dim):
      fu[:,d] = slope[d]*mc_coll.coll.nodes + intercept[d]
    fusub   = np.zeros((self.P,dim))
    for m in range(self.M):
      if m==0:
        ta = mc_coll.coll.tleft
      else:
        ta = mc_coll.coll.nodes[m-1]
      tb   = mc_coll.coll.nodes[m]
      int_m_mp1 = mc_coll.integrate_m_mp1(fu, m)
      int_p_pp1 = np.zeros((dim,1))
      for d in range(dim):
        fusub[:,d] = slope[d]*mc_coll.coll_sub[m].nodes + intercept[d]
      for p in range(self.P):
        int_p_pp1 += mc_coll.integrate_p_pp1_sub(fusub, m, p)
      err = np.linalg.norm(int_m_mp1 - int_p_pp1, np.inf)
      assert err<1e-14, ("Sum of integrate_p_pp1 does not match integrate_m_mp1. Error: %5.3e" % err)

  # Linear function given at standard nodes is integrated exactly over embedded sub step
  def test_integrate_p_p1_linear(self):
    mc_coll = multirateCollocation(self.M, self.P, self.tleft, self.tright, 1)
    # Make sure method throws if input data does not fit
    with self.assertRaises(Exception):
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

  # Linear function given at embedded nodes is integrated exactly over embedded sub steps
  def test_integrate_p_pp1_sub_linear(self):
    mc_coll   = multirateCollocation(self.M, self.P, self.tleft, self.tright, 1)
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

  # Linear function at standard nodes is integrated exactly over standard sub steps
  def test_integrate_m_mp1_linear(self):
    mc_coll   = multirateCollocation(self.M, self.P, self.tleft, self.tright, 1)
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

  # Linear function given at embedded nodes is integrated exactly over standard sub steps
  def test_integrate_m_mp1_sub_linear(self):
    mc_coll   = multirateCollocation(self.M, self.P, self.tleft, self.tright, 1)
    slope     = np.random.rand(1)
    intercept = np.random.rand(1)
    for m in range(self.M):
      fu      = slope*mc_coll.coll_sub[m].nodes + intercept
      if m==0:
        ta = mc_coll.coll.tleft
      else:
        ta = mc_coll.coll.nodes[m-1]
      tb   = mc_coll.coll.nodes[m]
      intval  = mc_coll.integrate_m_mp1_sub(fu, m)
      intex   = 0.5*slope*(tb**2-ta**2) + intercept*(tb - ta)
      err     =  abs(intval - intex)
      assert err<1e-14, ("Function integrate_m_mp1_sub failed to integrate linear function excatly. Error: %5.3e" % err)
