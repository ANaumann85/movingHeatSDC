import numpy as np
from Collocation import *
from CollocationClasses import *

class problem():

  def __init__(self, M, P):
    self.I_m_mp1 = np.zeros((M,1))
    self.I_p_pp1 = np.zeros((M,P))
    self.lambda_fast = -1.0
    self.lambda_slow = -0.0
    self.M = M
    self.P = P
    # Consider a single time step [0.0, 1.0]
    self.coll = CollGaussRadau_Right(num_nodes = M, tleft = 0.0, tright = 1.0)
    self.coll_fast = []
    self.coll_fast.append(CollGaussRadau_Right(num_nodes=P, tleft = 0.0, tright = self.coll.nodes[0]))
    for i in range(1,M):
      self.coll_fast.append(CollGaussRadau_Right(num_nodes=P, tleft = self.coll.nodes[i-1] , tright=self.coll.nodes[i]))

  def fexpl(self, u, t = 0.0):
    return self.lambda_slow*u
    
  def fimpl(self, u, t = 0.0):
    return self.lambda_fast*u

  def solve(self, a, b, t = 0.0):
    return b/(1.0 - a*self.lambda_fast)

  def update_I_m(self, u):
    assert np.shape(u)==(self.M,1), "u must have shape Mx1"
    S = self.coll.Smat
    S = S[1:, 1:]
    for i in range(self.M):
      self.I_m_mp1[i] = 0.0
      for j in range(self.M):
        self.I_m_mp1[i] += S[i,j]*( self.fexpl(u[j]) + self.fimpl(u[j]) )
    return None

  def end_value(self, u, u0):
    uend = u0
    for m in range(0,self.M):
      uend += self.coll.weights[m]*(self.fexpl(u[m]) + self.fimpl(u[m]))
    return uend

  def print_nodes(self):
    for i in range(0,self.M):
      for j in range(0,self.P):
        print "Fast node: %5.3f" % self.coll_fast[i].nodes[j]
      print "Slow node: %5.3f" % self.coll.nodes[i]
      print "Delta: %5.3f" % self.coll.delta_m[i]
      print ""
