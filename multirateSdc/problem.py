import numpy as np
from Collocation import *
from CollocationClasses import *

class problem():

  def __init__(self, M, P, tleft, tright, lambda_slow, lambda_fast):
    assert tleft<tright, "tleft must be smaller than tright"
    self.I_m_mp1 = np.zeros((M,1))
    self.I_p_pp1 = np.zeros((M,P))
    self.lambda_fast = lambda_fast
    self.lambda_slow = lambda_slow
    self.lamb = self.lambda_fast + self.lambda_slow
    self.dt = abs(tright - tleft)
    self.M = M
    self.P = P
    # Consider a single time step [0.0, 1.0]
    self.coll = CollGaussRadau_Right(num_nodes = M, tleft = tleft, tright = tright)
    self.coll_fast = []
    self.coll_fast.append(CollGaussRadau_Right(num_nodes=P, tleft = tleft, tright = self.coll.nodes[0]))
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

  def update_I_p(self, u):
    assert np.shape(u)==(self.M, self.P), "u must have shape MxP"
    for m in range(self.M):
      S = self.coll_fast[m].Smat
      S = S[1:,1:]
      for i in range(self.P):
        self.I_p_pp1[m,i] = 0.0 
        for j in range(self.P):
          self.I_p_pp1[m,i] += S[i,j]*( self.fexpl(u[m,j]) + self.fimpl(u[m,j]) )
    return None

  def end_value(self, u, u0, coll=None):
    if coll==None:
      coll = self.coll
    M = coll.num_nodes    
    assert np.shape(u)==(M,1), "u must have shape Mx1"
    for m in range(0,M):
      u0 += coll.weights[m]*(self.fexpl(u[m]) + self.fimpl(u[m]))    
    return u0

  def end_value_sub(self, u, u0, m):
    assert m in range(0,self.M), "Must have 0<=m<=M-1"
    assert np.shape(u)==(self.M, self.P), "u must have shape MxP"   
    for i in range(0,m+1): 
      u0 = self.end_value(np.reshape(u[i,:], (self.P, 1)), u0, self.coll_fast[i])
    return u0

  def end_value_sub_all(self, u, u0):
    assert np.shape(u)==(self.M, self.P), "u must have shape MxP"    
    for m in range(self.M):      
       u0 = self.end_value(np.reshape(u[m,:], (self.P, 1)), u0, self.coll_fast[m])
    return u0

  def get_residual(self, u, u0, coll=None):
    if coll==None:
      coll = self.coll
    M = coll.num_nodes
    Q = coll.Qmat
    Q = Q[1:,1:]
    assert np.shape(u)==(M,1), "u must be (M,1) shape"
    e = np.ones((M,1))
    return np.linalg.norm(u - u0*e - self.lamb*Q.dot(u), np.inf)

  def get_residual_sub(self, u, u0):
    assert np.shape(u)==(self.M, self.P), "u must have shape MxP"   
    res = np.zeros((self.M, 1))
    for m in range(0,self.M):
      res[m,0] = self.get_residual( np.reshape(u[m,:], (self.P,1)), u0, self.coll_fast[m])
      u0 = self.end_value(np.reshape(u[m,:], (self.P,1)), u0, self.coll_fast[m])
    return np.linalg.norm(res, np.inf)
 
  def get_coll_solution(self, u0, coll=None):
    if coll==None:
      coll = self.coll
    M = coll.num_nodes
    Q = coll.Qmat
    Q = Q[1:,1:]
    e = np.ones((M,1))
    I = np.eye(M)
    return np.linalg.solve(I - self.lamb*Q, e*u0)

  def get_coll_solution_sub(self, u0):
    u = np.zeros((self.M,self.P))
    for m in range(self.M):
        u[m,:] = np.reshape(self.get_coll_solution(u0, self.coll_fast[m]), (1, self.P))
        u0     = self.end_value(np.reshape(u[m,:], (self.P,1)), u0, self.coll_fast[m])
    return u

  def print_nodes(self):
    for i in range(0,self.M):
      for j in range(0,self.P):
        print "Fast node: %5.3f" % self.coll_fast[i].nodes[j]
      print "Slow node: %5.3f" % self.coll.nodes[i]
      print "Delta: %5.3f" % self.coll.delta_m[i]
      print ""
