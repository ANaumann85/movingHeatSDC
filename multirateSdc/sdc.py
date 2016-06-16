import numpy as np
import copy
from MultirateCollocation import multirateCollocation
from problem import problem

class sdc_step():

  '''
  '''
  def __init__(self, M, P, tstart, tend, problem):
    self.coll    = multirateCollocation(M, P, tstart, tend)
    self.I_m_mp1 = np.zeros(M)
    self.I_p_pp1 = np.zeros((M,P))
    self.prob    = problem

  '''
  '''
  def update_I_m_mp1(self, u, usub):
    for m in range(self.coll.M):
      self.I_m_mp1[m] = self.coll.integrate_m_mp1(self.prob.f1(u), m) + coll.integrate_m_mp1_sub(self.prob.f2(usub[m,:]), m)

  '''
  '''
  def update_I_p_pp1(self, u, usub):
    for m in range(self.coll.M):
      for p in range(self.coll.P):
        self.I_p_pp1[m,p] = self.coll.integrate_p_pp1( self.prob.f1(u), m, p ) + self.coll.integrate_p_pp1_sub( self.prob.f2(usub[m,:]), m, p )

  '''
  '''
  def predict(self, u0):
    u = np.zeros(self.coll.M) 
    u[0] = self.prob.solve_f1(self.coll.coll.delta_m[0], u0)
    for m in range(1,self.coll.M):
      u[m] = self.prob.solve_f1(self.coll.coll.delta_m[m], u[m-1])
    return u

  '''
  '''
  def sub_predict(self, u0, u):
    usub = np.zeros((self.coll.M,self.coll.P))
    for m in range(self.coll.M):
      usub[m,0] = u0 + self.coll.coll_sub[m].delta_m[0]*( self.prob.f1(u[m]) + self.prob.f2(u0) )
      for p in range(1,self.coll.P):
        usub[m,p] = usub[m,p-1] + self.coll.coll_sub[m].delta_m[p]*( self.prob.f1(u[m]) + self.prob.f2(usub[m,p-1]) )       
      u0 = usub[m,-1]
    return usub

  '''
  '''
  def sweep(self, u0, u):
    pass

  '''
  '''
  def sub_sweep(self, u0, u, usub):
    pass
