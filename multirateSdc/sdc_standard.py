import numpy as np
import copy
from MultirateCollocation import multirateCollocation
from problem import problem

class sdc_step():

  '''
  '''
  def __init__(self, M, P, tstart, tend, problem):
    self.coll    = multirateCollocation(M, P, tstart, tend, problem.dim)
    self.I_m_mp1 = np.zeros((M,problem.dim))
    self.prob    = problem
    self.dt      = tend-tstart

  '''
  '''
  def update_I_m_mp1(self, u):
    try:
      u     = np.reshape(u, (self.coll.M,self.prob.dim))
    except:
      raise
    
    fu = self.evaluate_f(u)
    
    for m in range(self.coll.M):
      self.I_m_mp1[m,:] = self.coll.integrate_m_mp1(fu, m)


  '''
  '''
  def evaluate_f(self, u, usub):
    fu     = np.zeros((self.coll.M,self.prob.dim))
    for m in range(self.coll.M):
      fu[m,:] = self.prob.f1(u[m,:]) + self.prob.f2(u[m,:])
    return fu
    
  '''
  '''
  def residual(self, u0, u):
    try:
      u = np.reshape(u, (self.coll.M,self.prob.dim))
    except:
      raise
    res = np.linalg.norm(u[0,:] - u0 - self.I_m_mp1[0,:], np.inf)
    for m in range(1,self.coll.M):
      resm = np.linalg.norm(u[m,:] - u[m-1,:] - self.I_m_mp1[m,:], np.inf)
      res  = max(res, resm)
    return res


  '''
  '''
  def predict(self, u0, u, usub):
    try:
      u    = np.reshape(u, (self.coll.M, self.prob.dim))
      u0   = np.reshape(u0, (self.prob.dim,))
    except:
      raise

    for m in range(self.coll.M):
      
        # standard step (explicit part)
        b = u0 # ... complete
        
        # standard step (implicit part)
        u[m,:] = self.prob.solve_f1(self.coll.coll.delta_m[m], b)


    return u


  '''
  '''
  def sweep(self, u0, u, u_):
    try:
      u     = np.reshape(u, (self.coll.M,self.prob.dim))
      u_    = np.reshape(u_, (self.coll.M,self.prob.dim))

    except:
      raise
    
    # update integral terms
    self.update_I_m_mp1(u_)

    for m in range(self.coll.M):
      
      # standard step (explicit part)
      rhs  = u0_step - self.coll.coll.delta_m[m]*( self.prob.f1(u_[m,:]) ) + self.I_m_mp1[m,:]
      u[m,:] = self.prob.solve_f1(self.coll.coll.delta_m[m], rhs)
      
      # standard step (implicit part)

    return 0

  '''
  '''
  def get_collocation_solution(self, u0):  
    Q = self.coll.coll.Qmat
    Q = Q[1:,1:]
    try:
      Mcoll = np.eye(self.coll.M) - Q*(self.prob.lambda_1 + self.prob.lambda_2)
    except:
      raise
    return np.reshape( np.linalg.inv(Mcoll).dot(u0*np.ones(self.coll.M)), (self.coll.M, 1))
