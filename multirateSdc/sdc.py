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
    self.I_p_pp1 = np.zeros((M,P,problem.dim))
    self.prob    = problem
    self.dt      = tend-tstart

  '''
  '''
  def update_I_m_mp1(self, u, usub):
    try:
      u     = np.reshape(u, (self.coll.M,self.prob.dim))
      usub  = np.reshape(usub, (self.coll.M, self.coll.P, self.prob.dim))
    except:
      raise
    
    for m in range(self.coll.M):
      self.I_m_mp1[m,:] = self.coll.integrate_m_mp1(self.prob.f1(u), m) + self.coll.integrate_m_mp1_sub(self.prob.f2(usub[m,:,:]), m)

  '''
  '''
  def update_I_p_pp1(self, u, usub):
    try:
      u     = np.reshape(u, (self.coll.M,self.prob.dim))
      usub  = np.reshape(usub, (self.coll.M, self.coll.P, self.prob.dim))
    except:
      raise
  
    for m in range(self.coll.M):
      for p in range(self.coll.P):
        self.I_p_pp1[m,p,:] = self.coll.integrate_p_pp1( self.prob.f1(u), m, p ) + self.coll.integrate_p_pp1_sub(self.prob.f2(usub[m,:,:]), m, p )

  '''
  '''
  def residual(self, u0, u):
    try:
      u = np.reshape(u, (self.coll.M,self.prob.dim))
    except:
      raise
    res = np.linalg.norm(u[0,:] - u0 - self.I_m_mp1[0,:], np.inf)
    for m in range(1,self.coll.M):
      resm = abs(u[m,:] - u[m-1,:] - self.I_m_mp1[m,:])
      res  = max(res, resm)
    return res

  '''
  '''
  def sub_residual(self, u0, usub):
    try:
      usub = np.reshape(usub, (self.coll.M, self.coll.P, self.prob.dim))
    except:
      raise
    res = 0.0
    for m in range(self.coll.M):
      res = max(res, self.sub_residual_m(u0, usub[m,:,:], m))
      u0 = usub[m,-1,:]
    return res

  '''
  '''
  def sub_residual_m(self, u0, u, m):
    try:
      u = np.reshape(u, (self.coll.P,self.prob.dim))
    except:
      raise TypeError("Failed to convert argument u into shape (P,1,dim)")
    res = np.linalg.norm(u[0,:] - u0 - self.I_p_pp1[m,0,:], np.inf)
    for p in range(1,self.coll.P):
      res = max(res, np.linalg.norm(u[p,:] - u[p-1,:] - self.I_p_pp1[m,p,:], np.inf))
    return res

  '''
  '''
  def predict(self, u0, u, usub):
    try:
      u = np.reshape(u, (self.coll.M, self.prob.dim))
      usub = np.reshape(usub, (self.coll.M, self.coll.P, self.prob.dim))
    except:
      raise

    u0_step = u0

    for m in range(self.coll.M):
      
        # standard step
        u[m,:] = self.prob.solve_f1(self.coll.coll.delta_m[m], u0_step)
        # embedded steps
        usub[m,0,:] = u0_step + self.coll.coll_sub[m].delta_m[0]*( self.prob.f1(u[m,:]) + self.prob.f2(u0_step) )
        for p in range(1,self.coll.P):
          usub[m,p,:] = usub[m,p-1,:] + self.coll.coll_sub[m].delta_m[p]*( self.prob.f1(u[m,:]) + self.prob.f2(usub[m,p-1,:]) )

        # overwrite standard value
        u0_step = usub[m,self.coll.P-1,:]
        u[m,:]  = usub[m,self.coll.P-1,:]

    return u, usub


  '''
  '''
  def sweep(self, u0, u, usub, u_, usub_):
    try:
      u     = np.reshape(u, (self.coll.M,self.prob.dim))
      u_    = np.reshape(u_, (self.coll.M,self.prob.dim))
      usub  = np.reshape(usub, (self.coll.M, self.coll.P, self.prob.dim))
      usub_ = np.reshape(usub_, (self.coll.M, self.coll.P, self.prob.dim))

    except:
      raise
    
    # update integral terms
    self.update_I_m_mp1(u_, usub_)
    self.update_I_p_pp1(u_, usub_)

    u0_step = u0

    for m in range(self.coll.M):
      
      # standard step
      rhs  = u0_step - self.coll.coll.delta_m[m]*( self.prob.f1(u_[m,:]) ) + self.I_m_mp1[m,:]
      u[m,:] = self.prob.solve_f1(self.coll.coll.delta_m[m], rhs)
      
      # embedded steps
      usub[m,0,:] = u0_step + self.coll.coll_sub[m].delta_m[0]*( self.prob.f1(u[m,:]) - self.prob.f1(u_[m,:]) + self.prob.f2(u0_step) - self.prob.f2(u0_step) ) + self.I_p_pp1[m,0,:]
      for p in range(1,self.coll.P):
        usub[m,p,:] = usub[m,p-1,:] + self.coll.coll_sub[m].delta_m[p]*( self.prob.f1(u[m,:]) - self.prob.f1(u_[m,:]) + self.prob.f2(usub[m,p-1,:]) - self.prob.f2(usub_[m,p-1,:]) ) + self.I_p_pp1[m,p,:]

      # overwrite standard value
      u0_step = usub[m,self.coll.P-1,:]
      u[m,:]  = usub[m,self.coll.P-1,:]

    return 0

  '''
  '''
  def get_collocation_solution(self, u0):  
    Q = self.coll.coll.Qmat
    Q = Q[1:,1:]
    try:
      Mcoll = np.eye(self.coll.M) - Q*self.prob.lambda_1
    except:
      raise
    assert self.prob.lambda_2==0.0, "Can only compute standard collocation solution analytically if lambda_2=0"
    return np.linalg.inv(Mcoll).dot(u0*np.ones(self.coll.M))

  def get_collocation_solution_sub(self, u0, m):  
    Q = self.coll.coll_sub[m].Qmat
    Q = Q[1:,1:]
    try:
      Mcoll = np.eye(self.coll.P) - Q*self.prob.lambda_2
    except:
      raise
    assert self.prob.lambda_1==0.0, "Can only compute embedded collocation solution analytically if lambda_1=0"
    return np.linalg.inv(Mcoll).dot(u0*np.ones(self.coll.P))