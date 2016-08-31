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
  def update_I_m_mp1(self, fu, fu_sub):
    try:
      fu     = np.reshape(fu, (self.coll.M,self.prob.dim))
      fu_sub  = np.reshape(fu_sub, (self.coll.M, self.coll.P, self.prob.dim))
    except:
      raise
    
    for m in range(self.coll.M):
      self.I_m_mp1[m,:] = self.coll.integrate_m_mp1(fu, m) + self.coll.integrate_m_mp1_sub(fu_sub[m,:,:], m)

  '''
  '''
  def update_I_p_pp1(self, fu, fu_sub):
    try:
      fu     = np.reshape(fu, (self.coll.M,self.prob.dim))
      fu_sub  = np.reshape(fu_sub, (self.coll.M, self.coll.P, self.prob.dim))
    except:
      raise
    
    for m in range(self.coll.M):
      for p in range(self.coll.P):
        self.I_p_pp1[m,p,:] = self.coll.integrate_p_pp1( fu, m, p ) + self.coll.integrate_p_pp1_sub( fu_sub[m,:,:], m, p )

  '''
  '''
  def evaluate_f(self, u, usub):
    fu     = np.zeros((self.coll.M,self.prob.dim))
    fu_sub = np.zeros((self.coll.M, self.coll.P, self.prob.dim))
    for m in range(self.coll.M):
      fu[m,:] = self.prob.f1(u[m,:])
      for p in range(self.coll.P):
        fu_sub[m,p,:] = self.prob.f2(usub[m,p,:], self.coll.coll_sub[m].nodes[p])
    return fu, fu_sub
    
  '''
  '''
  def collocation_update(self, u0, u, usub, m):
    try:
      usub  = np.reshape(usub, (self.coll.M, self.coll.P, self.prob.dim))
    except:
      raise
    
    # NOTE: this is inefficient since it updates ALL values
    fu, fu_sub = self.evaluate_f(u, usub)
    
    u0 += self.coll.integrate_m_mp1(fu, m) + self.coll.integrate_m_mp1_sub(fu_sub[m,:,:], m)
    return u0
    
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
      raise
    res = np.linalg.norm(u[0,:] - u0 - self.I_p_pp1[m,0,:], np.inf)
    for p in range(1,self.coll.P):
      res = max(res, np.linalg.norm(u[p,:] - u[p-1,:] - self.I_p_pp1[m,p,:], np.inf))
    return res

  '''
  '''
  def predict(self, u0, u, usub, fu, fu_sub):
    try:
      u      = np.reshape(u, (self.coll.M, self.prob.dim))
      usub   = np.reshape(usub, (self.coll.M, self.coll.P, self.prob.dim))
      fu     = np.reshape(fu, (self.coll.M, self.prob.dim))
      fu_sub = np.reshape(fu_sub, (self.coll.M, self.coll.P, self.prob.dim))
      u0     = np.reshape(u0, (self.prob.dim,))
    except:
      raise

    u0_step = u0

    for m in range(self.coll.M):
      
        # standard step
        u[m,:] = self.prob.solve_f1(self.coll.coll.delta_m[m], u0_step)
        
        # embedded steps
        t = self.coll.coll_sub[m].tleft
        usub[m,0,:] = u0_step + self.coll.coll_sub[m].delta_m[0]*( self.prob.f1(u[m,:]) + self.prob.f2(u0_step, t) )
        fu_sub[m,0,:] = self.prob.f2(usub[m,0,:], self.coll.coll_sub[m].nodes[0])
        
        for p in range(1,self.coll.P):
          t = self.coll.coll_sub[m].nodes[p-1]
          usub[m,p,:] = usub[m,p-1,:] + self.coll.coll_sub[m].delta_m[p]*( self.prob.f1(u[m,:]) + self.prob.f2(usub[m,p-1,:], t) )
          fu_sub[m,p,:] = self.prob.f2(usub[m,p,:], self.coll.coll_sub[m].nodes[p])
        
        # overwrite standard value
        u0_step = usub[m,self.coll.P-1,:]
        u[m,:]  = usub[m,self.coll.P-1,:]
        fu[m,:] = self.prob.f1(u[m,:])

    return 0


  '''
  '''
  def sweep(self, u0, u, usub, fu, fu_sub, fu_, fu_sub_):
    try:
      u     = np.reshape(u, (self.coll.M,self.prob.dim))
      usub  = np.reshape(usub, (self.coll.M, self.coll.P, self.prob.dim))
      fu    = np.reshape(fu, (self.coll.M,self.prob.dim))
      fu_   = np.reshape(fu_, (self.coll.M,self.prob.dim))
      fu_sub  = np.reshape(fu_sub, (self.coll.M, self.coll.P, self.prob.dim))
      fu_sub_  = np.reshape(fu_sub_, (self.coll.M, self.coll.P, self.prob.dim))
    
    except:
      raise
    
    
    # update integral terms
    self.update_I_m_mp1(fu_, fu_sub_)
    self.update_I_p_pp1(fu_, fu_sub_)
    
    for m in range(self.coll.M):
      
      if m==0:
        u_mm1 = u0
      
      # standard step
      rhs  = u_mm1 - self.coll.coll.delta_m[m]*( fu_[m,:] ) + self.I_m_mp1[m,:]
      
      fu_star = self.prob.f1(self.prob.solve_f1(self.coll.coll.delta_m[m], rhs))
      
      # embedded steps
      for p in range(self.coll.P):
      
        if p==0:
          t = self.coll.coll_sub[m].tleft
          if m==0:
            usub_mm1 = u0
          else:
            usub_mm1 = u[m-1,:]
          usub[m,p,:]   = usub_mm1 + self.coll.coll_sub[m].delta_m[p]*( fu_star - fu_[m,:]  ) + self.I_p_pp1[m,p,:]
        else:
          t = self.coll.coll_sub[m].nodes[p-1]
          usub_mm1 = usub[m,p-1,:]
          usub[m,p,:]   = usub_mm1 + self.coll.coll_sub[m].delta_m[p]*( fu_star - fu_[m,:] + fu_sub[m,p-1,:] - fu_sub_[m,p-1,:] ) + self.I_p_pp1[m,p,:]
        
        fu_sub[m,p,:] = self.prob.f2(usub[m,p,:], self.coll.coll_sub[m].nodes[p])
        
      # Prepare for next standard time step
      u[m,:]  = usub[m,self.coll.P-1,:]
      u_mm1   = u[m,:]
      fu[m,:] = self.prob.f1(u[m,:])

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
    return np.reshape( np.linalg.inv(Mcoll).dot(u0*np.ones(self.coll.M)), (self.coll.M, 1))

  def get_collocation_solution_sub(self, u0, m):  
    Q = self.coll.coll_sub[m].Qmat
    Q = Q[1:,1:]
    try:
      Mcoll = np.eye(self.coll.P) - Q*self.prob.lambda_2
    except:
      raise
    assert self.prob.lambda_1==0.0, "Can only compute embedded collocation solution analytically if lambda_1=0"
    return np.reshape( np.linalg.inv(Mcoll).dot(u0*np.ones(self.coll.P)), (self.coll.P, 1))
