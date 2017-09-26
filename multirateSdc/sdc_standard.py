import numpy as np
import copy
from Collocation import *
from CollocationClasses import *
from problem import problem

class sdc_standard_step():

  '''
  '''
  def __init__(self, M, tstart, tend, problem):
    self.coll    = CollGaussRadau_Right(M, tstart, tend)
    self.I_m_mp1 = np.zeros((M,problem.dim))
    self.prob    = problem
    self.dt      = tend-tstart

  '''
  '''
  def update_I_m_mp1(self, u):
    try:
      u     = np.reshape(u, (self.coll.num_nodes,self.prob.dim))
    except:
      raise
    
    fu = self.evaluate_f(u)
    
    for m in range(self.coll.num_nodes):
      self.I_m_mp1[m,:] = self.integrate_m_mp1(fu, m)

  def integrate_m_mp1(self, fu, m):
    #try:
    fu = np.reshape(fu, (self.coll.num_nodes,self.prob.dim))
    #except:
    #  raise TypeError("Failed to convert argument fu into shape Mx1")
    Smat = self.coll.Smat
    Smat = Smat[1:,1:]
    intvalue = np.zeros(self.prob.dim)
    for j in range(self.coll.num_nodes):
      intvalue += Smat[m,j]*fu[j,:]
    return intvalue

  '''
  '''
  def evaluate_f(self, u,):
    fu     = np.zeros((self.coll.num_nodes,self.prob.dim))
    for m in range(self.coll.num_nodes):
      fu[m,:] = self.prob.f1(u[m,:]) + self.prob.f2(u[m,:], self.coll.nodes[m])
    return fu
    
  '''
  '''
  def residual(self, u0, u):
    try:
      u = np.reshape(u, (self.coll.num_nodes,self.prob.dim))
    except:
      raise
    res = np.linalg.norm(u[0,:] - u0 - self.I_m_mp1[0,:], np.inf)
    for m in range(1,self.coll.num_nodes):
      resm = np.linalg.norm(u[m,:] - u[m-1,:] - self.I_m_mp1[m,:], np.inf)
      res  = max(res, resm)
    return res


  '''
  '''
  def predict(self, u0, u):
    try:
      u    = np.reshape(u, (self.coll.num_nodes, self.prob.dim))
      u0   = np.reshape(u0, (self.prob.dim,))
    except:
      raise
    
    for m in range(self.coll.num_nodes):
      
        # standard step (explicit part)
        if m==0:
          rhs = u0 + self.coll.delta_m[m]*self.prob.f2(u0, self.coll.tleft)
        else:
          rhs = u[m-1,:] + self.coll.delta_m[m]*self.prob.f2(u[m-1,:], self.coll.nodes[m-1])
        
        # standard step (implicit part)
        u[m,:] = self.prob.solve_f1(self.coll.delta_m[m], rhs)

    return u


  '''
  '''
  def sweep(self, u0, u, u_):
    try:
      u     = np.reshape(u, (self.coll.num_nodes,self.prob.dim))
      u_    = np.reshape(u_, (self.coll.num_nodes,self.prob.dim))

    except:
      raise
    
    # update integral terms
    self.update_I_m_mp1(u_)

    for m in range(self.coll.num_nodes):
      
      # standard step (explicit part)
      if m==0: # explicit term cancels out for first substep
        rhs = u0 - self.coll.delta_m[m]*( self.prob.f1(u_[m,:]) ) + self.I_m_mp1[m,:]
      else:
        rhs  = u[m-1,:] - self.coll.delta_m[m]*( self.prob.f1(u_[m,:]) ) \
                        + self.coll.delta_m[m]*( self.prob.f2(u[m-1,:], self.coll.nodes[m-1]) \
                              - self.prob.f2(u_[m-1,:], self.coll.nodes[m-1])) \
                        + self.I_m_mp1[m,:]

      # standard step (implicit part)
      u[m,:] = self.prob.solve_f1(self.coll.delta_m[m], rhs)

    return 0

  '''
  '''
  def get_collocation_solution(self, u0):  
    Q = self.coll.Qmat
    Q = Q[1:,1:]
    try:
      Mcoll = np.eye(self.coll.num_nodes) - Q*(self.prob.lambda_1 + self.prob.lambda_2)
    except:
      raise
    return np.reshape( np.linalg.inv(Mcoll).dot(u0*np.ones(self.coll.num_nodes)), (self.coll.num_nodes, 1))
