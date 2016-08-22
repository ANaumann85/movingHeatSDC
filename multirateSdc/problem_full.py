import numpy as np
import scipy.linalg as spla
import scipy.sparse as sp

class problem_full():

  '''
  '''
  def __init__(self, a, nu, alpha, v0, xaxis):
    self.a     = a
    self.nu    = nu
    self.alpha = alpha
    self.v0    = v0
    self.xaxis = xaxis
    self.dim   = np.size(xaxis)
    self.dx    = xaxis[1] - xaxis[0]
    e     = np.zeros(self.dim)
    e[0]  = -2.0
    e[1]  = 1.0
    e[-1] = 1.0
    self.S  = sp.csc_matrix(spla.circulant(e))
    self.S *= (self.nu/self.dx**2)
  
  '''
  '''
  def f1(self, u):
    try:
      u = np.reshape(u, ((self.dim,)))
    except:
      raise
    return self.S.dot(u)

  '''
  '''
  def solve_f1(self, alpha, b):
    try:
      b = np.reshape(b, ((self.dim,)))
    except:
      raise
    M = np.eye(self.dim) - alpha*self.S
    return np.linalg.solve(M, b)

  '''
  '''
  def f2(self, u, t):
    try:
      u = np.reshape(u, ((self.dim,)))
    except:
      raise
    b = self.get_hat(t)
    return np.multiply( self.alpha*(self.v0-u), b)

  '''
  '''
  def get_hat_trunc(self, t):
    u = 1.0/(2.0*np.pi) + 0.0*self.xaxis
    for k in range(1,3):
      u+= 1.0/np.pi*( np.cos(float(k)*self.a*t)*np.cos(float(k)*self.xaxis) + np.sin(float(k)*self.a*t)*np.sin(float(k)*self.xaxis) )
    return u

  '''
  '''
  def get_hat(self, t):
    x_c = self.a*t % self.xaxis[-1]
    u = np.zeros(self.dim)
    ind = np.argmin(np.abs(self.xaxis-x_c))
    u[ind] = 1.0/self.dx
    return u