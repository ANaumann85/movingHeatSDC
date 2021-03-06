import numpy as np

class problem_model():

  '''
  '''
  def __init__(self, a, nu, alpha, v0):
    self.a     = a
    self.nu    = nu
    self.alpha = alpha
    self.v0    = v0
    self.S     = -nu*np.diag([0.0, 1.0, 4.0, 1.0, 4.0])
    self.dim   = 5
  
  '''
  '''
  def get_mat1(self, t):
    A = np.zeros((5,5))
    # equation for a_0
    A[0,0] = 1.0/2.0
    
    # equation for a_1
    A[1,0] = np.cos(self.a*t)
    A[1,1] = np.cos(2.0*self.a*t)+1.0
    A[1,2] = np.cos(self.a*t)
    A[1,3] = np.sin(2.0*self.a*t)
    A[1,4] = np.sin(self.a*t)
    
    # equation for a_2
    A[2,0] = np.cos(2.0*self.a*t)
    A[2,1] = np.cos(self.a*t)
    A[2,2] = 1.0
    A[2,3] = -np.sin(self.a*t)
    
    # equation for b_1
    A[3,0] = np.sin(self.a*t)
    A[3,1] = np.sin(2.0*self.a*t)
    A[3,2] = -np.sin(self.a*t)
    A[3,3] = -np.cos(2.0*self.a*t)+1
    A[3,4] = np.cos(self.a*t)
    
    # equation for b_2
    A[4,0] = np.sin(2.0*self.a*t)
    A[4,1] = np.sin(self.a*t)
    A[4,3] = np.cos(self.a*t)
    A[4,4] = 1.0
    
    A *= -self.alpha/(2.0*np.pi)
    
    return A
  
  '''
  '''
  def get_b1(self, t):
    b = np.zeros((5,))
    b[0] = 1.0/2.0
    b[1] = np.cos(self.a*t)
    b[2] = np.cos(2.0*self.a*t)
    b[3] = np.sin(self.a*t)
    b[4] = np.sin(2.0*self.a*t)
    b *= self.alpha*self.v0/np.pi
    return b
  
  
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
  def solve_f1(self, c, b):
    try:
      b = np.reshape(b, ((self.dim,)))
    except:
      raise
    M = np.eye(self.dim) - c*self.S
    return np.linalg.solve(M, b)

  '''
  '''
  def f2(self, u, t):
    try:
      u = np.reshape(u, ((self.dim,)))
    except:
      raise
    A = self.get_mat1(t)
    b = self.get_b1(t)
    return A.dot(u) + b
    
  '''
  '''
  def f(self, u, t):
    return self.f1(u) + self.f2(u, t)

  '''
  '''
  def get_solution(self, u, x):
    uplot = 0.0*x + 0.5*u[0]
    uplot += u[1]*np.cos(1.0*x)
    uplot += u[3]*np.cos(2.0*x)
    uplot += u[2]*np.sin(1.0*x)
    uplot += u[4]*np.sin(2.0*x)
    return uplot