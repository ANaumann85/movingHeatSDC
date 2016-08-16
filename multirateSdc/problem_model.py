import numpy as np

class problem_model():

  def __init__(self, a, nu):
    self.a = a
    self.nu = nu
    self.f  = np.array([[nu, 0.0], [0.0, nu]])
    self.g  = np.array([[a, -a],[-a, a]])
    self.dim = 2
    
  def f1(self, u):
    try:
      u = np.reshape(u, ((2,)))
    except:
      raise
    return self.f.dot(u)

  def solve_f1(self, alpha, b):
    try:
      b = np.reshape(b, ((2,)))
    except:
      raise
    M = np.eye(2) - alpha*self.f
    return np.linalg.solve(M, b)

  def f2(self, u):
    try:
      u = np.reshape(u, ((2,)))
    except:
      raise
    return self.g.dot(u)