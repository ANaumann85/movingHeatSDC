import numpy as np

class problem_tcoeff():

  def __init__(self, nu, a):
    self.nu  = nu
    self.a   = a
    self.dim = 1
    
  def f1(self, u):
    return self.nu*u

  def solve_f1(self, a, b):
    return b/(1.0 - self.nu*a)

  def f2(self, u_, t):
    return self.a*np.cos(t)*u_