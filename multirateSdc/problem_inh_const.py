import numpy as np

class problem_inh_const():

  def __init__(self, nu, b):
    self.nu   = nu
    self.dim = 1
    self.b=b
    
  def f1(self, u):
    return self.nu*u

  def solve_f1(self, a, b):
    return b/(1.0 - self.nu*a)

  def constSrc(self):
    return self.b

  def f2(self, u_, t):
    return self.b
