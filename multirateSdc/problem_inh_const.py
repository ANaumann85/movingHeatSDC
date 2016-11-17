import numpy as np

class problem_inh_const():

  def __init__(self, nu1, nu2, b):
    self.nu1 = nu1
    self.nu2 = nu2
    #self.nu1=nu1+nu2
    #self.nu2=0
    self.dim = 1
    self.b=b
    
  def f1(self, u):
    return self.nu1*u

  def solve_f1(self, a, b):
    return b/(1.0 - self.nu1*a)

  def constSrc(self):
    return self.b

  def f2(self, u_, t):
    return self.nu2*u_+self.b
