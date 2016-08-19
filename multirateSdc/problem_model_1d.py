class problem():

  def __init__(self, nu, omega):
    self.nu    = nu
    self.omega = omega
    self.dim   = 1
    
  def f1(self, u):
    return self.nu*u

  def solve_f1(self, a, b):
    return b/(1.0 - self.lambda_1*a)

  def f2(self, u):
    return self.lambda_2*u
