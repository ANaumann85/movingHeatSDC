class problem():

  def __init__(self, lambda_1, lambda_2):
    self.lambda_1 = lambda_1
    self.lambda_2 = lambda_2

  def f1(self, u):
    return self.lambda_1*u

  def solve_f1(self, a, b):
    return b/(1.0 - self.lambda_1*a)

  def f2(self, u):
    return self.lambda_2*u
