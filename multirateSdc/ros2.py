import numpy as np
import copy

class ros2_step():

  '''
  '''
  def __init__(self, tstart, tend, problem):
    self.dt     = tend-tstart
    self.gamma  = 1.0 - np.sqrt(1.0/2.0)
    self.prob   = problem
    self.tstart = tstart
    self.tend   = tend

  def step(self, u0):
    r  = self.prob.f1(u0) + self.prob.f2(u0,self.tstart)
    k1 = self.prob.solve_f1(self.dt*self.gamma, r)
    r  = self.prob.f1(u0 + self.dt*k1) + self.prob.f2(u0 + self.dt*k1, self.tend) - 2.0*k1
    k2 = self.prob.solve_f1(self.dt*self.gamma, r)
    u0 += 0.5*(3.0*k1 + k2)*self.dt
    return u0