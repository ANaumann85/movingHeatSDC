import sys
from firedrake import *
import numpy as np

class problem():

  def __init__(self, nu, alpha, nx, ny):
    sys.path.append('../movingHeat')
    from heateq_moving3 import MovingHeat as mh

    self.mh = mh(nu, alpha, nx, ny)
    helper = Function(self.mh.V_B)
    self.dim = helper.vector().local_size()
    
  def f1(self, u):
    t=0.0
    uB = Function(self.mh.V_B, name="uB") 
    uB.vector().set_local(u)
    fslow=self.mh.solveM(self.mh.evalSlow(uB, t))
    return fslow.vector().get_local()

  def solve_f1(self, a, b):
    t=0.0
    bM = Function(self.mh.V_B, name="uB")
    bM.vector().set_local(b)
    uB = Function(self.mh.V_B, name="uB") 
    uB = self.mh.evalMmJ(uB, t, self.mh.evalM(bM), a)
    return uB.vector().get_local()

  def f2(self, u):
    t=0.0
    #uB = Function(self.mh.V_B, name="uB") 
    #uB.vector().set_local(u)
    #uB = self.mh.evalFast(uB, t)
    return np.array(range(0, self.dim))*0.0 #uB.vector().get_local()

  def getU0(self):
    uRet = Function(self.mh.V_B, name="u0")
    uRet.interpolate(Expression("x[0]*(x[0]-1)*x[1]*(x[1]-4)"))
    return uRet.vector().get_local()

  def startFile(self, fname):
      self.outFixed=File(fname+"_fixed.pvd")
      self.outMoving=File(fname+"_moving.pvd")

  def write(self, uVec):
      uFixed = Function(self.mh.V_B, name="uB")
      uFixed.vector().set_local(uVec)
      self.outFixed.write(uFixed)
