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
    t=1.0
    uB = Function(self.mh.V_B, name="uB") 
    uB.vector().set_local(u)
    fslow=self.mh.solveM(self.mh.evalSlow(uB, t))
    return fslow.vector().get_local()

  def solve_f1(self, a, b):
    t=1.0
    bM = Function(self.mh.V_B, name="uB")
    bM.vector().set_local(b)
    uB = Function(self.mh.V_B, name="uB") 
    uB = self.mh.evalMmJ(uB, t, self.mh.evalM(bM), a)
    return uB.vector().get_local()

  def f2(self, u):
    t=1.0
    uB = Function(self.mh.V_B, name="uB") 
    uB.vector().set_local(u)
    uB = self.mh.evalFast(uB, t)
    uB = self.mh.solveM(uB)
    return uB.vector().get_local()

  def getU0(self):
    uRet = Function(self.mh.V_B, name="u0")
    uRet.interpolate(Expression("0.0*x[0]*(x[0]-1)*x[1]*(x[1]-4)")) 
    return uRet.vector().get_local()

  def startFile(self, fname):
    self.outFixed=File(fname+"_fixed.pvd")
    self.outMoving=File(fname+"_moving.pvd")

  def write(self, uVec, time=0):
    uFixed = Function(self.mh.V_B, name="uB")
    uFixed.vector().set_local(uVec)
    self.outFixed.write(uFixed, time=time)
    if hasattr(self.mh,'u_A'):
      self.outMoving.write(self.mh.u_A, time=time)
