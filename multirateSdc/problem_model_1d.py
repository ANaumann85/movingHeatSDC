import numpy as np
class problem_model_1d():

  #periodic 1d heat with source al\delta(x-at)(v0-u(t,x)), i.e. 
  #\dot{u} = \Delta u + al\delta(x-at)(v0-u(t,x))
  #approximate \delta(x-at) as point or as triangle with unit area
  #stencil point:
  #\dot u_{0}=(u_{N}-2u_0+u_1)/h^2 +al\delta(x_0-at)(v_0-u_0)
  #\dot u_{j}=(u_{j-1}-2u_{j}+u_{j+1})/h^2+al\delta(x_j-at)(v_0-u_j), j=N-1..1
  #\dot u_{N}=(u_{N-1}-2u_N+u_1)/h^2+al\delta(x_N-at)(v_0-u_N)
  #stencil unit area:
  #\dot u_{0}=(u_{N}-2u_0+u_1)/h^2 +al*tri(x_0-at,h)(v_0-u_0)/h
  #\dot u_{j}=(u_{j-1}-2u_{j}+u_{j+1})/h^2+al*tri(x_j-at,h)(v_0-u_j)/h, j=N-1..1
  #\dot u_{N}=(u_{N-1}-2u_N+u_1)/h^2+al*tri(x_N-at,h)(v_0-u_N)/h
  #utilizing tri(x,h)=1-|x/h| 0<=x<h, 0 else
  def __init__(self, a, nu, alpha, v0):
    self.a     = a
    self.nu    = nu
    self.alpha = alpha
    self.v0    = v0
    self.dim   = 10
    N=self.dim
    self.S=np.diag(np.ones(N-1),-1)-2*np.diag(np.ones(N))+np.diag(np.ones(N-1),1)
    self.S[0,N-1]=1;
    self.S[N-1,0]=1;
    self.h=1.0/self.dim
    self.S=self.S*nu/(self.h*self.h)
    
  def get_mat1(self, t):
    jc=int(self.a*t/self.h)
    B=np.zeros((self.dim,self.dim))
    B[jc,jc] =-self.alpha/self.h
    return B

  def get_b1(self, t):
    jc=int(self.a*t/self.h)
    b=np.zeros((self.dim,))
    b[jc]=self.v0*self.alpha/self.h
    return b

  def f1(self, u):
    try:
      u = np.reshape(u, ((self.dim,)))
    except:
      raise
    return self.S.dot(u)

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
