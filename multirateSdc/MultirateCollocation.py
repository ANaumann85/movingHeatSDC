from Collocation import *
from CollocationClasses import *

class multirateCollocation(object):

  '''
  '''
  def __init__(self, M, P, tleft, tright):
    # approximations of integrals from t_m to t_mp1
    self.I_m_mp1 = np.zeros((M,1))
    # approximations of integrals from t_p to t_pp1
    self.I_p_pp1 = np.zeros((M,P))
    assert tleft<tright, "tleft must be smaller than tright"
    self.dt      = abs(tright - tleft)
    self.M       = M
    self.P       = P

    self.coll    = CollGaussRadau_Right(num_nodes = M, tleft = tleft, tright = tright)
    self.coll_sub = []
    self.coll_sub.append(CollGaussRadau_Right(num_nodes=P, tleft = tleft, tright = self.coll.nodes[0]))
    for i in range(1,M):
      self.coll_sub.append(CollGaussRadau_Right(num_nodes=P, tleft = self.coll.nodes[i-1] , tright=self.coll.nodes[i]))
    # NOTE: this assumes that the last node if self.coll is equal to tright; otherwise need one additional sub-step collocation class

    self.S_mnp = np.zeros((M,M,P))
    for n in range(M):
      fl = np.zeros(M)
      fl[n] = 1.0
      coeff = self.coll._poly_newton(fl)
      for m in range(M):
        for p in range(P):
          if p==0:
            t_m_p = self.coll_sub[m].tleft
          else:
            t_m_p = self.coll_sub[m].nodes[p-1]
          t_m_pp1 = self.coll_sub[m].nodes[p]
          [nodes, weights] = CollBase._GaussLegendre(M, t_m_p, t_m_pp1)
          flag    = self.coll._evaluate_horner(nodes, coeff)
          self.S_mnp[m,n,p] = CollBase.evaluate(weights, flag)

  '''
  '''
  def integrate_p_pp1(self, fu, m, p):
    try:
      fu = np.reshape(fu, (self.M,1))
    except:
      raise TypeError("Failed to convert argument fu into shape Mx1")
    intvalue = 0.0
    for n in range(self.M):
      intvalue += self.S_mnp[m,n,p]*fu[n]
    return intvalue

  '''
  '''
  def integrate_p_pp1_sub(self, fu, m, p):
    try:
      fu = np.reshape(fu, (self.P,1))
    except:
      raise TypeError("Failed to convert argument fu into shape Px1")
    Smat = self.coll_sub[m].Smat
    Smat = Smat[1:,1:]
    intvalue = 0.0
    for j in range(self.P):
      intvalue += Smat[p,j]*fu[j]
    return intvalue
      

  '''
  '''
  def integrate_m_mp1(self, fu, m):
    try:
      fu = np.reshape(fu, (self.M,1))
    except:
      raise TypeError("Failed to convert argument fu into shape Mx1")
