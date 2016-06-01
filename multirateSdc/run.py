from problem import *
import numpy as np
import copy

def imex_step(u0, dt, u0_, u_, b, prob):
  b += u0 + dt*( prob.fexpl(u0) - prob.fexpl(u0_) ) - dt*prob.fimpl(u_) 
  return prob.solve(dt, b)

M = 3 # SDC nodes for slow process
P = 2 # SDC nodes for fast process
K = 3 # Number of sweeps

tend = 1.0
prob = problem(M,P, 0.0, tend)
#prob.print_nodes()

u      = np.zeros((M,1))
u_     = np.zeros((M,1))
usub   = np.zeros((M,P))
u0     = 1.0
uex    = np.exp((prob.lambda_fast + prob.lambda_slow)*tend)

# Prediction sweep to populate u
u[0]   = imex_step(u0, prob.coll.delta_m[0], 0.0, 0.0, 0.0, prob)
for m in range(1,M):
  u[m] = imex_step(u[m-1], prob.coll.delta_m[m], 0.0, 0.0, 0.0, prob)
 
print abs(prob.end_value(u,1.0) - uex)/abs(uex)
#print abs(u[M-1,0] - uex)/abs(uex)

for k in range(0,K):

  u_ = copy.deepcopy(u)
  prob.update_I_m(u_)

  ### update k -> k+1 ###
  b    = prob.I_m_mp1[0]
  u[0] = imex_step(u0, prob.coll.delta_m[0], u0, u_[0], b, prob)
  for m in range(1,M):
    b    = prob.I_m_mp1[m]
    u[m] = imex_step(u[m-1], prob.coll.delta_m[m], u_[m-1], u_[m], b, prob)


  print ("**** Sweep k = %2i ****" % k)
  print ("End value error: %5.3e" % (abs( prob.end_value(u,1.0) - uex )/abs(uex)) )
  print ("Sub-step end value error: %5.3e" % (abs(prob.end_value_sub(usub, 1.0) -  uex)/abs(uex)) )
  print "************************"
