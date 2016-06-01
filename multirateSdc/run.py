from problem import *
import numpy as np
import copy

def imex_step(u0, dt, u0_, u_, b, prob):
  b += u0 + dt*( prob.fexpl(u0) - prob.fexpl(u0_) ) - dt*prob.fimpl(u_) 
  return prob.solve(dt, b)

def subcycling(u0, usub, usub_, prob, m):
  b         = prob.I_p_pp1[m,0]
  usub[m,0] = imex_step(u0, prob.coll_fast[m].delta_m[0], u0, usub_[m,0], b, prob)
  for p in range(1,prob.P):
    b         = prob.I_p_pp1[m,p]
    usub[m,p] = imex_step(usub[m,p-1], prob.coll_fast[m].delta_m[p], usub_[m,p-1], usub_[m,p], b, prob)
  return usub

M = 3 # SDC nodes for slow process
P = 3 # SDC nodes for fast process
K = 15 # Number of sweeps

tend = 1.0
prob = problem(M,P, 0.0, tend)
#prob.print_nodes()

u      = np.zeros((M,1))
u_     = np.zeros((M,1))
usub   = np.zeros((M,P))
usub_  = np.zeros((M,P))
u0     = 1.0
uex    = np.exp((prob.lambda_fast + prob.lambda_slow)*tend)

# Prediction sweep to populate u

print prob.get_residual(u, u0)
print prob.get_residual_sub(usub, u0)

for k in range(0,K):

  u_ = copy.deepcopy(u)
  usub_ = copy.deepcopy(usub)
  prob.update_I_m(u_)
  prob.update_I_p(usub_)

  ### update k -> k+1 ###
  b    = prob.I_m_mp1[0]
  u[0] = imex_step(u0, prob.coll.delta_m[0], u0, u_[0], b, prob)
  usub = subcycling(u0, usub, usub_, prob, 0)
  for m in range(1,M):
    b    = prob.I_m_mp1[m]
    u[m] = imex_step(u[m-1], prob.coll.delta_m[m], u_[m-1], u_[m], b, prob)
    usub0 = prob.end_value_sub(usub, u0, m-1)
    usub = subcycling(usub0, usub, usub_, prob, m)
  #print prob.get_residual(u, u0)
  #print prob.get_residual_sub(usub, u0)
  #print ""
  print ("**** Sweep k = %2i ****" % k)
  print ("End value error: %5.3e" % (abs( prob.end_value(u,1.0) - uex )/abs(uex)) )
  print ("Residual: %5.3e" % prob.get_residual(u, u0))
  print ("Sub-step end value error: %5.3e" % (abs(prob.end_value_sub_all(usub, 1.0) -  uex)/abs(uex)) )
  print ("Sub-step residual: %5.3e" % prob.get_residual_sub(usub, u0))
  print "************************"
