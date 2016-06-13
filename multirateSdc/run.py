from problem import *
import numpy as np
import copy

'''
'''
def imex_step(u0, dt, u0_, u_, b, prob):
  b += u0 + dt*( prob.fexpl(u0) - prob.fexpl(u0_) ) - dt*prob.fimpl(u_) 
  return prob.solve(dt, b)

'''
'''
def impeuler_step(u0, dt, u0_, u_, b, prob):
  b += u0 - dt*prob.fimpl(u_)
  return prob.solve(dt, b)

'''
'''
def expeuler_step(u0, dt, u0_, u_, b, prob, rhs):
  b+= u0 + dt*(prob.fexpl(u0) - prob.fexpl(u0_) + rhs) 
  return b

'''
'''
def subcycling_imex(u0, usub, usub_, prob, m):
  b         = prob.I_p_pp1[m,0]
  usub[m,0] = imex_step(u0, prob.coll_fast[m].delta_m[0], u0, usub_[m,0], b, prob)
  for p in range(1,prob.P):
    b         = prob.I_p_pp1[m,p]
    usub[m,p] = imex_step(usub[m,p-1], prob.coll_fast[m].delta_m[p], usub_[m,p-1], usub_[m,p], b, prob)
  return usub

'''
'''
def subcycling_mr(u0, usub, usub_, prob, m, rhs_slow):
  b = prob.I_p_pp1[m,0]
  # explicit Euler steps in fast process plus additional term from slow process as rhs_slow
  usub[m,0] = expeuler_step(u0, prob.coll_fast[m].delta_m[0], u0, usub_[m,0], b, prob, rhs_slow)
  for p in range(1,prob.P):
    b = prob.I_p_pp1[m,p]
    usub[m,p] = expeuler_step(usub[m,p-1], prob.coll_fast[m].delta_m[p], usub_[m,p-1], usub_[m,p], b, prob, rhs_slow)
  return usub

M = 4 # SDC nodes for slow process
P = 4 # SDC nodes for fast process
K = 5 # Number of sweeps

tend     = 0.5
prob     = problem(M, P, 0.0, tend, -0.25, -1.0)


#prob.print_nodes()

u       = np.zeros((M,1))
u_      = np.zeros((M,1))
usub    = np.zeros((M,P))
usub_   = np.zeros((M,P))
u0      = 1.0
uex     = np.exp((prob.lamb)*tend)
uex_sub = np.exp((prob.lamb)*tend)

### test collocation and sub-step collocation method for consistency ###
coll = prob.get_coll_solution(u0)
collsub = prob.get_coll_solution_sub(u0)
coll_rec = np.zeros((M,1))
for m in range(M):
  coll_rec[m,0] = prob.end_value_sub(collsub, u0, m)
print ("Defect between collocation and sub-step collocation solution: %5.3e" % np.linalg.norm(coll - coll_rec, np.inf))

#####

# initial step
u[0]     = impeuler_step(u0, prob.coll.delta_m[0], u0, u_[0], 0.0, prob)
rhs_slow = prob.fimpl(u[0])

# Sub-cycling in first sub-step
usub = subcycling_mr(u0, usub, usub_, prob, 0, rhs_slow)

# prediction loop
for m in range(1,M):
  u[m]     = impeuler_step(u[m-1], prob.coll.delta_m[m], u_[m-1], u_[m], 0.0, prob)
  usub0    = prob.end_value_sub(usub, u0, m-1)
  rhs_slow = prob.fimpl(u[m])
  usub     = subcycling_mr(usub0, usub, usub_, prob, m, rhs_slow)

#print prob.get_residual(u, u0)
#print prob.get_residual_sub(usub, u0)

for k in range(0,K):

  print ("**** Sweep k = %2i ****" % k)
  print ("Update: %5.3e" % np.linalg.norm(u - u_, np.inf))
  print ("Sub-step update: %5.3e" % np.linalg.norm(usub - usub_, np.inf))
  print ""
  u_    = copy.deepcopy(u)
  usub_ = copy.deepcopy(usub)

  prob.update_I_m(u_)
  prob.update_I_p(usub_)

  ### update k -> k+1 ###
  b        = prob.I_m_mp1[0]
  u[0]     = impeuler_step(u0, prob.coll.delta_m[0], u0, u_[0], b, prob)
  rhs_slow = prob.fimpl(u[0]) - prob.fimpl(u_[0])
  usub = subcycling_mr(u0, usub, usub_, prob, 0, rhs_slow)

  for m in range(1,M):
    b        = prob.I_m_mp1[m]
    u[m]     = impeuler_step(u[m-1], prob.coll.delta_m[m], u_[m-1], u_[m], b, prob)
    rhs_slow = prob.fimpl(u[m]) - prob.fimpl(u_[m])
    usub0    = prob.end_value_sub(usub, u0, m-1)
    usub     = subcycling_mr(u0, usub, usub_, prob, m, rhs_slow)
  
  # reconstruct coarse level values from usub can compare
  urec = np.zeros((M,1))
  for m in range(M):
    urec[m,0] = prob.end_value_sub(usub, u0, m)
  print ("Difference u and reconstructed u: %5.3e" % np.linalg.norm(u - urec, np.inf))
  print ""

  print ("End value error: %5.3e" % (abs( prob.end_value(u,1.0) - uex )/abs(uex)) )
  print ("Residual: %5.3e" % prob.get_residual(u, u0))
  print ("Sub-step end value error: %5.3e" % (abs(prob.end_value_sub_all(usub, 1.0) -  uex)/abs(uex)) )
  print ("Sub-step residual: %5.3e" % prob.get_residual_sub(usub, u0))
  print "************************"
