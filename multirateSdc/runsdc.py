from problem import problem
from sdc import sdc_step
import numpy as np
import copy

M = 8
P = 5
tstart = 0.0
tend   = 3.5
lambda_1 = -0.1
lambda_2 = -1.0
prob = problem(lambda_1, lambda_2)
sdc  = sdc_step(M, P, tstart, tend, prob)
K_iter = 15
u_    = np.zeros((M,1))
usub_ = np.zeros((M,P))
u     = u_
usub  = usub_
u0    = 1.0
u_ex  = u0*np.exp(tend*(lambda_1+lambda_2))

### PREDICTOR: populate u_ and usub_
u0_step = u0

for m in range(M):

  # standard step
  u_[m] = prob.solve_f1(sdc.coll.coll.delta_m[m], u0_step)

  # embedded steps
  usub_[m,0] = u0_step + sdc.coll.coll_sub[m].delta_m[0]*( prob.f1(u_[m]) + prob.f2(u0) )
  for p in range(1,P):
    usub_[m,p] = usub_[m,p-1] + sdc.coll.coll_sub[m].delta_m[p]*( prob.f1(u_[m]) + prob.f2(usub_[m,p-1]) )

  # overwrite standard value 
  u0_step = usub_[m,P-1]
  u_[m]    = usub_[m,P-1] # ???

### SDC iteration
for k in range(K_iter):
  
  # print
  print ("========== iteration %2i  ========== " % k)
  print ("standard residual: %5.3e " % sdc.residual(u0, u_))
  print ("embedded residual: %5.3e " % sdc.sub_residual(u0, usub_))
  print ("final error:       %5.3e" % abs(u_[M-1] - u_ex))

  u0_step = u0  

  for m in range(M):

    # standard step
    rhs  = u0_step - sdc.coll.coll.delta_m[m]*( prob.f1(u_[m]) ) + sdc.I_m_mp1[m]
    u[m] = prob.solve_f1(sdc.coll.coll.delta_m[m], rhs)
  
    # embedded steps
    usub[m,0] = u0_step + sdc.coll.coll_sub[m].delta_m[0]*( prob.f1(u[m]) - prob.f1(u_[m]) + prob.f2(u0_step) - prob.f2(u0_step) ) + sdc.I_p_pp1[m,0]
    for p in range(1,P):
      usub[m,p] = usub[m,p-1] + sdc.coll.coll_sub[m].delta_m[p]*( prob.f1(u[m]) - prob.f1(u_[m]) + prob.f2(usub[m,p-1]) - prob.f2(usub_[m,p-1]) ) + sdc.I_p_pp1[m,p]

    # overwrite standard value
    u0_step = usub[m,P-1]
    u[m]    = usub[m,P-1]

  # print magnitude of updates
  print ("standard update:   %5.3e " % np.linalg.norm(u - u_, np.inf))
  print ("emnedded update:   %5.3e " % np.linalg.norm(usub - usub_, np.inf))
  print ""
  # prepare next iteration
  u_ = copy.deepcopy(u)
  usub_ = copy.deepcopy(usub)
  sdc.update_I_m_mp1(u_, usub_)
  sdc.update_I_p_pp1(u_, usub_)
