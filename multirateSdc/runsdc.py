from problem import problem
from sdc import sdc_step
import numpy as np
import copy

M = 4
P = 3
tstart = 0.0
tend   = 1.0
lambda_1 = -0.1
lambda_2 = -1.0
prob = problem(lambda_1, lambda_2)
sdc  = sdc_step(M, P, tstart, tend, prob)
K_iter = 1
u_    = np.zeros((M,1))
usub_ = np.zeros((M,P))
u     = u_
usub  = usub_
u0 = 1.0

### SDC iterations
for k in range(K_iter):

  ### PREDICTOR
  
  # first standard step
  u[0] = u0
  rhs  = u[0]
  u[1] = prob.solve_f1(sdc.coll.coll.delta_m[0], rhs)
  
  # embedded steps
  usub[0,0] = u0
  for p in range(1,P):
    usub[0,p] = usub[0,p-1] + sdc.coll.coll_sub[0].delta_m[p]*( prob.f1(u[1]) + prob.f2(usub[0,p-1]) )

  for m in range(1,M):
    rhs = u[m-1]
    u[m] = prob.solve_f1(sdc.coll.coll.delta_m[m], rhs)
    
    # embedded steps
    usub[m,0] = u[m-1]
    for p in range(1,P):
      usub[m,p] = usub[m,p-1] + sdc.coll.coll_sub[m].delta_m[p]*( prob.f1(u[m]) + prob.f2(usub[0,p-1]) )

    # Overwrite u[m] with end value from embedded step

