from problem import problem
from sdc import sdc_step
import numpy as np
import copy

M = 8
P = 5
tstart = 0.0
tend   = 2.5
lambda_1 = -0.1
lambda_2 = -1.0
prob = problem(lambda_1, lambda_2)
sdc  = sdc_step(M, P, tstart, tend, prob)
K_iter = 15
u_    = np.zeros((M,1))
usub_ = np.zeros((M,P))
u     = np.zeros((M,1))
usub  = np.zeros((M,P))
u0    = 1.0
u_ex  = u0*np.exp(tend*(lambda_1+lambda_2))

### PREDICTOR: populate u_ and usub_
u_, usub_ = sdc.predict(u0)

### SDC iteration
for k in range(K_iter):
  
  # print
  print ("========== iteration %2i  ========== " % k)
  print ("standard residual: %5.3e " % sdc.residual(u0, u_))
  print ("embedded residual: %5.3e " % sdc.sub_residual(u0, usub_))
  print ("final error:       %5.3e" % abs(u_[M-1] - u_ex))

  sdc.sweep(u0, u, usub, u_, usub_)

  # print magnitude of updates
  print ("standard update:   %5.3e " % np.linalg.norm(u - u_, np.inf))
  print ("embedded update:   %5.3e " % np.linalg.norm(usub - usub_, np.inf))
  print ""

  # prepare next iteration
  u_ = copy.deepcopy(u)
  usub_ = copy.deepcopy(usub)

