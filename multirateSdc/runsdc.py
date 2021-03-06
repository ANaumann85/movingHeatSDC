from problem import problem
from sdc import sdc_step
import numpy as np
import copy

M = 5
P = 4
tstart = 0.0
tend   = 10.0
nsteps = 4
dt = (tend - tstart)/float(nsteps)

lambda_1 = -0.1
lambda_2 = -1.0
prob = problem(lambda_1, lambda_2)

K_iter = 12
u_    = np.zeros((M,1,1))
usub_ = np.zeros((M,P,1))
u     = np.zeros((M,1,1))
usub  = np.zeros((M,P,1))
u0    = 2.0
u_ex  = u0*np.exp(tend*(lambda_1+lambda_2))

for n in range(nsteps):
  tstart = float(n)*dt
  tend   = float(n+1)*dt
  sdc    = sdc_step(M, P, tstart, tend, prob)
  
  # reset buffers to zero
  u      = np.zeros((M,1,1))
  usub   = np.zeros((M,P,1))
  u_     = np.zeros((M,1,1))
  usub_  = np.zeros((M,P,1))
  fu     = np.zeros((M, 1))
  fu_sub = np.zeros((M,P,1))

  sdc.predict(u0, u_, usub_, fu, fu_sub)
  for k in range(K_iter):
    sdc.sweep(u0, u, usub, fu, fu_sub)
    u_    = copy.deepcopy(u)
    usub_ = copy.deepcopy(usub)

  print ("+++ step: %3i +++ " % n)
  print ("standard residual: %5.3e " % sdc.residual(u0, u))
  print ("embedded residual: %5.3e " % sdc.sub_residual(u0, usub))
  u0 = u[M-1]

###
print ("error:             %5.3e " % (abs(u0 - u_ex)/abs(u_ex)))
