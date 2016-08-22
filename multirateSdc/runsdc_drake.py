from firedrake_problem import problem as fdProb
from problem import problem as scalProb
from sdc import sdc_step
import numpy as np
import copy
import firedrake as fd
import sys
import os
from subprocess import call 

#lambda_1 = -0.1
#lambda_2 = -1.0
#prob = scalProb(lambda_1, lambda_2)

nu= 1.0e-3
alpha=1.0e-4
nx=10
ny=40
prob = fdProb(nu, alpha, nx, ny)

M = 3
P = 2
tstart = 0.0
tend   = 1000.0
nsteps = int(sys.argv[1]) 
dt = (tend - tstart)/float(nsteps)

K_iter = int(sys.argv[2])
u_    = np.zeros((M,prob.dim))
usub_ = np.zeros((M,P,prob.dim))
u     = np.zeros((M,prob.dim))
usub  = np.zeros((M,P,prob.dim))
u0    = prob.getU0()
#u0    = 2.0 
#u_ex  = u0*np.exp(tend*(lambda_1+lambda_2))

baseName="T_sdc_K_%d_%d" % (K_iter, nsteps)
prob.startFile(baseName)
prob.write(u0)
for n in range(nsteps):
  tstart = float(n)*dt
  tend   = float(n+1)*dt
  sdc    = sdc_step(M, P, tstart, tend, prob)
  
  # reset buffers to zero
  u     = np.zeros((M,prob.dim))
  usub  = np.zeros((M,P,prob.dim))
  u_    = np.zeros((M,prob.dim))
  usub_ = np.zeros((M,P,prob.dim))

  sdc.predict(u0, u_, usub_)
  for k in range(K_iter):
    sdc.sweep(u0, u, usub, u_, usub_)
    u_    = copy.deepcopy(u)
    usub_ = copy.deepcopy(usub)

  print ("+++ step: %3i +++ " % n)
  print ("standard residual: %5.3e " % sdc.residual(u0, u))
  print ("embedded residual: %5.3e " % sdc.sub_residual(u0, usub))
  u0 = u[M-1]
  prob.write(u0, tend)

#backup result at end
lastName=baseName+'_fixed_%d.vtu' % nsteps
os.rename(lastName, ('M_%d_P_%d/'+lastName)%(M, P))
#remove all other files beginning with basename
call('rm '+baseName+'*.vtu', shell=True)
call('rm '+baseName+'*.pvd', shell=True)
###
#print ("error:             %5.3e " % (abs(u0 - u_ex)/abs(u_ex)))

