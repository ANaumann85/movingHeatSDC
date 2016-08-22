from firedrake_problem import problem as fdProb
from problem import problem as scalProb
from sdc import sdc_step
import numpy as np
import copy
import firedrake as fd

#lambda_1 = -0.1
#lambda_2 = -1.0
#prob = scalProb(lambda_1, lambda_2)

nu= 1.0e-3 #1.0e-3
alpha=1.0e-4
nx=10
ny=40
prob = fdProb(nu, alpha, nx, ny)

tstart = 0.0
tend   = 1000.0
nsteps = 20000

dt = (tend - tstart)/float(nsteps)
u0    = prob.getU0()
#u0    = 2.0 
#u_ex  = u0*np.exp(tend*(lambda_1+lambda_2))

prob.startFile("T_noSrc_imex_%d" % nsteps)
prob.write(u0)
for n in range(nsteps):
  tstart = float(n)*dt
  tend   = float(n+1)*dt
  t=n*dt
  r=u0+dt*prob.f2(u0) #, t+dt)
  u0=prob.solve_f1(dt,r)
prob.write(u0,tend)
