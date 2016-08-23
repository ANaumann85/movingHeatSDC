from problem_firedrake import problem as fdProb
from problem import problem as scalProb
from sdc import sdc_step
import numpy as np
import copy
import firedrake as fd
import sys

if len(sys.argv) < 2:
    print "usage:" << sys.argv[0], " <nStep>"
    sys.exit(1)

nu= 1.0e-3 #1.0e-3
alpha=1.0e-4
nx=10
ny=40
prob = fdProb(nu, alpha, nx, ny)

tstart = 0.0
tend   = 20.0
nsteps = int(sys.argv[1]) 

dt = (tend - tstart)/float(nsteps)
u0    = prob.getU0()
gamma = 1-np.sqrt(1/2.0)
#u0    = 2.0 
#u_ex  = u0*np.exp(tend*(lambda_1+lambda_2))

prob.startFile("M_3_P_2_testMove/T_ros2_%d" % nsteps)
prob.write(u0)
for n in range(nsteps):
  tstart = float(n)*dt
  tend   = float(n+1)*dt
  t=n*dt
  r=prob.f1(u0)+prob.f2(u0, tstart)
  k1=prob.solve_f1(dt*gamma,r)
  r=prob.f1(u0+dt*k1)+prob.f2(u0+dt*k1, tend)-2*k1
  k2=prob.solve_f1(dt*gamma, r)
  u0 += (3*k1+k2)*dt/2
prob.write(u0,tend)
