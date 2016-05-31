from problem import *
import numpy as np

M = 4 # SDC nodes for slow process
P = 3 # SDC nodes for fast process
K = 4 # Number of sweeps

prob = problem(M,P)
#prob.print_nodes()

u      = np.zeros((M,1))
u_     = np.zeros((M,1))
u0     = 1.0

# Prediction sweep to populate u
u[0]   = prob.solve(prob.coll.delta_m[0], u0)
for m in range(1,M):
  b      = u[m-1] + prob.coll.delta_m[m]*prob.fexpl(u[m-1])
  u[m]   = prob.solve(prob.coll.delta_m[m], b)

print abs(prob.end_value(u,1.0) - np.exp((prob.lambda_fast + prob.lambda_slow)))

for k in range(0,K):
  prob.update_I_m(u)
  ### update k -> k+1 ###
  b    = u0 + prob.I_m_mp1[0] - prob.coll.delta_m[0]*prob.fimpl(u[0]) # explicit terms cancel out here
  u[0] = prob.solve(prob.coll.delta_m[0], b)
  for m in range(1,M):
    b    = u[m-1] + prob.I_m_mp1[m] - prob.coll.delta_m[m]*prob.fimpl(u[m]) + prob.coll.delta_m[m]*( prob.fexpl(u[m-1]) - prob.fexpl(u[m-1]) )
    u[m] = prob.solve(prob.coll.delta_m[m], b)
  print abs(prob.end_value(u,1.0) - np.exp((prob.lambda_fast + prob.lambda_slow)))

