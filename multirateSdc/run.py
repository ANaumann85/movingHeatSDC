from problem import *
import numpy as np

M = 4 # SDC nodes for slow process
P = 3 # SDC nodes for fast process
K = 4 # Number of sweeps

tend = 1.25
prob = problem(M,P, 0.0, tend)
prob.print_nodes()

u      = np.zeros((M,1))
u_     = np.zeros((M,1))
u0     = 1.0

uf     = np.zeros((M,P))
# Prediction sweep to populate u and uf
# populate uf at first
uf[0,0] = u0+prob.coll_fast.delta_m[0]*prob.fexpl(u0)
for m in range(1, P):
    uf[0,m]=uf[0,m-1]+prob.coll_fast.delta_m[m]*prob.fexpl(uf[0,m-1])
u[0]   = prob.solve(prob.coll.delta_m[0], u0+prob.coll.delta_m[0]*prob.fexpl(uf[0,P-1])) #u0 <-> u(delta_m[0])
for m in range(1,M):
  b      = u[m-1] + prob.coll.delta_m[m]*prob.fexpl(uf[m,P-1])
  u[m]   = prob.solve(prob.coll.delta_m[m], b)

print abs(prob.end_value(u,1.0) - np.exp((prob.lambda_fast + prob.lambda_slow)*tend))

for k in range(0,K):
  prob.update_I_m(u)
  ### update k -> k+1 ###
  b    = u0 + prob.I_m_mp1[0] - prob.coll.delta_m[0]*prob.fimpl(u[0]) # explicit terms cancel out here
  u[0] = prob.solve(prob.coll.delta_m[0], b)
  for m in range(1,M):
    b    = u[m-1] + prob.I_m_mp1[m] - prob.coll.delta_m[m]*prob.fimpl(u[m]) + prob.coll.delta_m[m]*( prob.fexpl(u[m-1]) - prob.fexpl(u[m-1]) )
    u[m] = prob.solve(prob.coll.delta_m[m], b)
  print abs(prob.end_value(u,1.0) - np.exp((prob.lambda_fast + prob.lambda_slow)*tend))
