from problem import *
import numpy as np

M = 4 # SDC nodes for slow process
P = 3 # SDC nodes for fast process
K = 8 # Number of sweeps

tend = 1.0
prob = problem(M,P, 0.0, tend, -1.0, -2.0) #lslow, lfast
#prob.print_nodes()

u      = np.zeros((M,1))
u_     = np.zeros((M,1))
u0     = 1.0

uf     = np.zeros((M,P))
# Prediction sweep to populate u and uf
# populate uf at first
uf[0,0] = u0+prob.coll_fast[0].delta_m[0]*prob.fexpl(u0)
for m in range(1, P):
  uf[0,m]=uf[0,m-1]+prob.coll_fast[0].delta_m[m]*prob.fexpl(uf[0,m-1])
u[0]   = prob.solve(prob.coll.delta_m[0], u0+prob.coll.delta_m[0]*prob.fexpl(uf[0,P-1])) #u0 <-> u(delta_m[0])
for m in range(1,M):
  for p in range(1, P):
    uf[m,p]=uf[m,p-1]+prob.coll_fast[m].delta_m[p]*prob.fexpl(uf[m,p-1])
  b      = u[m-1] + prob.coll.delta_m[m]*prob.fexpl(uf[m,P-1])
  u[m]   = prob.solve(prob.coll.delta_m[m], b)

#print 'diff predictor:',abs(u[M-1] - np.exp((prob.lambda_fast + prob.lambda_slow)*tend))

#oldG= np.zeros((M,1)) #not needed as vector..... really?
for k in range(0,K):
  #first finestep
  uf[0,0]=u0+prob.coll_fast[0].delta_m[0]*prob.fexpl(u0)
  #for p in range(0, M):
  oldG[p] = prob.fimpl(u[0])
  #integral term int_s^{s+1} using the fast parts only
  prob.update_I_p(uf)
  
  oldF=prob.fexpl(uf[0,0])
  for p in range(1, P):
    #integral term, uses u^{F,k} for u^k, see above with prob.update_I_p
    swap = prob.fexpl(uf[0,p])
    uf[0,p]=uf[0,p-1]+prob.coll_fast[0].delta_m[p]*(prob.fexpl(uf[0,p-1])-oldF+prob.fimpl(u[0])-oldG[0]) #impl terms should cancel u^k_0=u(0)
    uf[0,p] += prob.I_p_pp1[0,p] #integral part
    oldF = swap

  prob.update_I_m(u) #integral part on coarse level
  u[0]=prob.solve(prob.coll.delta_m[0], u0+prob.coll.delta_m[0]*(prob.fexpl(uf[0,P-1])-oldF-prob.fimpl(u[0]))+prob.I_m_mp1[0])
  #u[0]=prob.solve(prob.coll.delta_m[0], u0+prob.coll.delta_m[0]*(prob.fexpl(uf[0,P-1])-oldF[P-1]-prob.fimpl(uf[0,P-1]))+prob.I_m_mp1[0])
  for m in range(1, M):
    for p in range(0, P):
#      oldF[p] = prob.fexpl(uf[m,p])
      oldG[p] = prob.fimpl(u[m,p])
    #update the inner fast values
    oldF = prob.fexpl(uf[m,0])
    for p in range(1, P):
      #integral term, uses u^{F,k} for u^k, see above with prob.update_I_p
      swap = prob.fexpl(uf[m,p])
      uf[m,p] = uf[m,p-1]+prob.coll_fast[m].delta_m[p]*(prob.fexpl(uf[m,p-1])-oldF+prob.fimpl(u[m])-oldG) #impl terms should cancel u^k_0=u(0)
      uf[m,p] += prob.I_p_pp1[m,p] #integral part
      oldF = swap
   #update the slow values with new fast ones
    u[m]=prob.solve(prob.coll.delta_m[m], u[m-1]+prob.coll.delta_m[m]*(prob.fexpl(uf[m,P-1])-oldF-prob.fimpl(u[m]))+ prob.I_m_mp1[m])
    #u[m]=prob.solve(prob.coll.delta_m[m], uf[m-1,P-1]+prob.coll.delta_m[m]*(prob.fexpl(uf[m,P-1])-oldF[P-1]-prob.fimpl(uf[m,P-1]))+ prob.I_m_mp1[m])
#  print abs(prob.end_value(u,1.0) - np.exp((prob.lambda_fast + prob.lambda_slow)))
#  print 'diff:',k,abs(u[M-1] - np.exp((prob.lambda_fast + prob.lambda_slow)*tend)) #, u[M-1]
#  print 'res res/res_sub:',prob.get_residual(u, 1.0),prob.get_residual_sub(uf, 1.0)
print 'diff:',k,abs(u[M-1] - np.exp((prob.lambda_fast + prob.lambda_slow)*tend)) #, u[M-1]
#print 'diff-endvaluecorrection:',k,abs(prob.end_value(u, 1.0) - np.exp((prob.lambda_fast + prob.lambda_slow)*tend)) #, u[M-1]

#collSol=prob.get_coll_solution_sub(1.0)
#print 'diffCollSol:', abs(collSol[M-1,P-1]-np.exp((prob.lambda_fast+prob.lambda_slow)*tend))
#print 'diffFineSol:', abs(uf[M-1,P-1]-np.exp((prob.lambda_fast+prob.lambda_slow)*tend))
#print 'res res/res_sub:',prob.get_residual(u, 1.0),prob.get_residual_sub(uf, 1.0)
