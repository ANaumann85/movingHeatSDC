from problem import *
import numpy as np
import copy
from MultirateCollocation import multirateCollocation

lambda_1 = -1.0
lambda_2 = -0.0
u0       = 1.0

M = 4  # SDC nodes for slow process
P = 4  # SDC nodes for fast process
K = 4 # Number of sweeps


tend     = 0.5
nsteps   = 1
dt       = tend/float(nsteps)

I_m_mp1 = np.zeros(M)
I_p_pp1 = np.zeros((M,P))

def solve_f1(a, b):
  return b/(1.0 - a*lambda_1)

def f1(u):
  return lambda_1*u

def f2(u):
  return lambda_2*u

def uex(t):
  return u0*np.exp((lambda_1+lambda_2)*t)

def update_I_m_mp1(u, usub, coll):
  for m in range(M):
    I_m_mp1[m] = coll.integrate_m_mp1(f1(u), m) + coll.integrate_m_mp1_sub(f2(usub[m,:]), m)

def update_I_p_pp1(u, usub, coll):
  for m in range(M):
    for p in range(P):
      I_p_pp1[m,p] = coll.integrate_p_pp1(f1(u), m, p) + coll.integrate_p_pp1_sub(f2(usub[m,:]), m, p)


u  = np.zeros(M)
u_ = np.zeros(M)
# initialise for first time step
u0_step = u0

usub  = np.zeros((M,P))
usub_ = np.zeros((M,P))
u0sub_step = u0

for n in range(nsteps):
  coll = multirateCollocation(M, P, float(n)*dt, float(n+1)*dt)
  print "=============================="
  print ("Time step from %4.3f to %4.3f" % (coll.coll.tleft, coll.coll.tright))
  
  #### Run predictor ###
  for m in range(M):
   
    # single implicit step in f_1
    if m==0:
      u[0] = solve_f1(coll.coll.delta_m[0], u0_step)
    else:
      u[m] = solve_f1(coll.coll.delta_m[m], u[m-1])

    # multiple explicit steps in f_2 holding f_1 term constant
    for p in range(P):
      if p==0:
        # Sum up right hand side terms in SDC sweep
        usub[m,0] = u0sub_step + coll.coll_sub[m].delta_m[0]*( f1(u[m]) + f2(u0sub_step) )
      else:            
        usub[m,p] = usub[m,p-1] + coll.coll_sub[m].delta_m[p]*( f1(u[m]) + f2(usub[m,p-1]) )
    
    # update substep initial value
    u0sub_step = usub[m,P-1]

  ### Run SDC iteration ###
  for k in range(K):

    update_I_m_mp1(u, usub, coll)
    update_I_p_pp1(u, usub, coll)

    # Prepare for next iteration
    u_    = copy.deepcopy(u)
    usub_ = copy.deepcopy(usub)
    
    for m in range(M):
      if m==0:
        rhs  = u0_step - coll.coll.delta_m[0]*f1(u_[0]) + I_m_mp1[0]
        u[0] = solve_f1(coll.coll.delta_m[0], rhs)
      else:
        rhs  = u[m-1] - coll.coll.delta_m[m]*f1(u_[m]) + I_m_mp1[m]
        u[m] = solve_f1(coll.coll.delta_m[m], rhs)

      # Run small time steps; f1 term remains constant
      slow = f1(u[m]) - f1(u_[m])
      for p in range(P):
        if p==0:
          # f2 term cancels out here
          usub[m,0] = u0sub_step + coll.coll_sub[m].delta_m[0]*slow + I_p_pp1[m,0]
        else:
          usub[m,p] = usub[m,p-1] + coll.coll_sub[m].delta_m[p]*( slow + f2(usub[m,p-1]) - f2(usub_[m,p-1]) ) + I_p_pp1[m,p]

      # update substep initial value
      u0sub_step = usub[m,P-1]

  # update coarse step initial value for next time step
  u0_step = u[M-1]

print ("Error on coarse level: %4.3e" % abs(u0*np.exp(lambda_1*tend) - u[M-1]))
print ("Error on fine level:   %4.3e" % abs(uex(tend) - usub[M-1,P-1]))
