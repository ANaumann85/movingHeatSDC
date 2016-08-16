from problem_model import problem_model
from sdc import sdc_step
import numpy as np
import copy

from pylab import rcParams
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from subprocess import call

M = 3
P = 3
tstart = 0.0
tend   = 1.0

a_vec  = np.linspace(0.0,   2.0, 20)
nu_vec = np.linspace(-5.0,  0.0, 50)
stab   = np.zeros((np.size(a_vec), np.size(nu_vec)))

K_iter = 2

for kk in range(np.size(a_vec)):
  for ll in range(np.size(nu_vec)):

    prob = problem_model(a_vec[kk], nu_vec[ll])
    sdc    = sdc_step(M, P, 0.0, 1.0, prob)

    R = np.zeros((2,2))
    u0 = np.eye(2)
    
    for mm in range(2):
    
      # reset buffers to zero
      u     = np.zeros((M,2))
      usub  = np.zeros((M,P,2))
      u_    = np.zeros((M,2))
      usub_ = np.zeros((M,P,2))

      sdc.predict(u0[:,mm], u_, usub_)
      for k in range(K_iter):
        sdc.sweep(u0[:,mm], u, usub, u_, usub_)
        u_    = copy.deepcopy(u)
        usub_ = copy.deepcopy(usub)
      R[:,mm] = u[M-1]
  
    eval, evec = np.linalg.eig(R)
    stab[kk,ll] = np.linalg.norm(eval, np.inf)

#rcParams['figure.figsize'] = 1.5, 1.5
fs = 16
fig  = plt.figure()
levels = np.array([0.25, 0.5, 0.75, 0.9, 1.1])
CS1 = plt.contour(nu_vec, a_vec, np.absolute(stab), levels, colors='k', linestyles='dashed')
CS2 = plt.contour(nu_vec, a_vec, np.absolute(stab), [1.0],  colors='k')
plt.clabel(CS1, inline=True, fmt='%3.2f', fontsize=fs-2)
#plt.gca().add_patch(Polygon([[0, 0], [lambda_1[-1],0], [lambda_1[-1],lambda_1[-1]]], visible=True, fill=True, facecolor='.75',edgecolor='k', linewidth=1.0,  zorder=11))
plt.xlabel(r'$\nu$', fontsize=fs)
plt.ylabel(r'$a$', fontsize=fs)
filename = 'stability-K'+str(K_iter)+'-M'+str(M)+'-P'+str(P)+'.pdf'
fig.savefig(filename, bbox_inches='tight')
call(["pdfcrop", filename, filename])
#plt.show()
