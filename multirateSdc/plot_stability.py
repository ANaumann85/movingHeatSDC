from problem import problem
from sdc import sdc_step
import numpy as np
import copy

from pylab import rcParams
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from subprocess import call

M = 2
P = 2
tstart = 0.0
tend   = 1.0

lambda_1 = np.linspace(-1.0,  0.0, 25)
lambda_2 = np.linspace(-5.0,  0.0, 20)
stab     = np.zeros((np.size(lambda_1), np.size(lambda_2)))

K_iter = 2

u0    = 2.0

for kk in range(np.size(lambda_1)):
  for ll in range(np.size(lambda_2)):

    prob = problem(lambda_1[kk], lambda_2[ll])
    sdc    = sdc_step(M, P, 0.0, 1.0, prob)
    
    # reset buffers to zero
    u     = np.zeros((M,1))
    usub  = np.zeros((M,P))
    u_    = np.zeros((M,1))
    usub_ = np.zeros((M,P))

    sdc.predict(u0, u_, usub_)
    for k in range(K_iter):
      sdc.sweep(u0, u, usub, u_, usub_)
      u_    = copy.deepcopy(u)
      usub_ = copy.deepcopy(usub)

    stab[kk,ll] = u[M-1]

#rcParams['figure.figsize'] = 1.5, 1.5
fs = 16
fig  = plt.figure()
levels = np.array([0.25, 0.5, 0.75, 0.9, 1.1])
CS1 = plt.contour(lambda_2, lambda_1, np.absolute(stab), levels, colors='k', linestyles='dashed')
CS2 = plt.contour(lambda_2, lambda_1, np.absolute(stab), [1.0],  colors='k')
plt.clabel(CS1, inline=True, fmt='%3.2f', fontsize=fs-2)
#plt.gca().add_patch(Polygon([[0, 0], [lambda_1[-1],0], [lambda_1[-1],lambda_1[-1]]], visible=True, fill=True, facecolor='.75',edgecolor='k', linewidth=1.0,  zorder=11))
plt.xlabel(r'$\lambda_{\rm explicit}$', fontsize=fs)
plt.ylabel(r'$\lambda_{\rm implicit}$', fontsize=fs)
plt.show()
