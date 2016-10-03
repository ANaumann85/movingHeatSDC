from problem_model import problem_model
from sdc import sdc_step
from sdc_standard import sdc_standard_step
from ros2 import ros2_step
import numpy as np
import copy

#from pylab import rcParams

M = 3
tstart = 0.0
tend   = 1.0
nue=50
a = 1.0
alpha_vec = np.linspace(0.0, 35.0, 20)
nu_vec    = np.linspace(0.0, nue, 20)
P_vec = range(M,M+11)

stab     = np.zeros((np.size(alpha_vec), np.size(nu_vec)))

K_iter = 2

for P in P_vec:
  for kk in range(np.size(alpha_vec)):
    for ll in range(np.size(nu_vec)):
      prob = problem_model(a=a, nu=nu_vec[ll], alpha=alpha_vec[kk], v0=0.0)
      
      sdc     = sdc_step(M, P, 0.0, 1.0, prob)    
    
      R     = np.zeros((prob.dim,prob.dim))
      u0    = np.eye(prob.dim)
    
      for mm in range(prob.dim):
        #### run multirate SDC
        u     = np.zeros((M,prob.dim))
        usub  = np.zeros((M,P,prob.dim))
        u_    = np.zeros((M,prob.dim))
        usub_ = np.zeros((M,P,prob.dim))
        fu     = np.zeros((M,prob.dim))
        fu_sub  = np.zeros((M,P,prob.dim))

        sdc.predict(u0[:,mm], u, usub, fu, fu_sub)
        for k in range(K_iter):
          sdc.sweep(u0[:,mm], u, usub, fu, fu_sub)
        #res = sdc.residual(u0[:,mm], u)
        #res_sub = sdc.sub_residual(u0[:,mm], usub)
        #print ("standard residual: %5.3e" % res)
        #print ("embedded residual: %5.3e" % res_sub)
        R[:,mm] = u[M-1]
  
      eval, evec  = np.linalg.eig(R)
      stab[kk,ll] = np.linalg.norm(eval, np.inf)
  fname='stability-K-%d-M-%d-P-%d-nue_%d.txt' %(K_iter, M, P, nue)
  print "save:",fname
  np.savetxt(fname, stab)
    
#rcParams['figure.figsize'] = 1.5, 1.5
#fs = 16
#fig  = plt.figure()
#levels = np.array([0.25, 0.5, 0.75, 0.9, 1.1])
#
##CS1 = plt.contour(nu_vec, alpha_vec, np.absolute(stab), levels, colors='k', linestyles='dashed')
#CS2 = plt.contour(nu_vec, alpha_vec, np.absolute(stab), [1.0],  colors='k')
#CS3 = plt.contour(nu_vec, alpha_vec, np.absolute(stab_std), [1.0],  colors='k', linestyles='dashed')
#CS3 = plt.contour(nu_vec, alpha_vec, np.absolute(stab_ros), [1.0],  colors='k', linestyles='dotted')
#
##CS3 = plt.contour(nu_vec, alpha_vec, stab_c, [0.0],  colors='r', linestyle='dotted')
#
##plt.clabel(CS1, inline=True, fmt='%3.2f', fontsize=fs-2)
##plt.gca().add_patch(Polygon([[0, 0], [lambda_1[-1],0], [lambda_1[-1],lambda_1[-1]]], visible=True, fill=True, facecolor='.75',edgecolor='k', linewidth=1.0,  zorder=11))
#plt.xlabel(r'$\nu$', fontsize=fs)
##plt.ylabel(r'$a$', fontsize=fs)
#plt.ylabel(r'$\alpha$', fontsize=fs)
#filename = 'stability-K'+str(K_iter)+'-M'+str(M)+'-P'+str(P)+'-nue_'+str(nue)+'.pdf'
#fig.savefig(filename, bbox_inches='tight')
#call(["pdfcrop", filename, filename])
#plt.show()


