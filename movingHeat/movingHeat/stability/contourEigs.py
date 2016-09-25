import matplotlib.pyplot as plt
import numpy as np

data=np.loadtxt('stabEigs_ros2.mtx')
n=data.shape[0]
nu_vec=data[0,1:n]
alpha_vec=data[1:n,0]
stab_ros=data[1:,1:]

#print stab_ros
fig  = plt.figure()
#levels = np.array([0.25, 0.5, 0.75, 0.9, 1.1])
levels = np.array([0.8, 1.0, 1.2])
#CS=plt.contour(nu_vec, alpha_vec, (stab_ros), [1.0],  colors='k', linestyles='dotted')
CS=plt.contour(nu_vec, alpha_vec, (stab_ros), levels)
plt.clabel(CS,inline=1)
plt.xlabel(r'$\nu$')
plt.ylabel(r'$\alpha$')
fig.savefig("stability_ros2.pdf", bbox_inches='tight')
