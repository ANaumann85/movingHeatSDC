import matplotlib.pyplot as plt
import numpy as np

data=np.loadtxt('stabEigs_ros2.mtx')
n=data.shape[0]
nu_vec=data[0][0:n-1]
alpha_vec=data[0:n-1,0]
stab_ros=data[1:,1:]

fig  = plt.figure()
levels = np.array([0.25, 0.5, 0.75, 0.9, 1.1])
plt.contour(nu_vec, alpha_vec, np.absolute(stab_ros), [1.0],  colors='k', linestyles='dotted')
plt.xlabel('$\\nu$')
plt.ylabel('$\\alpha$')
fig.savefig("stability_ros2.pdf", bbox_inches='tight')
