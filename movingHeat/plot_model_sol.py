import numpy as np
import copy

from matplotlib import pyplot as plt
from pylab import rcParams
from matplotlib.ticker import ScalarFormatter
from subprocess import call

xaxis_u = np.linspace(-10, 0, 100)
xaxis_v = np.linspace(0, 10, 100)
dt     = 0.2
nsteps = 40
a      = 5.0
kappa  = [0, 1, 2]
ak0    = np.array([1.0, 0.1, 0.05])
bk0    = np.array([0.0, 0.0, 0.0])
ck0    = np.array([0.25, 0.2, 0.1])
dk0    = -bk0

fig = plt.figure()

for n in range(nsteps):
  t = dt*float(n)
  u = 0.0*xaxis_u
  v = 0.0*xaxis_v
  for k in kappa:
    ak = ak0[k]*np.exp(-float(k)**2*t)
    ck = ck0[k]*np.exp(-float(k)**2*t)
    if k==0.0:
      bk = 0.0
    else:
      bk = (a/float(k))*(ak0[k] - ck0[k])*(np.exp(-float(k)**2*t)-1.0) + bk0[k]
    dk = -bk
    u  += ak*np.cos(k*xaxis_u) + bk*np.sin(k*xaxis_u)
    v  += ck*np.cos(k*xaxis_v) + dk*np.sin(k*xaxis_v)
  fig.clear()
  plt.plot(xaxis_u, u, 'r')
  plt.plot(xaxis_v, v, 'g')
  plt.ylim([-1.0, 2.0])
  plt.title(("t = %3.2f" % t))
  plt.show(block=False)
  plt.pause(0.02)
