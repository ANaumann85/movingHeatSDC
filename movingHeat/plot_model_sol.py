import numpy as np
import copy

from matplotlib import pyplot as plt
from pylab import rcParams
from matplotlib.ticker import ScalarFormatter
from subprocess import call

xaxis_u = np.linspace(-10, 0, 100)
xaxis_v = np.linspace(0, 10, 100)
dt     = 0.1
nsteps = 80
a      = 0.2
kappa  = [0, 1, 2]
ak0    = np.array([1.0,    0.1, 0.00])
ck0    = np.array([ak0[0], 2.0, 0.0])
b10    = np.sum(ck0 - ak0)*a
bk0    = np.array([0.0, b10, 0.0])
dk0    = -bk0

fig = plt.figure()

for n in range(nsteps):
  t = dt*float(n)
  u = 0.0*xaxis_u
  v = 0.0*xaxis_v
  for k in kappa:
    ak = ak0[k]*np.exp(-float(k)**2*t)
    bk = bk0[k]*np.exp(-float(k)**2*t)
    ck = ck0[k]*np.exp(-float(k)**2*t)
    dk = dk0[k]*np.exp(-float(k)**2*t)
    
    u  += ak*np.cos(k*xaxis_u) + bk*np.sin(k*xaxis_u)
    v  += ck*np.cos(k*xaxis_v) + dk*np.sin(k*xaxis_v)
  fig.clear()
  plt.plot(xaxis_u, u, 'r')
  plt.plot(xaxis_v, v, 'g')
  plt.ylim([-1.0, 2.0])
  plt.title(("t = %3.2f" % t))
  plt.show(block=False)
  plt.pause(0.04)
