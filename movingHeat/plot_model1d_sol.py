import numpy as np
import copy

from matplotlib import pyplot as plt
from pylab import rcParams
from matplotlib.ticker import ScalarFormatter
from subprocess import call

nu    = 0.5
k     = 1.0
omega = 1.0
y_0   = 1.0
v_0   = 2.0

def f(y):
  return -nu*k*k*y

def solve_f(y, dt):
  return y/(1.0 + dt*nu*k*k)

def g(t):
  return np.exp(1j*omega*t)*(y_0 - v_0)/np.sqrt(2*np.pi)

tend   = 100.0
nsteps = 10
taxis  = np.linspace(0,tend,nsteps+1)
dt     = tend/float(nsteps)

c          = y_0 - (y_0 - v_0)/( np.sqrt(2*np.pi)*(k*k*nu + 1j*omega) )
taxis_plot = np.linspace(0,tend,1e5)
y_exact    = c*np.exp(-nu*k*k*taxis_plot) + (y_0 - v_0)*np.exp(1j*taxis_plot*omega)/( np.sqrt(2*np.pi)*(k*k*nu + 1j*omega))

y_euler = np.zeros(nsteps+1, dtype='complex')
y_euler[0] = y_0
for n in range(nsteps):
  t            = float(n)*dt
  y_temp       = y_euler[n] + dt*g(t)
  y_euler[n+1] = solve_f(y_temp, dt)

  # explicit Euler works correctly
  #y_euler[n+1] = y_euler[n] + dt*(f(y_euler[n]) + g(t))

print ("Error: %5.3e" % abs(y_euler[-1] - y_exact[-1]))
fig = plt.figure()
fig.clear()
plt.plot(taxis_plot, y_exact.real, 'b')
plt.plot(taxis, y_euler.real, 'r')
plt.ylim([-1.0, 2.0])
plt.show()