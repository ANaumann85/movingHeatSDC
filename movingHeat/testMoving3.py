from heateq_moving3 import *
from firedrake import *

p=MovingHeat(1.0e-3, 1.0e-4, 10, 40)
u0=Function(p.V_B, name="temperatureB")
u0.interpolate(Expression("0.0"))
outS=File("fslow.pvd")
outF=File("ffast.pvd")
for t in range(0,20):
    ffast=Function(p.evalFast(u0, t), name="ffast")
    outF.write(ffast)
    fslow=Function(p.evalSlow(u0, 0.0), name="fslow")
    outS.write(fslow)

#test with implicit euler
tend=20.0
nstep=40
dt=tend/nstep
t=0.0
outT=File("temperature.pvd")
for k in range(0,nstep):
    #r=Function(p.V_B, p.evalM(u0).vector()+dt*p.evalFast(u0, t+dt).vector())
    r=assemble(p.evalM(u0)+dt*p.evalFast(u0, t+dt))
    u0 = p.evalMmJ(u0, t+dt, r, dt)
    outT.write(u0)
    t += dt
