from firedrake import *

nx     = 20
ny     = 40

t      = 0.0
tend   = 2.0
nsteps = 5
timestep = (tend-t)/float(nsteps)

meshA = RectangleMesh(10, 10, 0.5, 0.5)
Vc    = meshA.coordinates.function_space()
f2    = Function(Vc).interpolate(Expression(("x[0]-0.5","x[1]+2.0")))
#move the first time
#meshA.coordinates.assign(f2)
meshA = Mesh(Function(f2))

out_A = File('meshA.pvd')
while t < tend:
  V_A = FunctionSpace(meshA, "CG", 1)
  u_A = Function(V_A, name="TemperatureA")
  u_A.interpolate(Expression("1.0*t", t=t))

  # Check value at one specific point, it is allways inside 
  print 'uA(-1.0e-15,2.2):', u_A.at((-1.0e-15, 2.2), dont_raise=True)
  print 'uA(0,2.2):', u_A.at((0, 2.2), dont_raise=True)
  out_A.write(u_A)
  Vc = meshA.coordinates.function_space()
  f2  = Function(Vc).interpolate(Expression(("x[0]", "x[1] - shift"), shift = 0.1*timestep ))
#  del  meshA.spatial_index
#  meshA.coordinates.assign(f2)  
#  del  meshA.spatial_index
  meshA = Mesh(Function(f2))

  t += timestep
