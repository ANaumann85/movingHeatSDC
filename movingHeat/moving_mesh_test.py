from firedrake import *

nu = 1e-3

nx     = 20
ny     = 40

t      = 0.0
tend   = 20.0
nsteps = 40
timestep = (tend-t)/float(nsteps)

meshA = RectangleMesh(10, 10, 0.5, 0.5)
meshB = RectangleMesh(nx, ny, 1.0, 4.0)

# Move meshA 0.5 to left and 2.0 upward to be next to meshB; from http://www.firedrakeproject.org/mesh-coordinates.html
Vc    = meshA.coordinates.function_space()
f2    = Function(Vc).interpolate(Expression(("x[0]-0.5","x[1]+2.0")))
meshA.coordinates.assign(f2)

# Set up problem to be solved on meshB
V_B = FunctionSpace(meshB, "CG", 1)
Vector_V_B = VectorFunctionSpace(meshB, "CG", 1)

phi = TrialFunction(V_B)
v   = TestFunction(V_B)

u_B   = Function(V_B, name="TemperatureB")
u_B_  = Function(V_B, name="TemperatureBOld") 

# The function on meshB to be filled with data from meshA
myf  = Function(V_B, name="TransferFunction")

a = (inner(phi, v) + timestep * nu * inner(grad(phi), grad(v)))*dx(domain=meshB)

boundary_nodes_B = DirichletBC(V_B, 0, 1).nodes
c = Function(Vector_V_B)
c.interpolate(SpatialCoordinate(meshB))
Xs = c.dat.data_ro[boundary_nodes_B, :]

out_A = File('meshA.pvd')
out_B = File('meshB.pvd')

while (t<tend):

  # The function on meshA -- later, this will also be the solution of some PDE
  # but for now just set it to some constant value
  V_A = FunctionSpace(meshA, "CG", 1)
  u_A = Function(V_A, name="TemperatureA")
  u_A.interpolate(Expression("1.0*t", t=t))

  # Evaluate u_A at meshB boundary nodes with Id 1
  node_vec = u_A.at(Xs, dont_raise=True)

  # Check value at one specific point
  print 'uA:', u_A.at((0.0, 2.2+0.1*t), dont_raise=True)
#  del  meshA.spatial_index

  # Remove None's 
  node_vec = [0.0 if x is None else x for x in node_vec] 

  myf.dat.data[boundary_nodes_B] = node_vec

  # Solve on meshB and use myf with data from meshA as Neumann boundary condition
  R = (inner(u_B_, v))*dx + timestep*(inner(myf, v)*ds(1,domain=meshB))
  solve(a == R, u_B)
  
  out_A.write(u_A)
  out_B.write(u_B, myf)

  u_B_.assign(u_B)

  # move meshA downward by timestep; from http://www.firedrakeproject.org/mesh-coordinates.html
  Vc = meshA.coordinates.function_space()
  f2  = Function(Vc).interpolate(Expression(("x[0]", "x[1] - shift"), shift = 0.1*timestep ))
#  del  meshA.spatial_index
#  meshA.coordinates.assign(f2)  
#  del  meshA.spatial_index
  meshA = Mesh(Function(f2))

  t += timestep
