from firedrake import *
import numpy as np

class MovingHeat:

    def __init__(self, nu, alpha, nx, ny):
        self.nu = nu
        self.alpha = alpha
        self.meshA = RectangleMesh(nx, ny, 0.5, 0.5)
        self.meshB = RectangleMesh(nx, 4*ny, 1.0, 4.0)
        #shift meshA to the right initial position
        Vc    = self.meshA.coordinates.function_space()
        f2    = Function(Vc).interpolate(Expression(("x[0]-0.5","x[1]+2.0")))
        self.meshA.coordinates.assign(f2)
        self.meshA.clear_spatial_index()
        self.V_B = FunctionSpace(self.meshB, "CG", 1)
        Vector_V_B = VectorFunctionSpace(self.meshB, "CG", 1)
        self.phi = TrialFunction(self.V_B)
        self.v   = TestFunction(self.V_B)
        self.myf  = Function(self.V_B, name="TransferFunction")
        self.boundary_nodes_B = DirichletBC(self.V_B, 0, 1).nodes
        c = Function(Vector_V_B)
        c.interpolate(SpatialCoordinate(self.meshB))
        self.Xs = c.dat.data_ro[self.boundary_nodes_B, :]
        self.Xs = np.array(map(lambda p: p-np.array([1.0e-10, 0.0]), self.Xs))

        self.lastTime=0.0

    def evalFast(self, u_B, t): 
        """evaluates the fast part
        """
        #move meshA to the current position
        timestep = t-self.lastTime
        self.lastTime = t
        Vc = self.meshA.coordinates.function_space()
        f2  = Function(Vc).interpolate(Expression(("x[0]", "x[1] - shift"), shift = 0.1*timestep ))
        if hasattr(self.meshA, 'spatial_index'):
            self.meshA.clear_spatial_index()
            #del self.meshA.spatial_index #???
        #self.meshA.coordinates.assign(f2)  
        self.meshA=Mesh(Function(f2))
        #create functionspace on the moving part
        V_A = FunctionSpace(self.meshA, "CG", 1)
        self.u_A = Function(V_A, name="TemperatureA")
        va = TestFunction(V_A)
        #set temperature of the moving part
        self.u_A.interpolate(Expression("5.0*t",t=t))
        # Evaluate u_A at meshB boundary nodes with Id 1
        node_vec = self.u_A.at(self.Xs, dont_raise=True)
        # Remove None's 
        node_vec = [0.0 if x is None else x for x in node_vec] 
        print node_vec
        self.myf.dat.data[self.boundary_nodes_B] = node_vec
        myfOut=File("myf.pvd")
        myfOut.write(self.myf, time=t)
        #define the linear form of the fast part
        #Rfast=(self.alpha*inner((self.myf-u_B), self.v)*ds(1,domain=self.meshB))
        Rfast=(self.alpha*inner((self.myf), self.v)*ds(1,domain=self.meshB))
        ass=assemble(Rfast)
        print np.sum(ass.vector().get_local())
        swap=assemble(inner(self.u_A, va)*dx)
        print np.sum(swap.vector().get_local())
        return ass

    def evalSlow(self, uB, t):
        """evaluate the slow part
        """
        Rslow=-self.nu*inner(grad(uB), grad(self.v))*dx(domain=self.meshB)
        return assemble(Rslow)

    def evalF(self, u, t):
        """ evaluate the sum of both parts
        """
        return evalSlow(u,t)+evalFast(u,t)

    def evalMmJ(self, u, t, r, a):
        """ evaluate (M-a*J)u=r
        """
        bif = (inner(self.phi, self.v) + a* self.nu * inner(grad(self.phi), grad(self.v)))*dx(domain=self.meshB)
        bifM=assemble(bif)
        u_B=Function(u, name="solEvalMmJ")
        solve(bifM, u_B, r)
        return u_B

    def evalM(self, u):
        """ evaluate M*u
        """
        aM=inner(u,self.v)*dx
        return assemble(aM)

    def solveM(self, r):
        """ evaluate u, such that Mu=r
        """
        bif = (inner(self.phi, self.v))*dx(domain=self.meshB)
        bifM=assemble(bif)
        u_B=Function(self.V_B, name="solM")
        solve(bifM, u_B, r)
        return u_B

