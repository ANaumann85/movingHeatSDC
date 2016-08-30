#include "Heat.h"
#include <sstream>

int main(int argc, char* argv[])
{
  //mpi-helper from dune
  MPIHelper::instance(argc, argv);
  Heat heat(10);
  Heat::VectorType y0, rhs, fVal;
  //heat.init(y0, [](auto x) { return (x[0]-0.5)*(x[0]-0.5)*(x[1]-2)*(x[1]-2); });
  heat.init(y0); y0 = 0.0;
  heat.init(rhs);
  heat.init(fVal);

  double t0(0.0), tend(20.0);
  unsigned nStep(40);

  double dt((tend-t0)/nStep);
  heat.startFile("heat_imexEuler", y0);
  heat.writeResult(t0);
  for(unsigned n(0); n < nStep; ++n) {
    heat.updateMatrix(t0, dt);
    heat.Mv(y0, rhs);
    fVal = 0.0;
    heat.fast(t0+dt, y0, fVal);
    fVal *= dt; rhs += fVal;
    //fVal.axpy(dt, rhs);
    heat.solveMaJ(rhs, y0);
    t0 +=dt;
    /*std::stringstream cfn; cfn << "heat_implEuler_" << t0; 
    heat.writeResult(cfn.str(), y0);*/
    heat.writeResult(t0);
  }
  return 0;
}
