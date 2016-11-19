#include "Heat.h"
#include "RKC.h"
#include <sstream>
#include "Timer.h"

namespace Dune
{
  template< typename A, typename B >
  inline void axpy(double a, const BlockVector<A, B>& x, BlockVector<A, B>& y)
  { y.axpy(a, x); }
}

struct NP
{ void operator()(const double& , const Heat::VectorType& ) const {} };
int main(int argc, char* argv[])
{
  if(argc < 3) {
    std::cerr << "usage: " << argv[0] << "<nStage> <nStep>\n";
    return 1;
  }
  //mpi-helper from dune
  MPIHelper::instance(argc, argv);
  double nu(1.0e-3), alpha(1.0e-3), v0(5.0), src(0.0);
  bool useLapl0(false), addConstRobin(false);
  Heat heat(10, nu, alpha, v0, src, useLapl0, addConstRobin);
  Heat::VectorType y0;
  heat.init(y0); 
  if(addConstRobin) {
    y0 = heat.getBVal();
    heat.setbAlph(1e-2);
  }else {
    y0 = 0.0;
  }

  unsigned nStage(2), order(2);
  {std::stringstream ss; ss << argv[1]; ss >> nStage; }
  RKC<Heat::VectorType > rkc([&heat](auto& in) { heat.init(in); }, nStage, order);
  //heat.init(y0, [](auto x) { return (x[0]-0.5)*(x[0]-0.5)*(x[1]-2)*(x[1]-2); });

  double t0(0.0), tend(20.0);
  unsigned nStep(40);
  {std::stringstream ss; ss << argv[2]; ss >> nStep; }

  std::stringstream ss; ss << "heat_rkc_" << nStep;
  heat.startFile(ss.str(), y0);
  heat.writeResult(t0);
  NP np;
  { MyTimer timer("rkc-solve");
  rkc.solve(heat, y0, t0, tend, nStep, np);
  }
  heat.writeResult(tend);
  return 0;
}
