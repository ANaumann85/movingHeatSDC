#include "Heat.h"
#include "MRSdc.h"
#include <sstream>

namespace Dune
{
  template< typename A, typename B >
  inline void axpy(double a, const BlockVector<A, B>& x, BlockVector<A, B>& y)
  { y.axpy(a, x); }

  template< typename A, typename B >
  inline void setValue(BlockVector<A, B>& x, const double& v)
  { x=v; }
}

template<int M, int P>
void solve(Heat& heat, unsigned k_iter, unsigned nStep)
{
  typedef MRSdc<Heat::VectorType, M,P> Method;
  Method::Init init([&heat](Heat::VectorType& d) { heat.init(d); });
  Method sdc(init, k_iter, "radau_right", "equi_noleft");
  Heat::VectorType y0;
  //heat.init(y0, [](auto x) { return (x[0]-0.5)*(x[0]-0.5)*(x[1]-2)*(x[1]-2); });
  heat.init(y0); y0 = 0.0;

  double t0(0.0), tend(20.0);

  std::stringstream ss; ss << "heat_sdc_M-" << M << "_P-" << P << "_nStep-" << nStep;
  heat.startFile(ss.str(), y0);
  heat.writeResult(t0);
  sdc.solve(heat, y0, t0, tend, nStep);
  heat.writeResult(tend);
}

int main(int argc, char* argv[])
{
  if(argc < 2) {
    std::cerr << "usage: " << argv[0] << " <nStep>\n";
    return 1;
  }
  //mpi-helper from dune
  MPIHelper::instance(argc, argv);
  Heat heat(10);

  unsigned nStep(40);
  {std::stringstream ss; ss << argv[1]; ss >> nStep; }
  solve<3,2>(heat, 5, nStep);
  return 0;
}

