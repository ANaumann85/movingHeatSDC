#include "Heat.h"
#include "Ros2.h"
#include <sstream>
#include "Timer.h"

namespace Dune
{
  template< typename A, typename B >
  inline void axpy(double a, const BlockVector<A, B>& x, BlockVector<A, B>& y)
  { y.axpy(a, x); }
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
  Ros2<Heat::VectorType > ros2([&heat](auto& in) { heat.init(in); });
  Heat::VectorType y0;
  //heat.init(y0, [](auto x) { return (x[0]-0.5)*(x[0]-0.5)*(x[1]-2)*(x[1]-2); });
  heat.init(y0); y0 = 0.0;

  double t0(0.0), tend(20.0);
  unsigned nStep(40);
  {std::stringstream ss; ss << argv[1]; ss >> nStep; }

  std::stringstream ss; ss << "heat_ros2_" << nStep;
  heat.startFile(ss.str(), y0);
  heat.writeResult(t0);
  { MyTimer timer("ros2-solve");
  ros2.solve(heat, y0, t0, tend, nStep);
  }
  heat.writeResult(tend);
  return 0;
}

