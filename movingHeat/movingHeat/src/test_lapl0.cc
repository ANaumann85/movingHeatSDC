#include "Heat.h"
#include <sstream>
#include "Timer.h"
#include <iostream>

namespace Dune
{
  template< typename A, typename B >
  inline void axpy(double a, const BlockVector<A, B>& x, BlockVector<A, B>& y)
  { y.axpy(a, x); }
}

using namespace std;

int main(int argc, char* argv[])
{

  MPIHelper::instance(argc, argv);
  Heat heatOrig(10, 1.0e-3, 0.0e-3, 5.0, 0.0, false );
  Heat heatLapl0(10, 1.0e-3, 0.0e-3, 5.0, 0.0, true);
  Heat::VectorType y0, fo, fl;
  heatOrig.init(y0);
  heatOrig.init(fo);
  heatOrig.init(fl);
  for(unsigned i(0); i < y0.size(); ++i)
    y0[i] = i;
    
  //y0=1.0;
  for(unsigned i(0); i < 20; ++i) {
    double t=i*0.1;
    heatOrig(t, y0, fo);
    heatLapl0(t, y0, fl);
    fo -= fl;
    std::cout << "rdiff adiff norm(fl):" << fo.two_norm()/fl.two_norm() << " " << fo.two_norm() << " " << fl.two_norm() << endl;
  }

  return 0;
}
