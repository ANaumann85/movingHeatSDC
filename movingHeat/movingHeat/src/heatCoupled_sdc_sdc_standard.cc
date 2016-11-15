#include "HeatCoupled.h"
#include "MRSdc.h"
#include "Sdc.h"
#include <sstream>
#include "Timer.h"

namespace Dune
{
  template< typename A, typename B >
  inline void axpy(double a, const BlockVector<A, B>& x, BlockVector<A, B>& y)
  { y.axpy(a, x); }

  template< typename A, typename B >
  inline void setValue(BlockVector<A, B>& x, const double& v)
  { x=v; }
}

namespace std 
{
  template< typename A, unsigned long s >
  inline void axpy(double a, const array<A, s>& x, array<A, s>& y)
  { for(unsigned i(0); i < s; ++i) axpy(a, x[i], y[i]); }
  
  template< typename A, unsigned long s >
  inline void setValue(array<A, s>& x, const double& v)
  { for(auto& d : x) setValue(d, v); }
}

template<int M, int P>
void solveMRSDC(HeatCoupled& heat, unsigned k_iter, unsigned nStep)
{
  typedef MRSdc<HeatCoupled::VectorType, M,P> Method;
  typename Method::Init init([&heat](HeatCoupled::VectorType& d) { heat.init(d); });
  Method sdc(init, k_iter, "radau_right", "radau_right");
  HeatCoupled::VectorType y0;
  //heat.init(y0, [](auto x) { return (x[0]-0.5)*(x[0]-0.5)*(x[1]-2)*(x[1]-2); });
  heat.init(y0); //y0 = 0.0;
  y0[0]=1.0;
  y0[1]=5.0;

  double t0(0.0), tend(20.0);

  /*std::stringstream ss; ss << "heat_sdc_M-" << M << "_P-" << P << "_nStep-" << nStep;
  heat.startFile(ss.str(), y0);
  heat.writeResult(t0);*/
  { MyTimer timer("mrsdc-solve");
  sdc.solve(heat, y0, t0, tend, nStep);
  }
  //heat.writeResult(tend);
}

template<int M>
void solveIMEX(HeatCoupled& heat, unsigned k_iter, unsigned nStep)
{
  typedef Sdc<HeatCoupled::VectorType, M> Method;
  typename Method::Init init([&heat](HeatCoupled::VectorType& d) { heat.init(d); });
  Method sdc(init, k_iter, "radau_right");
  HeatCoupled::VectorType y0;
  //heat.init(y0, [](auto x) { return (x[0]-0.5)*(x[0]-0.5)*(x[1]-2)*(x[1]-2); });
  heat.init(y0); //y0 = 0.0;
  y0[0]=1.0;
  y0[1]=5.0;

  double t0(0.0), tend(20.0);

  /*std::stringstream ss; ss << "heat_sdc_standard_M-" << M << "_nStep-" << nStep;
  heat.startFile(ss.str(), y0);
  heat.writeResult(t0);*/
  { MyTimer timer("sdc-solve");
  sdc.solve(heat, y0, t0, tend, nStep);
  }
  //heat.writeResult(tend);
}


int main(int argc, char* argv[])
{
  if(argc < 2) {
    std::cerr << "usage: " << argv[0] << " <nStep>\n";
    return 1;
  }
  static const unsigned M=3;
  static const unsigned P=8;
  static const unsigned kIter=2;
  //mpi-helper from dune
  MPIHelper::instance(argc, argv);
  HeatCoupled heat(20, 1.0e-3, 1.0e-3, 5.0, 0.0, false, false);

  unsigned nStep(40);
  {std::stringstream ss; ss << argv[1]; ss >> nStep; }

  solveMRSDC<M,P>(heat, kIter, nStep);
  solveIMEX<M>(heat, kIter, nStep);
  return 0;
}

