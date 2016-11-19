#include "Heat.h"
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

template<int M>
void solve(Heat& heat, Heat::VectorType& y0, unsigned k_iter, unsigned nStep)
{
  typedef Sdc<Heat::VectorType, M> Method;
  typename Method::Init init([&heat](Heat::VectorType& d) { heat.init(d); });
  Method sdc(init, k_iter, "radau_right");
  //heat.init(y0, [](auto x) { return (x[0]-0.5)*(x[0]-0.5)*(x[1]-2)*(x[1]-2); });

  double t0(0.0), tend(20.0);

  std::stringstream ss; ss << "heat_sdc_standard_M-" << M << "_nStep-" << nStep;
  heat.startFile(ss.str(), y0);
  heat.writeResult(t0);
  { MyTimer timer("sdc-solve");
  sdc.solve(heat, y0, t0, tend, nStep);
  }
  heat.writeResult(tend);
}

void solve(Heat& heat, Heat::VectorType& y0, unsigned k_iter, unsigned nStep, unsigned M)
{
  switch(M) {
   case 1:
     solve<1>(heat, y0, k_iter, nStep);
     break;
   case 2:
     solve<2>(heat, y0, k_iter, nStep);
     break;
   case 3:
     solve<3>(heat, y0, k_iter, nStep);
     break;
   case 4:
     solve<4>(heat, y0, k_iter, nStep);
     break;
   case 5:
     solve<5>(heat, y0, k_iter, nStep);
     break;
   case 6:
     solve<6>(heat, y0, k_iter, nStep);
     break;
   case 7:
     solve<7>(heat, y0, k_iter, nStep);
     break;
   case 8:
     solve<8>(heat, y0, k_iter, nStep);
     break;
  }
}

int main(int argc, char* argv[])
{
  if(argc < 5) {
    std::cerr << "usage: " << argv[0] << " <nStep> <M> <P> <kIter> <laplExplHExpo>\n";
    return 1;
  }
  //mpi-helper from dune
  MPIHelper::instance(argc, argv);
  double nu(1.0e-3), alpha(1.0e-3), v0(5.0), src(0.0);
  bool useLapl0(false), addConstRobin(false);
  unsigned laplHExpo(0);
  {std::stringstream ss ; ss << argv[4] ; ss >> laplHExpo; }

  unsigned nInter(10);
  double laplExplFac(0.0) ; //1.0/nInter);
  switch(laplHExpo) 
  {
    case 0:
      laplExplFac = 0.0;
      break;
    case 1: // h
      laplExplFac = 1.0/nInter;
      break;
    case 2: // h^2
      laplExplFac = 1.0/nInter;
      laplExplFac *= laplExplFac;
      break;
    case 3: // sqrt(h) 
      laplExplFac = 1.0/nInter;
      laplExplFac = std::sqrt(laplExplFac);
      break;
    case 100:
      laplExplFac = 1.0;
      break;
  }
  Heat heat(nInter, nu, alpha, v0, src, useLapl0, addConstRobin, laplExplFac);
  Heat::VectorType y0;
  heat.init(y0); 
  if(addConstRobin) {
    y0 = heat.getBVal();
    heat.setbAlph(1e-2);
  }else {
    y0 = 0.0;
  }

  unsigned nStep(40), kIter(2);
  {std::stringstream ss; ss << argv[1]; ss >> nStep; }
  {std::stringstream ss; ss << argv[3]; ss >> kIter; }
  unsigned M;
  {std::stringstream ss; ss << argv[2]; ss >> M; }

  solve(heat, y0, kIter, nStep, M);
  return 0;
}

