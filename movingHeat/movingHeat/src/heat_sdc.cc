#include "Heat.h"
#include "MRSdc.h"
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

template<int M, int P>
void solve(Heat& heat, Heat::VectorType& y0, unsigned k_iter, unsigned nStep, double tend)
{
  typedef MRSdc<Heat::VectorType, M,P> Method;
  typename Method::Init init([&heat](Heat::VectorType& d) { heat.init(d); });
  Method sdc(init, k_iter, "radau_right", "radau_right");
  //heat.init(y0, [](auto x) { return (x[0]-0.5)*(x[0]-0.5)*(x[1]-2)*(x[1]-2); });
  //heat.init(y0); y0 = heat.getBVal();

  double t0(0.0);

  std::stringstream ss; ss << "heat_sdc_M-" << M << "_P-" << P << "_nStep-" << nStep;
  heat.startFile(ss.str(), y0);
  heat.writeResult(t0);
  { MyTimer timer("mrsdc-solve");
  sdc.solve(heat, y0, t0, tend, nStep);
  }
  heat.writeResult(tend);
}

template<int M>
void solve(Heat& heat, Heat::VectorType& y0, unsigned k_iter, unsigned nStep, unsigned P, double tend)
{
  switch(P) 
  {
    case 1:
      solve<M,1>(heat, y0, k_iter, nStep, tend);
      break;
    case 2:
      solve<M,2>(heat, y0, k_iter, nStep, tend);
      break;
    case 3:
      solve<M,3>(heat, y0, k_iter, nStep, tend);
      break;
    case 4:
      solve<M,4>(heat, y0, k_iter, nStep, tend);
      break;
    case 5:
      solve<M,5>(heat, y0, k_iter, nStep, tend);
      break;
    case 6:
      solve<M,6>(heat, y0, k_iter, nStep, tend);
      break;
    case 7:
      solve<M,7>(heat, y0, k_iter, nStep, tend);
      break;
    case 8:
      solve<M,8>(heat, y0, k_iter, nStep, tend);
      break;
    default:
      throw std::runtime_error("not supported");
  } 
}

void solve(Heat& heat, Heat::VectorType& y0, unsigned k_iter, unsigned nStep, unsigned M, unsigned P, double tend)
{
  switch(M) {
    case 1:
      solve<1>(heat, y0, k_iter, nStep, P, tend);
      break;
    case 2:
      solve<2>(heat, y0, k_iter, nStep, P, tend);
      break;
    case 3:
      solve<3>(heat, y0, k_iter, nStep, P, tend);
      break;
    case 4:
      solve<4>(heat, y0, k_iter, nStep, P, tend);
      break;
    case 5:
      solve<5>(heat, y0, k_iter, nStep, P, tend);
      break;
    case 6:
      solve<6>(heat, y0, k_iter, nStep, P, tend);
      break;
    case 7:
      solve<7>(heat, y0, k_iter, nStep, P, tend);
      break;
    case 8:
      solve<8>(heat, y0, k_iter, nStep, P, tend);
      break;
  }
}

int main(int argc, char* argv[])
{
  if(argc < 6) {
    std::cerr << "usage: " << argv[0] << " <nStep> <M> <P> <kIter> <laplExplHExpo>\n";
    return 1;
  }
  //mpi-helper from dune
  MPIHelper::instance(argc, argv);
  double tend(20.0);
  double nu(1.0e-3), alpha(1.0e-3), v0(5.0), src(0.0);
  bool useLapl0(false), addConstRobin(false);

  unsigned nStep(40), kIter(2);
  {std::stringstream ss; ss << argv[1]; ss >> nStep; }
  {std::stringstream ss; ss << argv[4]; ss >> kIter; }
  unsigned M, P;
  {std::stringstream ss; ss << argv[2]; ss >> M; }
  {std::stringstream ss; ss << argv[3]; ss >> P; }

#if 0
  Heat::LaplTilde laplTilde;
  {std::stringstream ss ; ss << argv[5] ; ss >> laplTilde.mode; }
  laplTilde.fac = 1.0;
  switch(laplTilde.mode) {
    case 1:
    case 2:
      break;
    case 3:
    case 4:
      laplTilde.mode -= 2;
      laplTilde.fac = tend/nStep;
  }
#endif
  unsigned laplHExpo(0);
  {std::stringstream ss ; ss << argv[5] ; ss >> laplHExpo; }

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
  }
  Heat heat(nInter, nu, alpha, v0, src, useLapl0, addConstRobin, laplExplFac);//, laplTilde
  Heat::VectorType y0;
  heat.init(y0); 
  if(addConstRobin) {
    y0 = heat.getBVal();
    heat.setbAlph(1e-2);
  }else {
    y0 = 0.0;
  }

  solve(heat, y0, kIter, nStep, M,P, tend);
  return 0;
}

