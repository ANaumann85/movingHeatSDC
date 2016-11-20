#include "HeatCoupled.h"
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
void solve(HeatCoupled& heat, unsigned k_iter, double thetaFast, unsigned nStep)
{
  typedef MRSdc<HeatCoupled::VectorType, M,P> Method;
  typename Method::Init init([&heat](HeatCoupled::VectorType& d) { heat.init(d); });
  Method sdc(init, k_iter, "radau_right", "radau_right", thetaFast);
  HeatCoupled::VectorType y0;
  //heat.init(y0, [](auto x) { return (x[0]-0.5)*(x[0]-0.5)*(x[1]-2)*(x[1]-2); });
  heat.init(y0);
  y0[0] = heat.getBVal(); //0.0;
  y0[1] = 5.0; 

  double t0(0.0), tend(20.0);

  std::stringstream ss; ss << "heat_sdc_M-" << M << "_P-" << P << "_nStep-" << nStep;
  std::stringstream ss_mv; ss_mv << "heat_sdc_mv_M-" << M << "_P-" << P << "_nStep-" << nStep;
  heat.startFile(ss.str(), ss_mv.str(), y0);
  heat.writeResult(t0);
  { MyTimer timer("mrsdc-solve");
  sdc.solve(heat, y0, t0, tend, nStep);
  }
  heat.writeResult(tend);
}

template<int M>
void solve(HeatCoupled& heat, unsigned k_iter, double thetaFast, unsigned nStep, unsigned P)
{
 switch(P) 
 {
   case 1:
     solve<M,1>(heat, k_iter, thetaFast, nStep);
     break;
   case 2:
     solve<M,2>(heat, k_iter, thetaFast, nStep);
     break;
   case 3:
     solve<M,3>(heat, k_iter, thetaFast, nStep);
     break;
   case 4:
     solve<M,4>(heat, k_iter, thetaFast, nStep);
     break;
   case 5:
     solve<M,5>(heat, k_iter, thetaFast, nStep);
     break;
   case 6:
     solve<M,6>(heat, k_iter, thetaFast, nStep);
     break;
   case 7:
     solve<M,7>(heat, k_iter, thetaFast, nStep);
     break;
   case 8:
     solve<M,8>(heat, k_iter, thetaFast, nStep);
     break;
   default:
     throw std::runtime_error("not supported");
 } 
}

void solve(HeatCoupled& heat, unsigned k_iter, double thetaFast, unsigned nStep, unsigned M, unsigned P)
{
  switch(M) {
   case 1:
     solve<1>(heat, k_iter, thetaFast, nStep, P);
     break;
   case 2:
     solve<2>(heat, k_iter, thetaFast, nStep, P);
     break;
   case 3:
     solve<3>(heat, k_iter, thetaFast, nStep, P);
     break;
   case 4:
     solve<4>(heat, k_iter, thetaFast, nStep, P);
     break;
   case 5:
     solve<5>(heat, k_iter, thetaFast, nStep, P);
     break;
   case 6:
     solve<6>(heat, k_iter, thetaFast, nStep, P);
     break;
   case 7:
     solve<7>(heat, k_iter, thetaFast, nStep, P);
     break;
   case 8:
     solve<8>(heat, k_iter, thetaFast, nStep, P);
     break;
  }
}

int main(int argc, char* argv[])
{
  if(argc < 6) {
    std::cerr << "usage: " << argv[0] << " <nStep> <M> <P> <kIter> <laplTilde>\n";
    return 1;
  }
  //mpi-helper from dune
  MPIHelper::instance(argc, argv);
  int laplTilde(0);
  {std::stringstream ss ; ss << argv[5] ; ss >> laplTilde; }
  double nu(1.0e-3), alpha(1.0e-3), v0(5.0), src(0.0);
  bool useLapl0(false), addConstRobin(false);
  HeatCoupled heat(20, nu, alpha, v0, src, useLapl0, addConstRobin, laplTilde);
  //heat.setbAlph(1e-2);

  double thetaFast(1.0);
  unsigned nStep(40), kIter(2);
  {std::stringstream ss; ss << argv[1]; ss >> nStep; }
  {std::stringstream ss; ss << argv[4]; ss >> kIter; }
  unsigned M, P;
  {std::stringstream ss; ss << argv[2]; ss >> M; }
  {std::stringstream ss; ss << argv[3]; ss >> P; }

  solve(heat, kIter, thetaFast, nStep, M,P);
  return 0;
}

