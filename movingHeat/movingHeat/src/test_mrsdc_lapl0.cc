#include "Heat.h"
#include "MRSdc.h"
#include <sstream>
#include <fstream>
#include <algorithm>
#include "Timer.h"
#include "test_mrsdc_lapl0/refSol.h"

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
Heat::VectorType solve(Heat& heat, unsigned k_iter, unsigned nStep)
{
  typedef MRSdc<Heat::VectorType, M,P> Method;
  typename Method::Init init([&heat](Heat::VectorType& d) { heat.init(d); });
  Method sdc(init, k_iter, "radau_right", "radau_right");
  Heat::VectorType y0;
  //heat.init(y0, [](auto x) { return (x[0]-0.5)*(x[0]-0.5)*(x[1]-2)*(x[1]-2); });
  heat.init(y0); y0 = 0.0;

  double t0(0.0), tend(20.0);

  std::stringstream ss; ss << "heat_sdc_M-" << M << "_P-" << P << "_nStep-" << nStep;
  //heat.startFile(ss.str(), y0);
  //heat.writeResult(t0);
  { MyTimer timer("mrsdc-solve");
  sdc.solve(heat, y0, t0, tend, nStep);
  }
  //heat.writeResult(tend);
  return y0;
}

using namespace std;

void printSolution(Heat::VectorType& y, string fname)
{
  fstream file(fname, fstream::trunc | fstream::out);
  file.precision(18);
  for(auto& d: y)
    file << d << endl;
  file.close();
}

void getRefSol(Heat::VectorType& ref)
{
  char* asPtr = new char[sizeof(char)*(refSol_dat_len+1)];
  copy(refSol_dat, refSol_dat+refSol_dat_len, asPtr);
  std::string asStr(asPtr, refSol_dat_len);
  delete [] asPtr;
  stringstream ss; ss << asStr;
  double val;
  ss >> val;
  unsigned p(0);
  while(ss.good() && p < ref.size()) {
    ref[p] = val;
    ss >> val;
    ++p;
  }
}

int main(int argc, char* argv[])
{
  //mpi-helper from dune
  MPIHelper::instance(argc, argv);
  Heat heat(10, 1.0e-3, 1.0e-3, 5.0, 0.0, true); //true=lapl0
  auto res=solve<3,8>(heat, 2, 10);
  //printSolution(res, "refSol.dat");
  Heat::VectorType ref(res.size());
  getRefSol(ref);
  ref -= res;
  double madiff=ref.infinity_norm();
  cout.precision(18);
  cout << "max diff: " << madiff << std::endl;
  return madiff==0.0 ? 0 :1;
}
