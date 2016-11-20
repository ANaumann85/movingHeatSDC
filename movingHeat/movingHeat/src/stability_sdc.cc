#include "Heat.h"
#include "MRSdc.h"
#include <fstream>
#include <vector>
#include <stdexcept>
#include <cmath>

namespace Dune
{
  template< typename A, typename B >
  inline void axpy(double a, const BlockVector<A, B>& x, BlockVector<A, B>& y)
  { y.axpy(a, x); }

  template< typename A, typename B >
  inline void setValue(BlockVector<A, B>& x, const double& v)
  { x=v; }
}

#include "linspace.h"
#include "ColMat.h"

int main(int argc, char* argv[])
{
  /*if(argc < 2) {
    std::cerr << "usage: " << argv[0] << " <Values>\n";
    return 1;
    }*/
  //mpi-helper from dune
  MPIHelper::instance(argc, argv);
  double al = 1e-0, nu=1e-0;
  Heat heat(10, nu, al, 0.0);
  typedef MRSdc<Heat::VectorType, 3,2> Method;
  Method::Init init([&heat](Heat::VectorType& d) { heat.init(d); d=0.0; });
  Method sdc(init, 5, "radau_right", "radau_right", 1.0);
  Heat::VectorType y0;
  heat.init(y0); y0 = 0.0;

  unsigned nTests(20);
  auto alpha_vec=linspace(0.0, 1.0e-2, nTests);
  auto nu_vec=linspace(0.0, 3.0, nTests);
  //std::vector<double> nu_vec(1); nu_vec[0]=2.0;
  ColMat colMat(y0.size(), y0.size());
  ColMat maxEigs(alpha_vec.size()+1, nu_vec.size()+1);
#if 1
  for(unsigned nun(0); nun < nu_vec.size(); ++nun) {
    const double nu=nu_vec[nun];
    maxEigs(0,nun+1)=nu;
    bool hasSpecial(false);
    for(unsigned aln(0); aln < alpha_vec.size(); ++aln) {
      const double al = alpha_vec[aln];
      maxEigs(aln+1,0)=al;
      heat.setParam(nu, al);
#endif
      if(!hasSpecial) 
        for(unsigned i(0); i < y0.size(); ++i) {
          y0[i]=1.0;
          sdc.solve(heat, y0, 0.0, 1.0, 1);
          for(auto& d:y0)
            if(std::isnan(d) || std::isinf(d)) {
              hasSpecial=true;
              break;
            }
          if(hasSpecial)
            break;
          colMat.setColumn(i, y0);
          y0=0.0;
        }
#if 1
      if(hasSpecial)
        maxEigs(aln+1,nun+1)=10;
      else
        maxEigs(aln+1,nun+1)=getMaxAbsEig(colMat);
      std::cout << "maxEig: " << al << " " << nu << " " << maxEigs(aln+1, nun+1) << std::endl;

    }
  }
#endif

#if 0
  std::fstream mat("mat.mtx", std::ios_base::out | std::ios_base::trunc);
  mat << colMat;
  mat.close();
  std::cout << "maxEig:" <<getMaxAbsEig(colMat)  << std::endl;
#else
  std::fstream eigOut("stabEigs_sdc.mtx", std::ios_base::out | std::ios_base::trunc);
  eigOut << maxEigs;
  eigOut.close();
#endif

  return 0;
}


