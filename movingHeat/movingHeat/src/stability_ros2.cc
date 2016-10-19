#include "Heat.h"
#include "Ros2.h"
#include <fstream>
#include <vector>

namespace Dune
{
  template< typename A, typename B >
  inline void axpy(double a, const BlockVector<A, B>& x, BlockVector<A, B>& y)
  { y.axpy(a, x); }
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
  double al = 1e-5, nu=1e-3;
  Heat heat(10, nu, al, 0.0, 0.0);
  Ros2<Heat::VectorType > ros2([&heat](auto& in) { heat.init(in); in=0.0; });
  Heat::VectorType y0;
  heat.init(y0); y0 = 0.0;

  unsigned nTests(20);
  auto alpha_vec=linspace(0.0, 15.0, nTests);
  auto nu_vec=linspace(0.0, 2, nTests); //nTests);
  ColMat colMat(y0.size(), y0.size());
  ColMat maxEigs(alpha_vec.size()+1, nu_vec.size()+1);
#if 1
  for(unsigned aln(0); aln < alpha_vec.size(); ++aln) {
    const double al = alpha_vec[aln];
    maxEigs(aln+1,0)=al;
    for(unsigned nun(0); nun < nu_vec.size(); ++nun) {
      const double nu=nu_vec[nun];
      maxEigs(0,nun+1)=nu;
      heat.setParam(nu, al);
#endif
      for(unsigned i(0); i < y0.size(); ++i) {
        y0[i]=1.0;
        ros2.solve(heat, y0, 0.0, 1.0, 1);
        colMat.setColumn(i, y0);
        y0=0.0;
      }
#if 1
      maxEigs(aln+1,nun+1)=getMaxAbsEig(colMat);
      //std::cout << "maxEig:" << maxEigs(aln+1, nun+1) << std::endl;

    }
  }
#endif

#if 0
  std::fstream mat("mat.mtx", std::ios_base::out | std::ios_base::trunc);
  mat << colMat;
  mat.close();
  std::cout << "maxEig:" <<getMaxAbsEig(colMat)  << std::endl;
#endif
#if 1
  std::fstream eigOut("stabEigs_ros2.mtx", std::ios_base::out | std::ios_base::trunc);
  eigOut << maxEigs;
  eigOut.close();
#else
  std::cout << maxEigs << std::endl;
#endif

  return 0;
}


