#include <iostream>
#include "Model.h"
#include "MRSdc.h"
#include <fstream>
#include <vector>
#include <stdexcept>
#include <cmath>

#include "linspace.h"
#include "ColMat.h"

namespace std
{
	template<typename T, unsigned long s >
	void axpy(double fac, const array<T, s>& x, array<T, s>& y)
	{
		for(unsigned i(0); i < s; ++i)
		y[i] += fac*x[i];
	}

	template<typename T, unsigned long s >
	void setValue(array<T, s>& dest, const T& v)
	{ for(T& d:dest) d= v; }

	template<unsigned long s >
	double norm(const array<double, s >& a)
	{
		double ret(0.0);
		for(const auto& d:a) ret =max(ret, abs(d));
		return ret;
	}

	template<unsigned long s>
	array<double, s > operator-(const array<double, s >& l, const array<double, s >& r)
	{ array<double, s> ret; for(unsigned i(0); i < s ; ++i) ret[i] = l[i]-r[i]; return ret; }
}


int main(int argc, char* argv[])
{
  /*if(argc < 2) {
    std::cerr << "usage: " << argv[0] << " <Values>\n";
    return 1;
    }*/
  //mpi-helper from dune
  double al = 0.0, nu=3e-0;
  Model heat(1.0, nu, al, 0.0);
  typedef MRSdc<Model::VectorType, 3,2> Method;
  Method::Init init([&heat](Model::VectorType& d) { heat.init(d); for(auto& e:d) e=0.0; });
  Method sdc(init, 6, "radau_right", "radau_right", 1.0);
  Model::VectorType y0;
  heat.init(y0); 
  for(auto& d:y0) d=0.0;

  unsigned nTests(20);
  auto alpha_vec=linspace(0.0, 35, nTests);
  auto nu_vec=linspace(0.0, 5, nTests);
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
#else
    bool hasSpecial(false);
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
          for(auto& d:y0) d=0.0;
					//for(auto& d:y0) std::cout << " " << d; cout << endl;
        }
#if 1
      if(hasSpecial)
        maxEigs(aln+1,nun+1)=10;
      else
        maxEigs(aln+1,nun+1)=getMaxAbsEig(colMat);
      //std::cout << "maxEig: " << al << " " << nu << " " << maxEigs(aln+1, nun+1) << std::endl;

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
