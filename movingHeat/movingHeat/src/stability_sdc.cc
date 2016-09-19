#include "Heat.h"
#include "MRSdc.h"
#include <fstream>
#include <vector>

namespace Dune
{
  template< typename A, typename B >
  inline void axpy(double a, const BlockVector<A, B>& x, BlockVector<A, B>& y)
  { y.axpy(a, x); }

  template< typename A, typename B >
  inline void setValue(BlockVector<A, B>& x, const double& v)
  { x=v; }
}

std::vector<double > linspace(double v0, double ve, unsigned nValues)
{
	if(nValues < 2)
		nValues == 2;
	double dv=(ve-v0)/(nValues-1);
	std::vector<double > ret(nValues);
	ret[0]=v0;
	for(unsigned i(1); i < nValues; ++i)
		ret[i] = v0+i*dv;
	return ret;
}

struct ColMat
{
	typedef std::vector<double > V;
	V data;
	unsigned nr, nc;
	ColMat(unsigned nR, unsigned nC):
		data(nR*nC),
		nr(nR), nc(nC)
	{}

	template<typename VIn >
	void setColumn(unsigned c, const VIn& v)
	{
		unsigned s(c*nr);
		for(unsigned i(0); i < nr; ++i, ++s) {
			data[s] = v[i];
		}
	}

	inline const double& operator()(unsigned r, unsigned c) const
	{ return data[c*nr+r]; }

	inline double& operator()(unsigned r, unsigned c)
	{ return data[c*nr+r]; }

	double* getData() 
	{ return data.data(); }
};

void operator<<(std::fstream& f, const ColMat& mat)
{
	f.precision(15);
	for(unsigned r(0); r < mat.nr; ++r) {
			f<< mat(r,0);
		for(unsigned c(1); c < mat.nc; ++c)
			f<< " " << mat(r,c);
		f << std::endl;
	}
}

extern "C" {
	void dgeev_(const char* JOBVL, const char* JOBVR, const int* N,
                       const double* A, const int* LDA, double* WR,  double* WI,
                       double* VL, const int* LDVL, double* VR, const int* LDVR,
                       double* WORK, const int* LWORK, int* INFO);
}
double getMaxAbsEig(ColMat& cm)
{
	char job('N');
	int N(cm.nr), one(1), lwork(-1), info(0);	
	std::vector<double > lr(N), li(N), work(1);

	double vl;
	
#if 0
	{std::fstream mat("mat_in.mtx", std::ios_base::out | std::ios_base::trunc);
	mat << cm;
	mat.close();}
#endif
	//request memory
	dgeev_(&job, &job, &N, cm.getData(), &N, lr.data(), li.data(), &vl, &one, &vl, &one, work.data(), 
			&lwork, &info);

	lwork = (int) work[0];
	work.resize(lwork);
	dgeev_(&job, &job, &N, cm.getData(), &N, lr.data(), li.data(), &vl, &one, &vl, &one, work.data(), 
			&lwork, &info);

#if 1
	{std::fstream eigs("eigs.mtx",  std::ios_base::out | std::ios_base::trunc);
		eigs.precision(15);
	//for(auto& d:lr) eigs << d << std::endl;}
	for(unsigned i(0); i < N; ++i) eigs << std::complex<double>(lr[i], li[i]) << std::endl;}
#endif

	vl = lr[0]*lr[0]+li[0]*li[0];
	for(unsigned i(1); i < N; ++i) {
		double cur = lr[i]*lr[i]+li[i]*li[i];
		if( cur > vl)
			vl=cur;
	}
	return sqrt(vl);
}

int main(int argc, char* argv[])
{
  /*if(argc < 2) {
    std::cerr << "usage: " << argv[0] << " <Values>\n";
    return 1;
    }*/
  //mpi-helper from dune
  MPIHelper::instance(argc, argv);
  double al = 1e-4, nu=1e-3;
  Heat heat(10, nu, al);
  typedef MRSdc<Heat::VectorType, 3,2> Method;
  Method::Init init([&heat](Heat::VectorType& d) { heat.init(d); d=0.0; });
  Method sdc(init, 5, "radau_right", "equi_noleft");
  Heat::VectorType y0;
  heat.init(y0); y0 = 0.0;

  unsigned nTests(2);
  auto alpha_vec=linspace(0.0, 0.001, nTests);
  auto nu_vec=linspace(0.0, 5.0, 2);
  ColMat colMat(y0.size(), y0.size());
  ColMat maxEigs(alpha_vec.size()+1, nu_vec.size()+1);
#if 0
  for(unsigned aln(0); aln < alpha_vec.size(); ++aln) {
    const double al = alpha_vec[aln];
    maxEigs(aln+1,0)=al;
    for(unsigned nun(0); nun < nu_vec.size(); ++nun) {
      const double nu=nu_vec[nun];
      maxEigs(0,nun+1)=nu;
#endif
      //heat.setParam(nu, al);
      for(unsigned i(0); i < y0.size(); ++i) {
        y0[i]=1.0;
        sdc.solve(heat, y0, 0.0, 1.0, 1);
        colMat.setColumn(i, y0);
        y0=0.0;
      }
#if 0
      maxEigs(aln+1,nun+1)=getMaxAbsEig(colMat);
      std::cout << "maxEig:" << maxEigs(aln+1, nun+1) << std::endl;

    }
  }
#endif

#if 1
  std::fstream mat("mat.mtx", std::ios_base::out | std::ios_base::trunc);
  mat << colMat;
  mat.close();
  std::cout << "maxEig:" <<getMaxAbsEig(colMat)  << std::endl;
#endif
#if 0
  std::fstream eigOut("stabEigs_sdc.mtx", std::ios_base::out | std::ios_base::trunc);
  eigOut << maxEigs;
  eigOut.close();
#endif

  return 0;
}


