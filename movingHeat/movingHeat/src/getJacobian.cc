#include "Heat.h"
#include <fstream>
#include <vector>
#include <complex>

namespace Dune
{
  template< typename A, typename B >
  inline void axpy(double a, const BlockVector<A, B>& x, BlockVector<A, B>& y)
  { y.axpy(a, x); }
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
//typedef std::vector<std::complex<double> > VCD;
typedef std::vector<double > VCD;
VCD getEig(ColMat& cm)
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

#if 0
	if(0){std::fstream eigs("eigs.mtx",  std::ios_base::out | std::ios_base::trunc);
		eigs.precision(15);
	for(auto& d:lr) eigs << d << std::endl;}
#endif

	VCD ret(N);
	for(unsigned i(0); i < N; ++i) {
		//ret[i] = std::complex<double>(lr[i], li[i]);
		ret[i] = lr[i];
	}
	return ret;
}

int main(int argc, char* argv[])
{
	/*if(argc < 2) {
		std::cerr << "usage: " << argv[0] << " <Values>\n";
		return 1;
	}*/
	//mpi-helper from dune
	MPIHelper::instance(argc, argv);
	Heat heat(10);
	Heat::VectorType y0, fVal, fVal0;
	heat.init(y0); y0 = 0.0;
	heat.init(fVal); fVal=0.0;
	heat.init(fVal0); fVal0=0.0;
	double al=1e-2;
	double nu=1e-1;
	double delta=1e-8;
	heat.setParam(nu, al);
	ColMat colMat(y0.size(), y0.size());

	heat(0.0, y0, fVal0);
	for(unsigned i(0); i < y0.size(); ++i) {
		y0[i]=delta;
		heat(0.0, y0, fVal);
		fVal -= fVal0;
		fVal *= 1.0/delta;
		colMat.setColumn(i, fVal);
		y0[i]=0.0;
	}

	VCD eigs(getEig(colMat));

	std::fstream eigOut("eigsHeat.mtx", std::ios_base::out | std::ios_base::trunc);
	eigOut.precision(15);
	for(auto& d : eigs)
		eigOut << d << std::endl;
	eigOut.close();

	return 0;
}



