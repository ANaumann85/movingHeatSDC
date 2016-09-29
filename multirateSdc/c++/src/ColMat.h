#include <complex>
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
	if(info != 0)
		throw std::runtime_error("error in dgeev");

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
