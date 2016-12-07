#include "MRSdc.h"
#include <array>
#include <cmath>

struct Problem
{
	typedef std::array<double, 1u > Vec;
	double nu;
	double fac;

	Problem(double nu):
		nu(nu)
	{}

	inline void fast(double t, const Vec& in, Vec& out) const
	{ out[0] = sin(t); }

  inline void slowSrc(double, Vec& out) const
  { out[0] =0.0; }

  inline void slowExpl(double, const Vec& , Vec& out) const
  { out[0] = 0.0; }

  inline void slowImpl(double, const Vec& in, Vec& out) const
  { out[0] = nu*in[0]; }

	inline void slow(double t, const Vec& in, Vec& out) const
	{ out[0] = nu*in[0]; }

	inline void updateMatrix(double t, double a)
	{ fac = 1.0/(1.0-a*nu); }

	inline void solveMaJ(const Vec& in, Vec& out) const
	{ out[0] = fac*in[0]; }

	inline void Mv(const Vec& in, Vec& out) const
	{ out[0] = in[0]; }

	inline void MinvV(const Vec& in, Vec& out) const
	{ out[0] = in[0]; }

	double uex(double t, double u0)
	{
		double c1=u0 + 1.0/(nu*nu+1);
		return c1*exp(nu*t) - (nu*sin(t) + cos(t))/(nu*nu+1);
	}
};

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

using namespace std;
int main(int argc, char* argv[])
{
	double nu(-1.0);
	double t0 = 0.0;
	double te   = 0.5;
	unsigned nStep(1);
	unsigned nTest(10);
	Problem::Vec u0({1.0});

	Problem problem(nu);
	double u_ex  = problem.uex(te, u0[0]);

	unsigned kIter(15);
	typedef MRSdc<Problem::Vec, 5, 5> Method;
	Method sdc([](Problem::Vec& ) {}, kIter, "radau_right", "radau_right", 1.0);
	std::cout.precision(8);
	double errOld;
	sdc.solve(problem, u0, t0, te, nStep);
	errOld = abs(u0[0]-u_ex);
	cout << "error(" << nStep << "): " << u0[0] << " " << u_ex << " " << errOld << " " << endl;///abs(u_ex)
	nStep *=2;
	for( ; nTest > 0; --nTest, nStep *=2) {
		u0[0] = 1.0;
		sdc.solve(problem, u0, t0, te, nStep);
		double errNew = abs(u0[0]-u_ex);
		cout << "error(" << nStep << "): " << u0[0] << " " << u_ex << " " << errOld << " " << log(errOld/errNew)/log(2) << endl;///abs(u_ex)
		errOld = errNew;
	}
	if(errOld > 3e-15)
		return 1;
	return 0;
}
