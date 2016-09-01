#include "MRSdc.h"
#include <array>
#include <cmath>

struct Problem
{
	typedef std::array<double, 1u > Vec;
	double l1, l2;
	double fac;

	Problem(double l1, double l2):
		l1(l1), l2(l2)
	{}

	inline void fast(double t, const Vec& in, Vec& out)
	{ out[0] = l2*in[0]; }

	inline void slow(double t, const Vec& in, Vec& out)
	{ out[0] = l1*in[0]; }

	inline void updateMatrix(double t, double a)
	{ fac = 1.0/(1.0-a*l1); }

	inline void solveMaJ(const Vec& in, Vec& out)
	{ out[0] = fac*in[0]; }
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
	double l1(-0.1), l2(-1.0);
	double t0 = 0.0;
	double te   = 10.0;
	unsigned nStep(1);
	Problem::Vec u0({2.0});
	double u_ex  = u0[0]*exp((te-t0)*(l1+l2));

	Problem problem(l1, l2);

	unsigned kIter(12);
	typedef MRSdc<Problem::Vec, 5, 4> Method;
	Method sdc(kIter, "radau_right", "equi_noleft");
	std::cout.precision(8);
#if 0
	double dt=(te-t0)/nStep;
	for(unsigned s(0); s < nStep; ++s) {
		double t0_=t0+s*dt; 
		double te= t0_+dt;
		sdc.predict(problem, u0, t0_, te);
		for(unsigned k(0); k < kIter; ++k) {
			sdc.sweep(problem, u0, t0_, te);
		}
		cout << "standard residual: " << sdc.residual(u0) << std::endl;
		cout << "embedded residual: " << sdc.sub_residual(u0) << std::endl;
		u0[0] = sdc.us[Method::M-1][0];
	}
	cout << "error:" << abs(u0[0]-u_ex)/abs(u_ex) << endl;
#else
	unsigned nTest(10);
	for(unsigned n(0); n < nTest; ++n, nStep *=2) {
		u0[0] = 2.0;
		sdc.solve(problem, u0, t0, te, nStep);
		cout << "error(" << nStep << "):" << abs(u0[0]-u_ex)/abs(u_ex) << endl;
	}
#endif
	return 0;
}
