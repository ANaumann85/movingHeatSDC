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
}
using namespace std;
int main(int argc, char* argv[])
{
	double l1(-0.1), l2(-0.0);
	double tstart = 0.0;
	double tend   = 10.0;
	Problem::Vec u0({2.0});
	double u_ex  = u0[0]*exp(tend*(l1+l2));

	Problem problem(l1, l2);

	MRSdc<Problem::Vec, 2, 2> sdc;
	double te=1.0;
	double t0=tstart;
	sdc.predict(problem, u0, t0, te);
	sdc.sweep(problem, u0, t0, te);
	return 0;
}
