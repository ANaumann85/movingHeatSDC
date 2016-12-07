#include <cmath>
#include "MRSdc.h"
#include "Sdc.h"
#include <array>

struct ConstSource
{
	double operator()(double ) const
	{ return 2.0; }

	double uex(double t, double nu, double u0) const
	{ return u0*exp(nu*t) + 2.0/nu*(exp(nu*t)-1.0); }
};

struct SineSource
{
	double operator()(double t) const
	{ return sin(t); }

	double uex(double t, double nu, double u0) const
	{
		double c1=u0 + 1.0/(nu*nu+1);
		return c1*exp(nu*t) - (nu*sin(t) + cos(t))/(nu*nu+1);
	}
};

template<bool useSlowExpl, typename Source >
struct Problem
{
	typedef std::array<double, 1u > Vec;
	double nu;
	double fac;
	Source source;

	Problem(double nu):
		nu(nu)
	{}

	inline void fast(double t, const Vec& in, Vec& out)
	{ out[0] =  useSlowExpl ? 0.0 : source(t); }

	inline void slow(double t, const Vec& in, Vec& out)
	{ out[0] = nu*in[0] ; if(useSlowExpl) out[0] += source(t); }

	inline void slowExpl(double t, const Vec& in, Vec& out) const
	{ out[0] = useSlowExpl? source(t) : 0.0; } 

	inline void slowImpl(double t, const Vec& in, Vec& out) const
	{ out[0] = nu*in[0]; }

	inline void slowSrc(double t, Vec& dest) const
	{ dest[0] = 0.0; }

	inline void updateMatrix(double t, double a)
	{ fac = 1.0/(1.0-a*nu); }

	inline void solveMaJ(const Vec& in, Vec& out)
	{ out[0] = fac*in[0]; }

	inline void Mv(const Vec& in, Vec& out) const
	{ out[0] = in[0]; }

	inline void MinvV(const Vec& in, Vec& out) const
	{ out[0] = in[0]; }

	double uex(double t, double u0)
	{ return source.uex(t, nu, u0);	}
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
	double te   = 3.0;
	const unsigned nStep(10);

	typedef Problem<false, ConstSource>::Vec Vec;
	typedef MRSdc<Vec, 3, 8> MRSDC;
	typedef Sdc<Vec, 3> SDC;
	Vec u0({1.0});

	Problem<false, SineSource> problem(nu);
	Problem<true, SineSource> problem2(nu);
	double u_ex  = problem.uex(te, u0[0]);

	for(unsigned kIter(0); kIter < 5; ++kIter) {

		MRSDC mrsdc([](Vec& ) {}, kIter, "radau_right", "radau_right", 1.0);
		SDC sdc([](Vec& ) {}, kIter, "radau_right");

		Vec u0M, u0S;
		u0M[0]=u0[0]; u0S[0] = u0[0];
		mrsdc.solve(problem2, u0M, t0, te, nStep);
		sdc.solve(problem, u0S, t0, te, nStep);
		std::cout << "diff(" << kIter << "): " << abs(u0S[0]-u0M[0]) << " " << abs(u_ex-u0S[0]) << " " << abs(u_ex-u0M[0]) <<   " " << u_ex << std::endl;
	}
	return 0;
}
