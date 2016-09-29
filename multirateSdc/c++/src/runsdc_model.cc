#include <iostream>
#include "Model.h"
#include "MRSdc.h"
#include <array>
#include <cmath>
#include "stdio.h"

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
	
	template<unsigned long s>
	ostream& operator<<(ostream& out, const array<double, s>& d)
	{
		out << "["; for(auto& e:d) out << " " << e; out << "]";
		return out;
	}
}

using namespace std;
int main(int argc, char* argv[])
{
	double nu(3.0);
	double al(0.1);
	double v0(0.0);
	double t0 = 0.0;
	double te   = 1.0;
	unsigned nStep(1);
	unsigned nTest(10);
	Model::Vec u0({1.0, 0.0, 0.0, 0.0, 0.0});
	Model::Vec f;

	Model problem(1.0, nu, al, v0);

#if 1
	unsigned kIter(6);
	typedef MRSdc<Model::Vec, 3, 2> Method;
	Method sdc([](Model::Vec& ) {}, kIter, "radau_right", "radau_right");
	std::cout.precision(20);
	double errOld;
	sdc.solve(problem, u0, t0, te, nStep);
	//cout << u0 <<endl;
	for(auto& d:u0) printf(" %25.20f",d); cout << endl;
#endif
#if 0
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
#endif
	return 0;
}
