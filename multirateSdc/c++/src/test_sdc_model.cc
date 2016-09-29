#include <iostream>
using namespace std;
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

	cout.precision(15); cout << scientific;
	unsigned kIter(6);
	typedef MRSdc<Model::Vec, 3, 2> Method;
	Method sdc([](Model::Vec& ) {}, kIter, "radau_right", "radau_right");
	sdc.solve(problem, u0, t0, te, nStep);
	for(auto& d:u0) printf(" %25.20f",d); cout << endl;
	//array<double, 5 > pySol({0.99207383188965780896,   -0.00360178492349851097 ,  0.00044055218247540922 ,  -0.00319568211431540657 ,  -0.00128755226472033351});
	array<double, 5 > pySol({0.99207383188965780896,   -0.00360030904090934351 ,  0.00044316085811095860 ,  -0.00319581108458978722 ,  -0.00128580982270810971});

	double errOld(abs(pySol[0]-u0[0])/abs(pySol[0]));
	cout << errOld;
	for(unsigned i(1); i < 5; ++i) {
		double relErr=abs(pySol[i]-u0[i])/abs(pySol[0]);
		cout << " " << relErr;
		errOld = max(errOld, relErr);
	}
	cout << endl;

	return errOld < 6e-16 ? 0 : 1;
}

