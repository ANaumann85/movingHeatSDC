#ifndef MODEL_H
#define MODEL_H

#include <array>
#include <cmath>

using namespace std;

namespace std
{
	template<unsigned long s>
	array<double, s> operator*(const array<array<double, s>, s>& m, const array<double, s>& v)
	{
		array<double, s> ret;
		for(unsigned i(0); i < s; ++i) {
			ret[i]=m[i][0]*v[0];
			for(unsigned j(1); j < s; ++j) {
				ret[i] += m[i][j]*v[j];
			}
		}
		return ret;
	}

	template<unsigned long s >
	array<double, s> operator+=(array<double, s>& d, const array<double,s>& r)
	{
		for(unsigned i(0); i < s; ++i)
			d[i]+=r[i];
		return d;
	}

}

struct Model
{
	typedef std::array<double, 5> VectorType;
	typedef std::array<VectorType, 5> MT;
	typedef VectorType A;
	typedef VectorType Vec;

	static constexpr unsigned dim = 5;
	//static constexpr double pi = 3.141592653589793; //wiki
	static constexpr double pi = 3.14159265359; //numpy 

	double a,nu,alpha, v0;
	A sd, fac;

	MT getFastMat(double t) const
	{
		MT ret;
		ret[0][0] = 0.5;
		for(unsigned i(1); i < 5; ++i)
			ret[0][i] = 0.0;
		// equation for a_1
		ret[1][0] = cos(a*t);
		ret[1][1] = cos(2.0*a*t)+1.0;
		ret[1][2] = cos(a*t);
		ret[1][3] = sin(2.0*a*t);
		ret[1][4] = sin(a*t);

		// equation for a_2
		ret[2][0] = cos(2.0*a*t);
		ret[2][1] = cos(a*t);
		ret[2][2] = 1.0;
		ret[2][3] = -sin(a*t);
		ret[2][4] = 0.0;

		// equation for b_1
		ret[3][0] = sin(a*t);
		ret[3][1] = sin(2.0*a*t);
		ret[3][2] = -sin(a*t);
		ret[3][3] = -cos(2.0*a*t)+1;
		ret[3][4] = cos(a*t);

		// equation for b_2
		ret[4][0] = sin(2.0*a*t);
		ret[4][1] = sin(a*t);
		ret[4][2] = 0.0;
		ret[4][3] = cos(a*t);
		ret[4][4] = 1.0;

		for(auto& d:ret)
			for(auto& e:d)
				e*= -alpha/(2*pi);
		return ret;
	}

	A getFastB(double t) const
	{
		A ret({0.5, cos(a*t), cos(2*a*t), sin(a*t), sin(2*a*t)});
		for(auto& d:ret) d *= alpha*v0/pi;
		return ret;
	}

	Model(double a, double nu, double al, double v0):
		a(a), nu(nu), alpha(al), v0(v0),
		sd({0.0,1.0,4.0,1.0,4.0})
	{
		for(auto& d:sd) d*= -nu;
	}

	void setParam(double nu, double al)
	{
		alpha=al;
		sd={0.0,1.0,4.0,1.0,4.0};
		this->nu=-nu;
		for(auto& d:sd) d*= -nu;
	}
	
	inline void init(A& ) const {}
	void slow(double, const A& in, A& out) const
	{
		for(unsigned i(0); i < dim; ++i)
			out[i] = in[i]*sd[i];
	}

	void fast(double t, const A& in, A& out) const
	{
		//A part
		auto fm=getFastMat(t);
		/*std::cout << "fm: \n[" ;
	       	for(const auto& d: getFastMat(t)) {
			std::cout << "[";
		       for(const auto& d2:d)
		             std::cout << " " << d2; 
		       std::cout << "]\n";
		}
		std::cout << "]\n";*/
		out = fm*in;
		//b part
		out += getFastB(t);
		//std::cout << "fast: [" ; for(const auto& d: getFastB(t)) std::cout << " " << d; std::cout << "]\n";
	}

	void Mv(const A& v, A& out) const
	{ out = v; }

	void MinvV(const A& v, A& out) const
	{ out =v; }

	void updateMatrix(double t, double a)
	{ 
		for(unsigned i(0); i < dim; ++i)
		 fac[i] = 1.0/(1-a*sd[i]);
	}

	void solveMaJ(const A& r, A& x) const
	{
		for(unsigned i(0); i < dim; ++i)
			x[i] = fac[i]*r[i];
	}
};
#endif
