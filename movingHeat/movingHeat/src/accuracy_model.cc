#include "Model.h"
#include "MRSdc.h"
#include "Ros2.h"
#include "Sdc.h"
#include <fstream>
#include <vector>
#include <stdexcept>
#include <cmath>

#include "linspace.h"
//#include "ColMat.h"

template<typename VT >
struct Matrix
{
	typedef std::vector<std::vector< VT > > Data;

	Data data;
	Matrix(unsigned nR, unsigned nC):
		data(nR)
	{
		for(auto& d : data)
			d.resize(nC);
	}

	VT& operator()(unsigned r, unsigned c)
	{ return data[r][c]; }

	const VT& operator()(unsigned r, unsigned c) const
	{ return data[r][c]; }
};

template<typename VT >
void operator<<(std::fstream& f, const Matrix<VT>& mat)
{
	f.precision(15);
	for(unsigned r(0); r < mat.data.size(); ++r) {
			f<< mat(r,0);
		for(unsigned c(1); c < mat.data[0].size(); ++c)
			f<< " " << mat(r,c);
		f << std::endl;
	}
}

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
	ostream& operator<<(ostream& o, const array<double, s >& l)
	{ 
		o << "{"; for(const auto& d: l) o << " " << d; o << "}";
		return o;
	}
}

static const int M =3;
static const int P =2;
int main(int argc, char* argv[])
{
	unsigned nStepSDC(1), nIter(1);
	double thetaFast(1.0);

  auto alpha_vec=linspace(0.0, 30, 60);
  auto nu_vec=linspace(0.0, 100, 200);

	double tend(1.0);
	unsigned nStepRef(100);
	Model::VectorType u0;
	for(auto& d:u0) d=0.0;
	u0[0]=1.0;

	std::string nodeName("radau_right");
	Matrix<Model::VectorType> refSol(alpha_vec.size(), nu_vec.size());
	Matrix<Model::VectorType> mrsdcSol(alpha_vec.size(), nu_vec.size());
	Matrix<Model::VectorType> sdcSol(alpha_vec.size(), nu_vec.size());
	Matrix<double > mrsdcAcc(alpha_vec.size(), nu_vec.size()), sdcAcc(alpha_vec.size(), nu_vec.size());
	Matrix<double > relSdcMrsdc(alpha_vec.size(), nu_vec.size());

	for(unsigned i(0); i < alpha_vec.size(); ++i) {
		for(unsigned j(0); j < nu_vec.size(); ++j) {
			Model model(1.0, nu_vec[j], alpha_vec[i], 5.0);
			auto initFu = [&model](auto& d) {model.init(d); for(auto& v:d) v=0.0; };

			Ros2<Model::VectorType > ros2(initFu);
			Model::VectorType uRos2; uRos2=u0;
			ros2.solve(model, uRos2, 0.0, tend, nStepRef);
			refSol(i,j) = uRos2;

			MRSdc<Model::VectorType, M, P > mrsdc(initFu, nIter, nodeName, nodeName, thetaFast);
			Model::VectorType uMRSDC; uMRSDC =u0;
			mrsdc.solve(model, uMRSDC, 0.0, tend, nStepSDC);
			mrsdcSol(i,j) = uMRSDC;

			Sdc<Model::VectorType, M > sdc(initFu, nIter, nodeName);
			Model::VectorType uSDC; uSDC =u0;
			sdc.solve(model, uSDC, 0.0, tend, nStepSDC);
			sdcSol(i,j) = uSDC;

			axpy(-1.0, uRos2, uMRSDC);
			mrsdcAcc(i,j) = norm(uMRSDC);
			axpy(-1.0, uRos2, uSDC);
			sdcAcc(i,j) = norm(uSDC);
			relSdcMrsdc(i,j) = sdcAcc(i,j)/mrsdcAcc(i,j);
		}
	}
	{
		std::fstream sdcFile("acc_sdc.dat", std::fstream::out);
		sdcFile << sdcAcc;
	}
	{
		std::fstream sdcFile("acc_mrsdc.dat", std::fstream::out);
		sdcFile << mrsdcAcc;
	}
	{
		std::fstream sdcFile("relSdcMrsdc.dat", std::fstream::out);
		sdcFile << relSdcMrsdc;
	}
#if 0
	{
		std::fstream sdcFile("sol_ros2.dat", std::fstream::out);
		sdcFile << refSol;
	}
#endif
	
  return 0;
}
