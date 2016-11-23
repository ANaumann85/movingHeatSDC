#include "Model.h"
#include "MRSdc.h"
#include "Ros2.h"
#include "Sdc.h"
#include <fstream>
#include <vector>
#include <stdexcept>
#include <cmath>

#include "linspace.h"

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

double tend(1.0);
unsigned nIter(1), nStepSDC(1);
std::string nodeName("radau_right");

template<unsigned M, typename InitFu>
void solveSDC(Model& model, Model::VectorType& uSDC, const InitFu& initFu)
{
			Sdc<Model::VectorType, M > sdc(initFu, nIter, nodeName);
			sdc.solve(model, uSDC, 0.0, tend, nStepSDC);
}

template<typename InitFu>
void solveSDC(Model& model, Model::VectorType& uSDC, const InitFu& initFu, unsigned M)
{
	switch(M) {
		case 2:
			solveSDC<2>(model, uSDC, initFu);
			break;
		case 3:
			solveSDC<3>(model, uSDC, initFu);
			break;
		case 4:
			solveSDC<4>(model, uSDC, initFu);
			break;
		case 5:
			solveSDC<5>(model, uSDC, initFu);
			break;
		case 6:
			solveSDC<6>(model, uSDC, initFu);
			break;
		case 7:
			solveSDC<7>(model, uSDC, initFu);
			break;
		case 8:
			solveSDC<8>(model, uSDC, initFu);
			break;
	}
}

template<unsigned M, unsigned P, typename InitFu>
void solveMRSDC(Model& model, Model::VectorType& uSDC, const InitFu& initFu)
{
	typedef MRSdc<Model::VectorType, M, P > Method;
	double thetaFast(1.0);
	typename Method::Init subInit = initFu;
	//MRSdc<Model::VectorType, M, P > mrsdc(initFu, nIter, nodeName, nodeName, thetaFast);
	Method mrsdc(subInit, nIter, nodeName, nodeName, thetaFast);
	mrsdc.solve(model, uSDC, 0.0, tend, nStepSDC);
}

template<unsigned P, typename InitFu>
void solveMRSDC(Model& model, Model::VectorType& uSDC, const InitFu& initFu, unsigned M)
{
	switch(M) {
		case 2:
			solveMRSDC<2,P>(model, uSDC, initFu);
			break;
		case 3:
			solveMRSDC<3,P>(model, uSDC, initFu);
			break;
		case 4:
			solveMRSDC<4,P>(model, uSDC, initFu);
			break;
		case 5:
			solveMRSDC<5,P>(model, uSDC, initFu);
			break;
		case 6:
			solveMRSDC<6,P>(model, uSDC, initFu);
			break;
		case 7:
			solveMRSDC<7,P>(model, uSDC, initFu);
			break;
		case 8:
			solveMRSDC<8,P>(model, uSDC, initFu);
			break;
	}
}

template<typename InitFu>
void solveMRSDC(Model& model, Model::VectorType& uSDC, const InitFu& initFu, unsigned M, unsigned P)
{
	switch(P) {
		case 2:
			solveMRSDC<2>(model, uSDC, initFu, M);
			break;
		case 3:
			solveMRSDC<3>(model, uSDC, initFu, M);
			break;
		case 4:
			solveMRSDC<4>(model, uSDC, initFu, M);
			break;
		case 5:
			solveMRSDC<5>(model, uSDC, initFu, M);
			break;
		case 6:
			solveMRSDC<6>(model, uSDC, initFu, M);
			break;
		case 7:
			solveMRSDC<7>(model, uSDC, initFu, M);
			break;
		case 8:
			solveMRSDC<8>(model, uSDC, initFu, M);
			break;
	}
}
	
bool constant_jacobian(const Model& ) { return true; } 
int main(int argc, char* argv[])
{
	if(argc < 5) {
		std::cout << "usage: " << argv[0] << "<M> <P> <nStepSDC> <nIter>\n";
		return 1;
	}
	unsigned M(3), P(2);
	{ std::stringstream ss; ss << argv[1] ; ss >> M; }
	{ std::stringstream ss; ss << argv[2] ; ss >> P; }
	{ std::stringstream ss; ss << argv[3] ; ss >> nStepSDC; }
	{ std::stringstream ss; ss << argv[4] ; ss >> nIter; }
	

  auto alpha_vec=linspace(0.0, 30, 60);
  auto nu_vec=linspace(0.0, 100, 200);

	unsigned nStepRef(100);
	Model::VectorType u0;
	for(auto& d:u0) d=0.0;
	u0[0]=1.0;

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

			Model::VectorType uMRSDC; uMRSDC =u0;
			solveMRSDC(model, uMRSDC, initFu, M, P);
			mrsdcSol(i,j) = uMRSDC;

			Model::VectorType uSDC; uSDC =u0;
			solveSDC(model, uSDC, initFu, M);
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
