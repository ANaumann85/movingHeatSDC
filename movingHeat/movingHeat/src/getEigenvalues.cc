#include "Heat.h"
#include <string>
#include <dune/istl/matrixmarket.hh>
#include <stdio.h>
#include <fstream>

template<typename M >
void writeFile(const M& m, std::string fname)
{
	std::fstream file(fname, std::ios_base::trunc | std::ios_base::out);
	file.precision(15);
	writeMatrixMarket(m, file);
	file.close();
}
int main(int argc, char* argv[])
{
	MPIHelper::instance(argc, argv);
	Heat heat(10, 1.0e-1, 1.0);
	const Heat::MatrixType& L(heat.getLaplacian());
	const Heat::MatrixType& M(heat.getMass());
	writeFile(L, "laplace.mm");
	writeFile(M, "mass.mm");
	unsigned nT(5);
	double te=15;
	double dt=te/nT;
	for(unsigned i(0); i < nT; ++i) {
		double t=i*dt;
		Heat::MatrixType lmv(heat.getLaplaceWithMove(t));
		char* strData;
		asprintf(&strData, "lmv-%0.2f.mm", t);
		std::string ss(strData);
		writeFile(lmv, ss);
		free(strData);
	}
	return 0;
}
