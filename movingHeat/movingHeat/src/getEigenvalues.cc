#include "Heat.h"
#include <string>
#include <dune/istl/matrixmarket.hh>
#include <stdio.h>

int main(int argc, char* argv[])
{
	MPIHelper::instance(argc, argv);
	Heat heat(10, 1.0e-3, 1.0);
	const Heat::MatrixType& L(heat.getLaplacian());
	const Heat::MatrixType& M(heat.getMass());
	storeMatrixMarket(L, "laplace.mm");
	storeMatrixMarket(M, "mass.mm");
	unsigned nT(5);
	double te=15;
	double dt=te/nT;
	for(unsigned i(0); i < nT; ++i) {
		double t=i*dt;
		Heat::MatrixType lmv(heat.getLaplaceWithMove(t));
		char* strData;
		asprintf(&strData, "lmv-%0.2f.mm", t);
		std::string ss(strData);
		storeMatrixMarket(lmv, ss);
		free(strData);
	}
	return 0;
}
