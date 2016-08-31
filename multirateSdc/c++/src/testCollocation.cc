inline void setValue(double& d, double v) { d=v; }
inline void axpy(double a, double x, double& y) { y += a*x ;}

#include "MultirateCollocation.hh"
#include <array>
#include <cmath>
#include <iostream>

using namespace std;

template<unsigned Deg >
struct PolynomMonomBase
{
	static const unsigned degree = Deg;
	array<double, Deg+1> coeffs;

	//create zero polynom
	PolynomMonomBase()
	{
		for(double& d : coeffs)
			d= 0.0;
	}

	PolynomMonomBase(const array<double, Deg+1>& c):
		coeffs(c)
	{}

	double operator()(const double& t) const
	{
		double ret(coeffs[0]);
		for(unsigned i(1); i <= Deg; ++ i)
			ret = ret*t + coeffs[i];
		return ret;
	}

	double integrate(double t0, double t1) const
	{
		PolynomMonomBase<degree+1> stamm;
		for(unsigned i(0); i < degree; ++i) {
			stamm.coeffs[i] = coeffs[i]/(degree-i+1);
		}
		stamm.coeffs[degree] = coeffs[degree];
		stamm.coeffs[degree+1] = 0.0;
		return stamm(t1)-stamm(t0);
	}

};

//test standard node-to-node
template<typename Poly, unsigned M, unsigned P >
void test_m_mp1(const Poly& poly, MultirateCollocation<M, P>& coll)
{
	auto fVals(coll.coll.evalAtNodes(poly));
	//cerr << "fVals:["; for(auto& d:fVals) cerr << " " << d ; cerr << "]" << endl;
	for(unsigned m(0); m < M-1; ++m) {
		double ih;
	       	coll.integrate_m_mp1(fVals, m, ih);
		double ta = m == 0? coll.coll.tleft : coll.coll.nodes[m-1];
		double iex=poly.integrate(ta, coll.coll.nodes[m]);
		if(abs(ih-iex) > 1e-7) {
			std::cerr << "error in node-to-node:\n";
			cerr << "nodes m, node[m], node[m+1]: " << m << " " << coll.coll.nodes[m] << " " << coll.coll.nodes[m+1] << endl;
			cerr << "ih iex " << ih << " " << iex << endl;
			cerr << "polynom: [";	for(auto& d : poly.coeffs) cerr << " " << d ;	cerr << "]" << endl;
			cerr << "diff:" << abs(ih-iex) << endl;
			throw std::runtime_error("wrong node-to-node m_mp1 integration");
		}
	}
}

//test standard node-to-node integration on sub polynom
template<typename Poly, unsigned M, unsigned P >
void test_m_mp1_sub(const Poly& poly, MultirateCollocation<M, P>& coll)
{
	//cerr << "fVals:["; for(auto& d:fVals) cerr << " " << d ; cerr << "]" << endl;
	for(unsigned m(0); m < M-1; ++m) {
		auto fVals(coll.coll_sub[m].evalAtNodes(poly));
		double ih;
	       	coll.integrate_m_mp1_sub(fVals, m, ih);
		double ta = m == 0? coll.coll.tleft : coll.coll.nodes[m-1];
		double iex=poly.integrate(ta, coll.coll.nodes[m]);
		if(abs(ih-iex) > 1e-7) {
			std::cerr << "error in node-to-node sub:\n";
			cerr << "nodes m, node[m], node[m+1]: " << m << " " << coll.coll.nodes[m] << " " << coll.coll.nodes[m+1] << endl;
			cerr << "ih iex " << ih << " " << iex << endl;
			cerr << "polynom: [";	for(auto& d : poly.coeffs) cerr << " " << d ;	cerr << "]" << endl;
			cerr << "diff:" << abs(ih-iex) << endl;
			throw std::runtime_error("wrong node-to-node m_mp1_sub integration");
		}
	}
}

//test sub-node to sub-node integration on standard polynom
template<typename Poly, unsigned M, unsigned P >
void test_p_pp1(const Poly& poly, MultirateCollocation<M, P>& coll)
{
	auto fVals(coll.coll.evalAtNodes(poly));
	//cerr << "fVals:["; for(auto& d:fVals) cerr << " " << d ; cerr << "]" << endl;
	for(unsigned m(0); m < M; ++m) {
		double ih;
		for(unsigned p(0); p < P; ++p) {
			coll.integrate_p_pp1(fVals, m, p, ih);
			double ta = p == 0 ? coll.coll_sub[m].tleft : coll.coll_sub[m].nodes[p-1];
			double iex=poly.integrate(ta, coll.coll_sub[m].nodes[p]);
			if(abs(ih-iex) > 1e-7) {
				std::cerr << "error in node-to-node:\n";
				cerr << "m, p, node[p-1], node[p]: " << m << " " << p << " " << ta << " " << coll.coll_sub[m].nodes[p] << endl;
				cerr << "ih iex " << ih << " " << iex << endl;
				cerr << "polynom: [";	for(auto& d : poly.coeffs) cerr << " " << d ;	cerr << "]" << endl;
				cerr << "diff:" << abs(ih-iex) << endl;
				throw std::runtime_error("wrong node-to-node p_pp1 integration");
			}
		}
	}
}

//test sub-node to sub-node integration on standard polynom
template<typename Poly, unsigned M, unsigned P >
void test_p_pp1_sub(const Poly& poly, MultirateCollocation<M, P>& coll)
{
	//cerr << "fVals:["; for(auto& d:fVals) cerr << " " << d ; cerr << "]" << endl;
	for(unsigned m(0); m < M; ++m) {
		double ih;
		auto fVals(coll.coll_sub[m].evalAtNodes(poly));
		for(unsigned p(0); p < P; ++p) {
			coll.integrate_p_pp1_sub(fVals, m, p, ih);
			double ta = p == 0 ? coll.coll_sub[m].tleft : coll.coll_sub[m].nodes[p-1];
			double iex=poly.integrate(ta, coll.coll_sub[m].nodes[p]);
			if(abs(ih-iex) > 1e-7) {
				std::cerr << "error in node-to-node:\n";
				cerr << "m, p, node[p-1], node[p]: " << m << " " << p << " " << ta << " " << coll.coll_sub[m].nodes[p] << endl;
				cerr << "ih iex " << ih << " " << iex << endl;
				cerr << "polynom: [";	for(auto& d : poly.coeffs) cerr << " " << d ;	cerr << "]" << endl;
				cerr << "diff:" << abs(ih-iex) << endl;
				throw std::runtime_error("wrong node-to-node p_pp1 integration");
			}
		}
	}
}

void testPoly()
{
	array<double, 2 > nodes={1,2};
	Helper::LagrangePoly<2> lPoly0(nodes, 0);
	Helper::MonomPoly monPoly0(toMonom(lPoly0));
	assert(monPoly0.coeffs[0]==-1);
	assert(monPoly0.coeffs[1]==2);
	Helper::LagrangePoly<2> lPoly1(nodes, 1);
	Helper::MonomPoly monPoly1(toMonom(lPoly1));
	assert(monPoly1.coeffs[0]==1);
	assert(monPoly1.coeffs[1]==-1);

	assert(monPoly0.integrate(0.0, 1.0)==1.5);
}

template<unsigned deg >
PolynomMonomBase<deg> getPoly()
{
	array<double, deg+1> coeffs;
	for(unsigned i(0); i <= deg; ++i)
		coeffs[i] = rand() % 10+1.0;
	return PolynomMonomBase<deg>(coeffs);
}

template<unsigned M, unsigned P>
void testMRC(double t0, double t1)
{
	MultirateCollocation<M, P> coll;
	auto poly(getPoly<1>());
	test_m_mp1(poly, coll); //, 0.0, 1.0);
	test_m_mp1_sub(poly, coll); //, 0.0, 1.0);
	test_p_pp1(poly, coll);
	test_p_pp1_sub(poly, coll);
}

int main(int argc, char* argv[])
{
	testPoly();
	testMRC<2,2>(0,1);
	testMRC<3,2>(0,1);
	testMRC<2,2>(1,4);
	testMRC<3,2>(1,4);
	testMRC<2,3>(1,4);
	testMRC<4,5>(1,4);
	cout << "OK" << endl;
	return 0;
}
