#ifndef MULTIRATE_COLLOCATION_HH
#define MULTIRATE_COLLOCATION_HH

#include <array>
#include <fstream>
#include <sstream>
#include <string>
#include <cassert>
#include <iostream>

template<unsigned M >
struct Collocation
{
	typedef std::array<std::array<double, M>, M> Mat;
	typedef std::array<std::array<double, M>, M+2> Mat2;
	typedef std::array<double, M> Vec;

	Vec delta_m, nodes;
	Mat sMat;
	double tleft;

	void set(const Mat2& data, double t0, double t1)
	{
		double dt=t1-t0;
		assert(dt > 0);
		for(unsigned mr(0); mr < M; ++mr) {
			for(unsigned mc(0); mc < M; ++mc) {
				sMat[mr][mc] = dt*data[mr][mc];
			}
			nodes[mr] = t0+dt*data[M][mr];
			delta_m[mr] = mr == 0 ? nodes[mr] : nodes[mr]-nodes[mr-1];
		}
		tleft = t0;
	}

	template<typename F >
	Vec evalAtNodes(const F& f) const
	{
		Vec ret;
		for(unsigned i(0); i < M; ++i)
			ret[i] = f(nodes[i]);
		return ret;
	}
};

template< unsigned M, unsigned P >
struct MultirateCollocation
{
	typedef std::array<std::array<double, M>, M> Mat;
	typedef std::array<std::array<double, M>, M+2> MatMp2;
	typedef std::array<std::array<double, P>, P+2> MatPp2;

	template<std::size_t S>
	void readMatrix(std::array<std::array<double, S>, S+2>& dest, std::string fname)
	{
		std::fstream file(fname);
		std::string line;
		unsigned r(0);
		getline(file, line);
		while(file.good() || line.length() > 0) {
			std::stringstream sl;
			sl << line;
			double col;
			unsigned j(0); 
			sl >> col;
			while(sl.good()) {
				dest[r][j] = col;
				++j;
				sl >> col;
			}
			assert(j == S);
			++r;
			line.clear();
			getline(file, line);
		}
		assert(r == S+2);
		file.close();
	}

#if 0
	template<std::size_t S >
	void print(std::array<std::array<double, S>, S+2>& mat)
	{
		std::cout << "[";
		for(auto& d:mat) {
			std::cout << "["; for(auto& d2:d) std::cout << " " << d2;
			std::cout << "]\n";
		}
		std::cout << "]\n";
	}
#endif


	MultirateCollocation()
	{
		MatMp2 sMat_M;
		MatPp2 sMat_P;
		readMatrix(sMat_M, "sdc_quad_weights/radau_right-M2.dat");
		readMatrix(sMat_P, "sdc_quad_weights/equi_noleft-M2.dat");
#if 0
		print(sMat_M);
		print(sMat_P);
#endif
		coll.set(sMat_M, 0.0, 1.0);

		Shat_mp[0] = { 0.33333333 , 0.};
		Shat_mp[1] = { 0.66666667,  0.};

		S_mnp[0][0]={ 0.22916667,  0.1875    };
		S_mnp[0][1]={-0.0625  ,   -0.02083333};
		S_mnp[1][0]={ 0.25   ,     0.08333333};
		S_mnp[1][1]={ 0.08333333 , 0.25      };
	

		coll_sub[0].set(sMat_P, coll.tleft, coll.nodes[0]);
		for(unsigned m(1); m < M; ++m) {
			coll_sub[m].set(sMat_P, coll.nodes[m-1], coll.nodes[m]);
		}
			
		/*coll_sub[0].sMat[0]={  0.25, -0.08333333};
		coll_sub[0].sMat[1]={  0.08333333 , 0.08333333};
		coll_sub[0].nodes={1.0/6.0, 1.0/3.0};
		coll_sub[0].tleft=0.0;

		coll_sub[1].sMat[0]={  0.5   ,     -0.16666667};
		coll_sub[1].sMat[1]={  0.16666667 , 0.16666667};
		coll_sub[1].nodes={2.0/3.0, 1.0};
		coll_sub[1].tleft=1.0/3.0;*/

		/*coll.delta_m={0.33333333 , 0.66666667};
		coll_sub[0].delta_m={0.16666667 , 0.16666667};
		coll_sub[1].delta_m={0.33333333 , 0.33333333};*/
	}
	

	template<typename Vec >
	void integrate_m_mp1(const std::array<Vec, M>& fu, unsigned m, Vec& res)
	{
		setValue(res, 0.0);
		for(unsigned j(0); j < M; ++j) {
			axpy(coll.sMat[m][j], fu[j], res);
		}
	}

	template<typename Vec >
	void integrate_m_mp1_sub(const std::array<Vec, P>& fu_sub, unsigned m, Vec& iVal)
	{
		setValue(iVal, 0.0);
		for(unsigned p(0); p < P; ++p)
			axpy(Shat_mp[m][p], fu_sub[p], iVal);
	}

	template<typename Vec >
	void integrate_p_pp1(std::array<Vec, M>& fu, unsigned m, unsigned p, Vec& iVal)
	{
		setValue(iVal, 0.0);
		for(unsigned n(0); n < M; ++n) {
			axpy(S_mnp[m][n][p] ,fu[n], iVal);
		}

	}

	template<typename Vec >
	void integrate_p_pp1_sub(const std::array<Vec, P>& fu_sub , unsigned m, unsigned p, Vec& iVal)
	{
		setValue(iVal, 0.0);
		const typename Collocation<M>::Mat& sMat(coll_sub[m].sMat);
		for(unsigned j(0); j < P; ++j) {
			axpy( sMat[p][j], fu_sub[j], iVal);
		}
	}

	Collocation<M> coll;
	std::array<Collocation<P>, M> coll_sub;
	private:
	Mat Shat_mp;
	std::array<Mat, M> S_mnp;


};
#endif
