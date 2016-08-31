#ifndef MULTIRATE_COLLOCATION_HH
#define MULTIRATE_COLLOCATION_HH

#include <array>
template<unsigned M >
struct Collocation
{
	typedef std::array<std::array<double, M+1>, M+1> Mat;
	typedef std::array<double, M> Vec;
	Vec delta_m, nodes;
	Mat sMat;
	double tleft;

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

	MultirateCollocation()
	{
		Shat_mp[0] = { 0.33333333 , 0.};
		Shat_mp[1] = { 0.66666667,  0.};

		S_mnp[0][0]={ 0.22916667,  0.1875    };
		S_mnp[0][1]={-0.0625  ,   -0.02083333};
		S_mnp[1][0]={ 0.25   ,     0.08333333};
		S_mnp[1][1]={ 0.08333333 , 0.25      };
	
		coll.sMat[0]={ 0. , 0. ,         0.        }; 
		coll.sMat[1]={ 0. , 0.41666667, -0.08333333};
		coll.sMat[2]={ 0. , 0.33333333,  0.33333333};
		coll.nodes = { 1.0/3.0, 1.0};
		coll.delta_m = {1.0/3.0, 2.0/3.0};
		coll.tleft = 0.0;

		coll_sub[0].sMat[0]={ 0.  , 0.  ,  0.        };
		coll_sub[0].sMat[1]={ 0.  , 0.25, -0.08333333};
		coll_sub[0].sMat[2]={ 0.  , 0.08333333 , 0.08333333};
		coll_sub[0].nodes={1.0/6.0, 1.0/3.0};
		coll_sub[0].tleft=0.0;

		coll_sub[1].sMat[0]={ 0. , 0.    ,      0.        };
		coll_sub[1].sMat[1]={ 0. , 0.5   ,     -0.16666667};
		coll_sub[1].sMat[2]={ 0. , 0.16666667 , 0.16666667};
		coll_sub[1].nodes={2.0/3.0, 1.0};
		coll_sub[1].tleft=1.0/3.0;

		coll.delta_m={0.33333333 , 0.66666667};
		coll_sub[0].delta_m={0.16666667 , 0.16666667};
		coll_sub[1].delta_m={0.33333333 , 0.33333333};
	}
	

	template<typename Vec >
	void integrate_m_mp1(const std::array<Vec, M>& fu, unsigned m, Vec& res)
	{
		setValue(res, 0.0);
		for(unsigned j(0); j < M; ++j) {
			axpy(coll.sMat[m+1][j+1], fu[j], res);
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
			axpy( sMat[p+1][j+1], fu_sub[j], iVal);
		}
	}

	Collocation<M> coll;
	std::array<Collocation<P>, M> coll_sub;
	private:
	Mat Shat_mp;
	std::array<Mat, M> S_mnp;


};
#endif