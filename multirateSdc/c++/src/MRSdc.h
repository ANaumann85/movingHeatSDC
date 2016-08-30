#ifndef MRSDC_H
#define MRSDC_H

#include <array>
#include <iostream>
#include <string>

template<unsigned M >
struct Collocation
{
	typedef std::array<std::array<double, M+1>, M+1> Mat;
	std::array<double, M> delta_m, nodes;
	Mat sMat;
	double tleft;
};

template< unsigned M, unsigned P >
struct MultirateCollocation
{
	typedef std::array<std::array<double, M>, M> Mat;

	Mat Shat_mp;
	std::array<Mat, M> S_mnp;

	Collocation<M> coll;
	std::array<Collocation<P>, M> coll_sub;

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

		coll_sub[0].sMat[0]={ 0.  , 0.  ,  0.        };
		coll_sub[0].sMat[1]={ 0.  , 0.25, -0.08333333};
		coll_sub[0].sMat[2]={ 0.  , 0.08333333 , 0.08333333};

		coll_sub[1].sMat[0]={ 0. , 0.    ,      0.        };
		coll_sub[1].sMat[1]={ 0. , 0.5   ,     -0.16666667};
		coll_sub[1].sMat[2]={ 0. , 0.16666667 , 0.16666667};

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

};

template<typename Vec, unsigned M, unsigned P>
struct MRSdc
{
	typedef std::array<Vec, M > US;
	typedef	std::array<std::array<Vec, P>, M > UE;

	US us, I_m_mp1;
	Vec fVal, rhs;

	UE ue, I_p_pp1;

	MultirateCollocation<M, P> coll;

	template<typename F >
	void evaluate_f(F& f, const US& u, US& fu)
	{
		for(unsigned m(0); m < M; ++m) {
			//todo: set correct time t
			const double t = 0.0;
			f.slow(t, u[m], fu[m]);
		}
	}
	
	template<typename F >
	void evaluate_f(F& f, const UE& ue, UE& fu_sub)
	{
		for(unsigned m(0); m < M; ++m) {
			for(unsigned p(0); p < P; ++p) {
				const double t = coll.coll_sub[m].nodes[p];
				f.fast(t, ue[m][p], fu_sub[m][p]);
			}
		}
	}

	template<typename F >
	void update_I_m_mp1(F& f, const US& u, const UE& usub)
	{
		US fu;
		UE fu_sub;
		Vec	iVal;
		evaluate_f(f, u, fu);
		evaluate_f(f, usub, fu_sub);

		for(unsigned m(0); m < M; ++m) {
			coll.integrate_m_mp1(fu, m, I_m_mp1[m]);
			print(I_m_mp1, "I_m_mp1[m]");
			coll.integrate_m_mp1_sub(fu_sub[m], m, iVal);
			std::cout << "I_m_mp1_sub[" << m << "]:" << iVal[0] << std::endl;
			axpy(1.0, iVal, I_m_mp1[m]);
		}
	}

	template<typename F >
	void update_I_p_pp1(F& f, const US& u, const UE& usub)
	{
		US fu;
		UE fu_sub;
		Vec	iVal;
		evaluate_f(f, u, fu);
		evaluate_f(f, usub, fu_sub);
		for(unsigned m(0); m < M; ++m) 
			for(unsigned p(0); p < P; ++p) {
				coll.integrate_p_pp1( fu, m, p, I_p_pp1[m][p]);
				coll.integrate_p_pp1_sub( fu_sub[m] , m, p, iVal);
				axpy(1.0, iVal, I_p_pp1[m][p]);
			}
	}

	void print(const US& us, std::string pref="us") const
	{
		std::cout << pref << ":[";
		for(auto& d: us)
			for(auto& d2:d)
				std::cout << " " << d2;
		std::cout << "]\n";

	}
	void print(const UE& ue, std::string pref="ue") const
	{
		std::cout << pref << ":[";
		for(auto& d: ue)
			for(auto& d2:d)
				for(auto& d3:d2)
					std::cout << " " << d3;
		std::cout << "]\n";

	}
	template<typename F >
	void predict(F& f, const Vec& u0, double t0, double te)
	{
		const double dt = te-t0;
		Vec u0_step; u0_step=u0;
		for(unsigned m(0); m < M; ++m)
		{
			f.updateMatrix(t0, coll.coll.delta_m[m]);
			f.solveMaJ(u0_step, us[m]);
		        ue[m][0] = u0_step ;

			double t = coll.coll_sub[m].tleft; 
			f.slow(t, us[m], fVal); 
		        axpy(coll.coll_sub[m].delta_m[0], fVal, ue[m][0]);
			f.fast(t, u0_step, fVal);
		        axpy(coll.coll_sub[m].delta_m[0], fVal, ue[m][0]);

			for(unsigned p(1); p < P; ++p) {
				t = coll.coll_sub[m].nodes[p-1];
				ue[m][p] = ue[m][p-1];
				f.slow(t, us[m], fVal); //same as before?
				axpy(coll.coll_sub[m].delta_m[p], fVal, ue[m][p]);
				f.fast(t, ue[m][p-1], fVal);
				axpy(coll.coll_sub[m].delta_m[p], fVal, ue[m][p]);
			}
			u0_step = ue[m][P-1];
			us[m]  = ue[m][P-1];
		}
		print(us, "us_predict");
		print(ue, "ue_predict");
		//std::cout << "ue:" << ue <<std::endl;
	}

	template<typename F >
	void sweep(F& f, Vec& u0, double t0, double te)
	{
		update_I_m_mp1(f, us, ue);
		update_I_p_pp1(f, us, ue);

		print(I_m_mp1, "I_m_mp1");
		print(I_p_pp1, "I_p_pp1");

		//TODO: should not need M new values
		std::array<Vec, M> us_new;
		//TODO: should not need M*P new values
		UE ue_new;
		Vec u0_step; u0_step = u0;
		for(unsigned m(0); m < M; ++m) {
			rhs = u0_step;
			f.slow(t0, us[m], fVal);
			axpy(- coll.coll.delta_m[m], fVal, rhs);
			axpy(1.0, I_m_mp1[m], rhs);
			f.updateMatrix(t0, coll.coll.delta_m[m]);
			f.solveMaJ(rhs, us_new[m]);

			double t = coll.coll_sub[m].tleft;
			ue_new[m][0] = u0_step;

			f.slow(t, us_new[m], fVal);
			axpy(coll.coll_sub[m].delta_m[0], fVal, ue_new[m][0]);
			f.slow(t, us[m], fVal);
			axpy(-coll.coll_sub[m].delta_m[0], fVal, ue_new[m][0]);

			//TODO: check!
			f.fast(t, u0_step, fVal);
			axpy(coll.coll_sub[m].delta_m[0], fVal, ue_new[m][0]);
			f.fast(t, u0_step, fVal);
			axpy(-coll.coll_sub[m].delta_m[0], fVal, ue_new[m][0]);
			axpy(1.0, I_p_pp1[m][0], ue_new[m][0]);
			for(unsigned p(1); p < P; ++p) {
				t = coll.coll_sub[m].nodes[p-1];
				ue_new[m][p]=ue_new[m][p-1];

				f.slow(t,us_new[m], fVal);
				axpy(coll.coll_sub[m].delta_m[p], fVal, ue_new[m][p]);
				f.slow(t,us[m], fVal);
				axpy(-coll.coll_sub[m].delta_m[p], fVal, ue_new[m][p]);

				f.fast(t, ue_new[m][p-1], fVal);
				axpy(coll.coll_sub[m].delta_m[p], fVal, ue_new[m][p]);
				f.fast(t, ue[m][p-1], fVal);
				axpy(-coll.coll_sub[m].delta_m[p], fVal, ue_new[m][p]);

				axpy(1.0, I_p_pp1[m][p], ue_new[m][p]);
			}
			u0_step = ue_new[m][P-1];
			us_new[m] = ue_new[m][P-1];
		}
		print(us_new, "us_new");
		print(ue_new, "ue_new");
		us = us_new;
		ue = ue_new;

	}
};
#endif
