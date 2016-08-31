#ifndef MRSDC_H
#define MRSDC_H

#include "MultirateCollocation.hh"
#include <array>
#include <iostream>
#include <string>

template<typename Vec, unsigned M_, unsigned P_>
struct MRSdc
{
	static const unsigned M = M_;
	static const unsigned P = P_;
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
			coll.integrate_m_mp1_sub(fu_sub[m], m, iVal);
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

#if 1
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
#endif

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
	}

	template<typename F >
	void sweep(F& f, Vec& u0, double t0, double te)
	{
		update_I_m_mp1(f, us, ue);
		update_I_p_pp1(f, us, ue);

#if 0
		print(I_m_mp1, "I_m_mp1");
		print(I_p_pp1, "I_p_pp1");
#endif

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
#if 0
		print(us_new, "us_new");
		print(ue_new, "ue_new");
#endif
		us = us_new;
		ue = ue_new;

	}

	double residual(const Vec& u0)
	{
		double ret=norm(us[0]-u0-I_m_mp1[0]);
		for(unsigned m(1); m < M; ++m) {
			ret = std::max(ret, norm(us[m]-us[m-1]-I_m_mp1[m]));
		}
		return ret;
	}

	double sub_residual(const Vec& u0)
	{
		double ret(0.0);
		ret = sub_residual_m(u0, ue[0], 0);
		for(unsigned m(1); m < M; ++m) {
			ret = std::max(ret, sub_residual_m(ue[m-1][P-1], ue[m], m));
		}
		return ret;
	}

	double sub_residual_m(const Vec& u0, const std::array<Vec, P>& uem, unsigned m)
	{
		double ret(norm(uem[0] - u0 - I_p_pp1[m][0]));
		for(unsigned p(1); p < P; ++p) {
			ret = std::max(ret, norm(uem[p] - uem[p-1] - I_p_pp1[m][p]));
		}
		return ret;
	}
};
#endif
