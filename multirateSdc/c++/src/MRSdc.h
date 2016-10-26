#ifndef MRSDC_H
#define MRSDC_H

#include "MultirateCollocation.hh"
#include <array>
#include <iostream>
#include <string>
#include <functional>

template<typename Vec, unsigned M_, unsigned P_>
struct MRSdc
{
	static const unsigned M = M_;
	static const unsigned P = P_;
	typedef std::array<Vec, M > US;
	typedef	std::array<std::array<Vec, P>, M > UE;
	typedef std::function<void(Vec&)> Init;

	US us, I_m_mp1, fus;
	Vec fVal, rhs;

	UE ue, I_p_pp1, fue;

	MultirateCollocation<M, P> coll;
	unsigned nIter;

	const Init& init;

	MRSdc(const Init& init, unsigned nIter, std::string mQuad, std::string pQuad):
		coll(mQuad, pQuad),
		nIter(nIter),
		init(init)
	{
		init(fVal); init(rhs);
		for(unsigned m(0); m < M; ++m) {
			init(us[m]);
			init(I_m_mp1[m]);
			init(fus[m]);
		}
		for(unsigned m(0); m < M; ++m) 
			for(unsigned p(0); p < P; ++p) {
				init(ue[m][p]);
				init(I_p_pp1[m][p]);
				init(fue[m][p]);
			}
	}

	template<typename F >
	void evaluate_f(F& f, const US& u, US& fu)
	{
		for(unsigned m(0); m < M; ++m) {
			//todo: set correct time t
			const double t = coll.coll.nodes[m];
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
		evaluate_f(f, u, fus);
		evaluate_f(f, usub, fue);
		update_I_m_mp1();

	}

	void update_I_m_mp1()
	{
		for(unsigned m(0); m < M; ++m) {
			coll.integrate_m_mp1(fus, m, I_m_mp1[m]);
			coll.integrate_m_mp1_sub(fue[m], m, fVal);
			axpy(1.0, fVal, I_m_mp1[m]);
		}
	}

	template<typename F >
	void update_I_p_pp1(F& f, const US& u, const UE& usub)
	{
		evaluate_f(f, u, fus);
		evaluate_f(f, usub, fue);

		update_I_p_pp1();
	}

	void update_I_p_pp1()
	{
		for(unsigned m(0); m < M; ++m) 
			for(unsigned p(0); p < P; ++p) {
				coll.integrate_p_pp1( fus, m, p, I_p_pp1[m][p]);
				coll.integrate_p_pp1_sub( fue[m] , m, p, fVal);
				axpy(1.0, fVal, I_p_pp1[m][p]);
			}
	}
#if 0
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
	void predict(F& f, const Vec& u0, double t0, double te, bool setInter=true)
	{
		if(setInter)
			coll.setInterval(t0, te);
		Vec u0_step; u0_step=u0;

		//Vec swap; swap = u0;
		f.fast(coll.coll_sub[0].tleft, u0_step, fVal);
		for(unsigned m(0); m < M; ++m)
		{
			f.Mv(u0_step, rhs);//rhs=Mu^0_{m}
			//compute u^*_{m+1}, w.r. (M-dtm*J(t_{m+1}))u^*_{m+1}=Mu^0_{m}
			f.updateMatrix(coll.coll.nodes[m], coll.coll.delta_m[m]);
			f.solveMaJ(rhs, us[m]);//us[m]=u^*_{m+1}
			//compute Mu^0_{m,p-1}+dt_{m+1,p} f(u^*_{m+1})+dt_{m+1,p}g(u^0_{m,p-1})
			f.Mv(u0_step, rhs);//rhs=Mu^0_{m,0} 

			//f.fast(coll.coll_sub[m].tleft, u0_step, fVal);//fVal=g(t_m, u^0_m) //use fast val from line 135 and line 166
			axpy(coll.coll_sub[m].delta_m[0], fVal, rhs);//rhs=Mu^0_{m,0}+dt_{m,0}*g(t_m,u^0_m)
			f.slow(coll.coll.nodes[m], us[m], fVal); //fVal=f(t_m, u^*_m)
			axpy(coll.coll_sub[m].delta_m[0], fVal, rhs);//rhs=Mu^0_{m,0}+dt_{m,0}*g(t_m,u^0_m)+dt_{m,0}f(t_m,u^*_m)
			f.MinvV(rhs, ue[m][0]); //ue[m][0]=u^0_{m,0}+M^{-1}(dt_{m,0}f(t_m,u^*_m)+dt_{m,0}*g(t_m,u^0_m))
			f.fast(coll.coll_sub[m].nodes[0], ue[m][0], fue[m][0]);//store fue[m][0]=g(u^1_{m,0})

			for(unsigned p(1); p < P; ++p) {
				//compute Mu^0_{m,p-1}+dt_{m+1,p} f(u^*_{m+1})+dt_{m+1,p}g(u^0_{m,p-1})
				f.Mv(ue[m][p-1], rhs); //rhs=Mu^0_{m,p-1}
				f.slow(coll.coll.nodes[m], us[m], fVal); //fVal=f(t_m, u^*_m) //same as before!
				axpy(coll.coll_sub[m].delta_m[p], fVal, rhs);//rhs=Mu^0_{m,0}+dt_{m,p}f(t_m,u^*_m)
				//f.fast(coll.coll_sub[m].nodes[p-1], ue[m][p-1], fVal);//fVal=g(t_{m,p-1}, u^0_{m,p-1})
				//axpy(coll.coll_sub[m].delta_m[p], fVal, rhs); //rhs=Mu^0_{m,0}+dt_{m,p}f(t_m,u^*_m)+dt_{m,p}g(t_{m,p-1}, u^0_{m,p-1})

				axpy(coll.coll_sub[m].delta_m[p], fue[m][p-1], rhs); //rhs=Mu^0_{m,0}+dt_{m,p}f(t_m,u^*_m)+dt_{m,p}g(t_{m,p-1}, u^0_{m,p-1})
				f.MinvV(rhs, ue[m][p]); //ue[m][p]=u^0_{m,p-1}+M^{-1}(dt_{m,p}f(t_m,u^*_m)+dt_{m,p}*g(t_m,u^0_{m,p-1}))
				f.fast(coll.coll_sub[m].nodes[p], ue[m][p], fue[m][p]);//store fue[m][p]=g(u^1_{m,p})
			}
			u0_step = ue[m][P-1];
			us[m]  = ue[m][P-1]; //update u^1_{m}=u^1_{m,P}
			fVal = fue[m][P-1];
			f.slow(coll.coll.nodes[m], us[m], fus[m]); //store fus[m]=f(u^1_{m})
		}
	}

	template<typename F >
	void sweep(F& f, Vec& u0, double t0, double te, bool setInter=true)
	{
		if(setInter)
			coll.setInterval(t0, te);
		update_I_m_mp1();
		update_I_p_pp1();

#if 0
		print(I_m_mp1, "I_m_mp1");
		print(I_p_pp1, "I_p_pp1");
#endif

		Vec u0_step; u0_step = u0;
		Vec fuStar, fuPm1Old;
		init(fuStar); init(fuPm1Old);
		for(unsigned m(0); m < M; ++m) {
			//standard step for u*
			f.Mv(u0_step, rhs); //rhs=Mu^{k+1}_{m}
#if 0
			f.slow(coll.coll.nodes[m], us[m], fVal); //fVal=f(t_m, u^k_m)
			axpy(- coll.coll.delta_m[m], fVal, rhs); //rhs=Mu^{k+1}_{m}-dt_{m}f(t_m,u^k_m)
#else
			axpy(- coll.coll.delta_m[m], fus[m], rhs); //rhs=Mu^{k+1}_{m}-dt_{m}f(t_m,u^k_m)
#endif
			axpy(1.0, I_m_mp1[m], rhs); //rhs=Mu^{k+1}_{m}-dt_{m}f(t_m,u^k_m)+I^m_{m}
			f.updateMatrix(coll.coll.nodes[m], coll.coll.delta_m[m]); //set (M-dt_m J(t_m))
			f.solveMaJ(rhs, fuPm1Old); //u^*_m=fuPm1Old=(M-dt_m J(t_m))^{-1}(Mu^{k+1}_{m}-dt_{m}f(t_m,u^k_m)+I^m_{m})
			f.slow(coll.coll.nodes[m], fuPm1Old, fuStar); //fuStar=f(t_m, u^*_m)

			//embedded steps
			f.Mv(u0_step, rhs); //rhs=M u^{k+1}_{m}
			axpy(coll.coll_sub[m].delta_m[0], fuStar, rhs); //rhs=M u^{k+1}_{m}+dt_{m,0}f(t_m, u^*_m)
			axpy(-coll.coll_sub[m].delta_m[0], fus[m], rhs); //rhs = Mu^{k+1}_m+dt_{m,0}f(t_m, u^*_m)-dt_{m,0}f(t_m, u^k_m)


			//TODO: check!
			/*f.fast(t, u0_step, fVal);
			  axpy(coll.coll_sub[m].delta_m[0], fVal, ue_new[m][0]);
			  f.fast(t, u0_step, fVal);
			  axpy(-coll.coll_sub[m].delta_m[0], fVal, ue_new[m][0]);*/

			axpy(1.0, I_p_pp1[m][0], rhs); //rhs = Mu^{k+1}_m+dt_{m,0}f(t_m, u^*_m)-dt_{m,0}f(t_m, u^k_m)+I^p_{m,0}
			f.MinvV(rhs, ue[m][0]); //ue[m][0] = u^{k+1}_m+M^{-1}(dt_{m,0}f(t_m, u^*_m)-dt_{m,0}f(t_m, u^k_m)+I^{1}_{m,0})
			fuPm1Old = fue[m][0];
			f.fast(coll.coll_sub[m].nodes[0], ue[m][0], fue[m][0]); //store g(u^{k+1}_{m,p})

			for(unsigned p(1); p < P; ++p) {
				f.Mv(ue[m][p-1], rhs); //rhs=M u^{k+1}_{m,p-1}
				axpy(coll.coll_sub[m].delta_m[p], fuStar, rhs); //rhs=M u^{k+1}_{m,p-1}+dt_{m,p}f(t_m, u^*_m)
				axpy(-coll.coll_sub[m].delta_m[p], fus[m], rhs); //rhs=M u^{k+1}_{m,p-1}+dt_{m,p}f(t_m, u^*_m)-dt_{m,p}f(t_m,u^k_m)

				axpy(coll.coll_sub[m].delta_m[p], fue[m][p-1], rhs); //rhs=M u^{k+1}_{m,p-1}+dt_{m,p}f(t_m, u^*_m)-dt_{m,p}f(t_m,u^k_m)+dt_{m,p}g(t_{m,p-1}, u^{k+1}_{m,p-1})
				axpy(-coll.coll_sub[m].delta_m[p], fuPm1Old, rhs); //rhs=M u^{k+1}_{m,p-1}+dt_{m,p}f(t_m, u^*_m)-dt_{m,p}f(t_m,u^k_m)+dt_{m,p}g(t_{m,p-1}, u^{k+1}_{m,p-1})-dt_{m,p}g(t_{m,p-1}, u^k_{m,p-1})

				axpy(1.0, I_p_pp1[m][p], rhs);//rhs=M u^{k+1}_{m,p-1}+dt_{m,p}f(t_m, u^*_m)-dt_{m,p}f(t_m,u^k_m)+dt_{m,p}g(t_{m,p-1}, u^{k+1}_{m,p-1})-dt_{m,p}g(t_{m,p-1}, u^k_{m,p-1})+I^{p+1}_{m,p}
				f.MinvV(rhs, ue[m][p]); //ue[m][p]=u^{k+1}_{m,p-1}+M^{-1}(dt_{m,p}f(t_m, u^*_m)-dt_{m,p}f(t_m,u^k_m)+dt_{m,p}g(t_{m,p-1}, u^{k+1}_{m,p-1})-dt_{m,p}g(t_{m,p-1}, u^k_{m,p-1})+I^{p+1}_{m,p});
				fuPm1Old = fue[m][p]; //backup g(u^{k}_{m,p})
				f.fast(coll.coll_sub[m].nodes[p], ue[m][p], fue[m][p]); //store fue[m][p]=g(u^{k+1}_{m,p})
			}
			u0_step = ue[m][P-1]; //u0_step=u^{k+1}_{m,P}
			us[m] = ue[m][P-1]; //us[m]=u^{k+1}_{m,P}
			f.slow(coll.coll.nodes[m], us[m], fus[m]); //fus[m]=f(t_m, u^k_m)
		}
#if 0
		print(us_new, "us_new");
		print(ue_new, "ue_new");
#endif
		//us = us_new;
		//ue = ue_new;

	}

	template<typename F >
	void solve(F& f, Vec& u0, double t0, double te, unsigned nStep)
	{
		double dt=(te-t0)/((double) nStep);
		for(unsigned s(0); s < nStep; ++s) { //, t0+=dt
			double t0_ = t0+dt*s;
			predict(f, u0, t0_, t0_+dt);
			for(unsigned k(0); k < nIter; ++k) {
				sweep(f, u0, t0_, t0_+dt, false);
			}
			u0 = us[M-1];
		}
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
