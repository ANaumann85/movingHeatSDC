#ifndef SDC_H
#define SDC_H

#include "MultirateCollocation.hh"
#include <array>
#include <iostream>
#include <string>
#include <functional>

template<typename Vec, unsigned M_>
struct Sdc
{
	static const unsigned M = M_;
	//static const unsigned P = 2;
	typedef std::array<Vec, M > US;
	typedef std::function<void(Vec&)> Init;

	US us, I_m_mp1, fuSlow, fuFast;
	Vec fVal, rhs;

	MultirateCollocation<M, 2> coll;
	unsigned nIter;

	const Init& init;

	Sdc(const Init& init, unsigned nIter, std::string mQuad):
		coll(mQuad, "equi_noleft"),
		nIter(nIter),
		init(init)
	{
		init(fVal); init(rhs);
		for(unsigned m(0); m < M; ++m) {
			init(us[m]);
			init(I_m_mp1[m]);
			init(fuSlow[m]);
			init(fuFast[m]);
		}
		//fuFast[M-1][0]=0.0;
	}

	template<typename F >
	void evaluate_f(F& f, const US& u)
	{
		for(unsigned m(0); m < M; ++m) {
			//todo: set correct time t
			const double t = coll.coll.nodes[m];
			f.slow(t, u[m], fuSlow[m]);
			f.fast(t, u[m], fuFast);
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
#endif
	template<typename F >
	void update_I_m_mp1(F& f, const US& u)
	{
		evaluate_f(f, u);
		update_I_m_mp1();
	}

	void update_I_m_mp1()
	{
		for(unsigned m(0); m < M; ++m) {
			coll.integrate_m_mp1(fuSlow, m, I_m_mp1[m]); 
			coll.integrate_m_mp1(fuFast, m, fVal);
			axpy(1.0, fVal, I_m_mp1[m]);
		}
	}

	template<typename F >
	void predict(F& f, const Vec& u0, double t0, double te, bool setInter=true)
	{
		if(setInter)
			coll.setInterval(t0, te);
		//Vec u0_step; u0_step=u0;

		//Vec swap; swap = u0;
		for(unsigned m(0); m < M; ++m)
		{
			//compute u^*_{m+1}, w.r. (M-dtm*J(t_{m+1}))u^*_{m+1}=Mu^0_{m}
			f.updateMatrix(coll.coll.nodes[m], coll.coll.delta_m[m]);
			if(m==0) {
				f.Mv(u0, rhs);//rhs=Mu^0_{m}
				f.fast(coll.coll.tleft, u0, fVal); //fVal=f(t_m, u^*_m)
				axpy(coll.coll.delta_m[m], fVal, rhs);
			} else {
				f.Mv(us[m-1], rhs);//rhs=Mu^0_{m}
				f.fast(coll.coll.nodes[m-1], us[m-1], fuFast[m-1]); //fuFast=f(t_m, u^*_m)
				axpy(coll.coll.delta_m[m], fuFast[m-1], rhs);
			}
			f.solveMaJ(rhs, us[m]);//us[m]=u^*_{m+1}
			f.slow(coll.coll.nodes[m], us[m], fuSlow[m]);
		}
		f.fast(coll.coll.nodes[M-1], us[M-1], fuFast[M-1]);
		//print(us);
	}

	template<typename F >
	void sweep(F& f, Vec& u0, double t0, double te, bool setInter=true)
	{
		if(setInter)
			coll.setInterval(t0, te);
		update_I_m_mp1(); //use fuSlow and fuFast
		//print(I_m_mp1, "I_m_mp1");

#if 0
		print(I_m_mp1, "I_m_mp1");
#endif

		/*Vec u0_step; u0_step = u0;
		Vec fuStar, fuPm1Old;
		init(fuStar); init(fuPm1Old);*/
		for(unsigned m(0); m < M; ++m) {
			f.updateMatrix(coll.coll.nodes[m], coll.coll.delta_m[m]);
			if(m==0) {
				f.Mv(u0, rhs); //rhs=Mu^{k+1}_{m}
				axpy(-coll.coll.delta_m[m], fuSlow[m], rhs); //use fuSlow
				axpy(1.0, I_m_mp1[m], rhs);
				//std::cout << "rhs0: " << rhs[0] << std::endl;
			} else {
				f.Mv(us[m-1], rhs);
				axpy(-coll.coll.delta_m[m], fuFast[m-1], rhs); //use fuFast
				f.fast(coll.coll.nodes[m-1], us[m-1], fuFast[m-1]); //update fuFast
				axpy(coll.coll.delta_m[m], fuFast[m-1], rhs);
				axpy(-coll.coll.delta_m[m], fuSlow[m], rhs); //use fuSlow
				axpy(1.0, I_m_mp1[m], rhs);
			}

			f.solveMaJ(rhs, us[m]);
			f.slow(coll.coll.nodes[m], us[m], fuSlow[m]); //update fuSlow
		}
		f.fast(coll.coll.nodes[M-1], us[M-1], fuFast[M-1]);
		//print(us, "us_new");
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
};
#endif
