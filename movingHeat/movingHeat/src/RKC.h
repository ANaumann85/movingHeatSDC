#ifndef RKC_H
#define RKC_H

#include <vector>
#include <stdexcept>

using namespace std;

template< typename VT >
class RKC
{
  static constexpr double eps = 0.15;

  double om0, om1;
  vector< double > b, c, mu, mub, nu, gamb;

  VT wn1, wn2, swap;
  VT curF, curF0;
  VT *cwn1, *cwn2;

  inline void shift()
  {
    VT* swap = cwn1;
    cwn1=cwn2;
    cwn2=swap;
  }

  template< typename F >
  void onestep(VT& wn, F& f, double tn, double h)
  {
    const unsigned s(b.size()-1);
    wn2 = wn;
    f(tn, wn, swap);
    f.MinvV(swap, curF0);
    wn1 = wn ;
    axpy((h*mub[1]), curF0, wn1);
    f(tn+h*c[0], wn1, swap);
    f.MinvV(swap, curF);
    cwn1=&wn1; cwn2=&wn2;
    for(unsigned j(2); j < s; ++j) {
      //(*cwn2) = (1-mu[j]-nu[j])*wn+mu[j]*(*cwn1)+nu[j]*(*cwn2)+(mub[j]*h)*curF+(gamb[j]*h)*curF0;
      (*cwn2) *= nu[j];
      axpy((1-mu[j]-nu[j]), wn, *cwn2);
      axpy(mu[j], *cwn1, *cwn2);
      axpy(mub[j]*h, curF, *cwn2);
      axpy(gamb[j]*h, curF0, *cwn2);
      f(tn+h*c[j-1], (*cwn2), swap);
      f.MinvV(swap, curF);
      shift();
    }
    //wn=(1-mu[s]-nu[s])*wn+mu[s]*(*cwn1)+nu[s]*(*cwn2)+(mub[s]*h)*curF+(gamb[s]*h)*curF0;
    wn *= (1-mu[s]-nu[s]);
    axpy(mu[s], *cwn1, wn);
    axpy(nu[s], *cwn2, wn);
    axpy(mub[s]*h, curF, wn);
    axpy(gamb[s]*h, curF0, wn);
  }

  inline double evalTRec(double x, unsigned j)
  { return j == 0 ? 1.0 : (j==1 ? x : 2*x*evalTRec(x,j-1)-evalTRec(x,j-2)); }

  inline double evalT(double x, unsigned j)
  { return abs(x) <= 1 ? cos(j*acos(x)) : (x > 1 ?  cosh(j*acosh(x)) : evalTRec(x, j)); }
  //{ return evalTRec(x, j); }

  inline double evalTP(double x, unsigned j)
    //{ return j==0 ? 0.0 : (j== 1 ? 1.0 : 2*evalT(x, j-1)+2*x*evalTP(x,j-1)-evalTP(x,j-2)); }
  { 
    double ret(j==0 ? 0.0 : (j== 1 ? 1.0 : 0));
    if(j > 1) {
      if(abs(x) <= 1)
        ret = j*sin(j*acos(x))/sin(acos(x));
      else if( x > 1) 
        ret = j*sinh(j*acosh(x))/sinh(acosh(x));
      else 
        throw std::runtime_error("no....");
    }
    return  ret; 
  }

  inline double evalTPP(double x, unsigned j)
  {
    double ret(0.0);
    switch(j)
    {
      case 0:
      case 1:
        break;
      default:
        //ret = 4*evalTP(x, j-1)+2*x*evalTPP(x,j-1)-evalTPP(x,j-2);
        if(abs(x) <= 1) {
          double s(1-x);
          ret = -j*j*(sqrt(1-x*x)-acos(x)*x)/sqrt(s*s);
        } else if( x > 1) {
          double s(x*x-1);s=sqrt(s);
          ret = j*(j*cosh(j*acosh(x))*s-sinh(j*acosh(x))*x)/(s*s*s);
          //ret = 4*evalTP(x, j-1)+2*x*evalTPP(x,j-1)-evalTPP(x,j-2);
        }
        else 
          throw std::runtime_error("no....");
    }
    return ret;
  }

  void initParams(unsigned s, unsigned order)
  {
    om0 = 1+eps/(s*s);
    if(order == 1) {
      om1 = evalT(om0,s)/evalTP(om0,s);
      for(unsigned i(0); i <= s; ++i)
        b[i] = 1.0/evalT(om0, i);
    } else {
      om1 = evalTP(om0,s)/evalTPP(om0,s);
      for(unsigned i(2); i <= s; ++i) {
        double tP(evalTP(om0, i));
        b[i] = evalTPP(om0, i)/(tP*tP);
      }
      b[0] = b[2];  b[1] = b[2];
    } 
    vector<double> a(s+1);
    a[0]=0;
    for(unsigned i(1); i <=s; ++i)
      a[i] =1.0-b[i]*evalT(om0,i);
    mub[1] = b[1]*om1;
    for(unsigned i(2); i <= s; ++i) {
      mu[i] = 2*b[i]*om0/(b[i-1]);
      nu[i] = -b[i]/b[i-2];
      mub[i] = 2*b[i]*om1/b[i-1];
      gamb[i] = -a[i-1]*mub[i];
    }
    if(order == 1) {
      for(unsigned i(0); i < s-1; ++i)
        c[i] = om1*evalTP(om0, i+1)/evalT(om0, i+1);
    }else {
      for(unsigned i(1); i < s-1; ++i)
        c[i] = om1*evalTPP(om0, i+1)/evalTP(om0, i+1);
      c[0] = c[1]/evalTP(om0,2);
    }
#if 0
    cout << "b: ";
    for(auto& d : b)
      cout << d << " " ;
    cout << endl;
    cout << "c: ";
    for(auto& d : c)
      cout << d << " " ;
    cout << endl;
    cout << "mu: ";
    for(auto& d : mu)
      cout << d << " " ;
    cout << endl;
    cout << "mub: ";
    for(auto& d : mub)
      cout << d << " " ;
    cout << endl;
    cout << "nu: ";
    for(auto& d : nu)
      cout << d << " " ;
    cout << endl;
    cout << "gamb: ";
    for(auto& d : gamb)
      cout << d << " " ;
    cout << endl;
#endif
  }
  public:
  template<typename Init >
  RKC(const Init& init, unsigned s, unsigned order):
    b(s+1), c(s-1), mu(s+1),
    mub(s+1), nu(s+1), gamb(s+1)
  {
    if(order != 1 && order != 2)
      throw std::runtime_error("only first and second order");
    initParams(s, order);
    init(wn1);
    init(wn2);
    init(curF);
    init(curF0);
    init(swap);
  }

  template< typename F, typename P>
  void solve(const F& f, VT& y0, double t0, double te, unsigned nStep,  P& p)
  {
    double h((te-t0)/nStep);
    for(unsigned n(0); n < nStep; ++n) {
      onestep(y0, f, t0, h);
      t0 += h;
      p(t0, y0);
    }
  }
};
#endif

