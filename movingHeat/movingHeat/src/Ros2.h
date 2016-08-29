#ifndef ROS2_H
#define ROS2_H

template< typename Vec >
class Ros2
{
  Vec k1, k2, rhs;
  const double gamma = 1+sqrt(0.5);

  template<typename F>
  void onestep(F& f,  Vec& y0, double t0, double dt)
  {
    f.updateMatrix(t0, dt*gamma);
    f(t0, y0, rhs);
    f.solveMaJ(rhs, k1);
    //r=f(t0+dt, y0+dt*k1)-2*M*k2
    axpy(dt, k1, y0);
    f(t0+dt, y0, rhs);
    f.Mv(k1, k2);
    axpy(2.0, k2, rhs);
    f.solveMaJ(rhs, k2);
    //y0 += 0.5*dt*(k1+k2); //0.5*3*dt*k1=dt*k1+0.5*dt*k1
    axpy(0.5*dt, k1, y0);
    axpy(0.5*dt, k2, y0);
  }

  public:
  template<typename Init >
  Ros2(const Init& init)
  {
    init(k1);
    init(k2);
    init(rhs);
  }

  template<typename F>
  void solve(F& f, Vec& y0, double t0, double te, unsigned nStep)
  {
    double dt=(te-t0)/nStep;
    for(unsigned i(0); i < nStep; ++i, t0+=dt) {
      onestep(f, y0, t0, dt);
      f.writeResult(t0+dt);
    }
    
  }
};
#endif
