#ifndef HEAT_HH
#define HEAT_HH

#include <config.h>
#include <vector>
#include <string>

#include <dune/common/function.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/geometrygrid.hh>

// { include_matrix_vector_begin }
#include <dune/istl/matrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
// { include_matrix_vector_end }
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/umfpack.hh>
//for writing
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include "CoarseLapl.h"
using namespace Dune;
/** This class represents the heat equation 
 *    \dot{u} = nu\Delta u 
 *    with moving source on boundary
 *      n\cdot \partial u=&\alpha (v-u)
 * We discretize it with (linear) finite elements and obtain the ODE
 *  M \dot{U} = Lu + B(t)*u+b(t)
 */
#define WITH_SLOW_SRC
class Heat
{
  static const int dim = 2;
  struct MovingFunction : AnalyticalCoordFunction< double, dim, dim,  MovingFunction>
  {
    std::array<double, dim > dx;

    MovingFunction()
    {
      for(unsigned i(0); i < dim; ++i)
        dx[i] = 0.0;
    }

    void evaluate ( const DomainVector &x, RangeVector &y ) const
    {
      for(unsigned i(0); i < dim ; ++i)
        y[i] = x[i]+dx[i];
    }
  };
  public:
  typedef YaspGrid<dim > GridType;
  typedef GridType::LeafGridView GridView;
  typedef BlockVector<FieldVector<double,1> > VectorType;
  typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;
  typedef Functions::PQkNodalBasis<GridView,1> Basis;
  private:
  typedef YaspGrid<dim,EquidistantOffsetCoordinates<double,dim> > HostGridType_MV;
  typedef GeometryGrid<HostGridType_MV, MovingFunction > GridType_MV;
  typedef GridType_MV::LeafGridView GridView_MV;
  typedef Functions::PQkNodalBasis<GridView_MV,1> Basis_MV;

  const FieldVector< double, dim > L, lower_mv, upper_mv; 
  std::array<int, dim > s;
  std::shared_ptr< GridType> grid;
  mutable MovingFunction mf;
  HostGridType_MV hgtmv;
  std::shared_ptr< GridType_MV> grid_mv; //, grid_mvb; 
  GridView gridView; 
  Basis basis;

  MatrixType mass, lapl, lapl0, constRobM;
  MatrixType mMaJ;
  MatrixType laplTilde, laplExpl;
  VectorType constRobB;

  typedef UMFPack<MatrixType > MSolver;
  std::shared_ptr<MSolver> mSolver, mmaJSolver;

  typedef VTKSequenceWriter< GridView > PvdWriter;
  typedef VTKSequenceWriter< GridView_MV > PvdWriter_MV;
  std::shared_ptr<PvdWriter > pvdWriter;
  std::shared_ptr<PvdWriter_MV > pvdWriter_mv;
  std::shared_ptr<CoarseLapl > coarseLapl;

  double nu, alpha, v0, sourceVal;
  double bAlph, bVal;
  unsigned nInter;
  bool useLapl0, addConstRobin; /*, useLaplTilde;*/
  bool useLaplExpl;
    
  void fillMatrices();
  void fillMatricesZeroLapl();
  void setLaplZero();
  void buildMatrices(double h);
  void setConstRobin();

  public:
#if 0
  struct LaplTilde
  {
    int mode;
    double fac;

    LaplTilde():
      mode(0), fac(0.0)
    {}
  };
#endif

  Heat(int nInter, double nu=1.0e-3, double alpha=1.0e-4, double v0=5.0, double source=100, 
      bool useLapl0=false, bool addConstRobin=false, double laplExpl=0.0);

  Heat(int nInter, unsigned nRef, double nu=1.0e-3, double alpha=1.0e-4, double v0=5.0, double source=100, 
      bool useLapl0=false, bool addConstRobin=false, double laplExpl=0.0);

  //sets nu and alpha to the new values and updates matrices
  void setParam(double nu, double alpha);

  void setbAlph(double balph)
  { bAlph = balph; if(addConstRobin) setConstRobin(); }
  double getBVal() const 
  { return bVal; }

  //updates M-a*J(t) with given a
  void updateMatrix(double t, double a)
  { 
    mMaJ = mass ; 
    mMaJ.axpy( -a, lapl); 
    if(addConstRobin) 
      mMaJ.axpy(-a, constRobM); 
    //mmaJSolver = make_shared<MSolver>(mMaJ);
    /*if(useLaplTilde)
      mMaJ.axpy(a, laplTilde);*/
  }

  //solves x, such that (M-aJ)x=rhs
  template<typename V >
  void solveMaJ(V& rhs, V& x)
  {
#if 0
    //std::cout << "build ilu with cg\n";
    MatrixAdapter<MatrixType,VectorType,VectorType> linearOperator(mMaJ);
    SeqILU0<MatrixType,VectorType,VectorType> preconditioner(mMaJ,1.0);
    CGSolver<VectorType> cg(linearOperator, preconditioner, 1e-10, 500,0);
#else
    UMFPack<MatrixType > cg(mMaJ);
    //mmaJSolver->apply(x, rhs, statistics);
#endif
    InverseOperatorResult statistics;
    cg.apply(x, rhs, statistics);
  }

  //computes out = J*yIn+B(t)*yIn+b(t)
  void operator()(double t, const VectorType& yIn, VectorType& out) const
  {
    lapl.mv(yIn, out);
    fastAdd(t, yIn, out);
    if(addConstRobin) {
      out += constRobB;
      constRobM.umv(yIn, out);
    }
    if(useLaplExpl)
      laplExpl.umv(yIn, out);
  }

  //adds the fast term, i.e. out += B(t)*yIn+b(t)
  void fastGrid(double t, const VectorType& yIn, VectorType& out) const;
  void fastRect(double t, const VectorType& yIn, VectorType& out) const;
  void fastFull(double t, const VectorType& yIn, VectorType& out) const;
  void fastAdd(double t, const VectorType& yIn, VectorType& out) const  
  { 
#ifdef HEAT_USE_FAST_RECT
    fastRect(t, yIn, out);
#else
    fastGrid(t, yIn, out); 
#endif
    /*if(useLaplTilde)
      laplTilde.umv(yIn, out);*/
    if(useLaplExpl)
      laplExpl.umv(yIn, out);

  }
  //{ fastGridZeroLapl(t, yIn, out); }
  //{ fastFull(t, yIn, out); }

  template<typename F >
  void fastBoundary(const VectorType& yIn, const F& flux, VectorType& out) const;

  void fastGridZeroLapl(double t, const VectorType& yIn, VectorType& out) const;

  //sets the fast term, i.e. out = B(t)*yIn+b(t)
  void fast(double t, const VectorType& yIn, VectorType& out) const
  { 
    out = 0.0; fastAdd(t, yIn, out); 
    if(addConstRobin) {
      out += constRobB; 
      //constRobM.umv(yIn, out);
    }
  }

  //compute the slow term, i.e. out = L*yIn
  void slow(double t, const VectorType& yIn, VectorType& out) const
  { 
    lapl.mv(yIn, out); 
    if(addConstRobin) { 
      constRobM.umv(yIn, out); 
    } 
    /*if(useLaplTilde) {
      laplTilde.mmv(yIn, out);
    }*/
  }

  void slowSrc(double , VectorType& out) const
  { 
#ifdef WITH_SLOW_SRC
    if(addConstRobin)  
      out = constRobB; 
    else
      out = 0.0;
#else
    out = 0.0;
#endif
  }

  //computes out = M*in
  void Mv(const VectorType& in, VectorType& out) const
  { mass.mv(in, out); }

  //computes out = M^{-1}in
  void MinvV(VectorType& in, VectorType& out) const
  { 
    InverseOperatorResult statistics;
    mSolver->apply(out, in, statistics);
  }

  void writeResult(std::string fname, const VectorType& sol);
  void startFile(std::string fname, VectorType& data)
  {
    shared_ptr< VTKWriter<GridView>> vtkWriter(new VTKWriter<GridView>(gridView));
    vtkWriter->addVertexData(data, "T");
    pvdWriter=std::shared_ptr<PvdWriter >(new PvdWriter(vtkWriter, fname));

    shared_ptr<VTKWriter<GridView_MV>> vtkWriter_mv(new VTKWriter<GridView_MV>(grid_mv->leafGridView()));
    pvdWriter_mv=std::shared_ptr<PvdWriter_MV>(new PvdWriter_MV(vtkWriter_mv, "grid_mv"));
  }

  void writeResult(double t) const
  { pvdWriter->write(t); pvdWriter_mv->write(t); }

  //initialize dest with func(x,y)
  template<typename F >
  inline void init(VectorType& dest, const F& func) const
  { Functions::interpolate(basis, dest, func); }

  inline void init(VectorType& dest) const
  { dest.resize(basis.size()); }

  const MatrixType& getLaplacian() const
  { return lapl; }
  const MatrixType& getMass() const
  { return mass; }
  MatrixType getLaplaceWithMove(double t);
  void addFastMatrix(double t, MatrixType& dest) const;
};
#endif
