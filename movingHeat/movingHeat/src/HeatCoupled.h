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

#include <array>
using namespace Dune;
/** This class represents the heat equation 
 *    \dot{u} = nu\Delta u 
 *    with moving source on boundary
 *      n\cdot \partial u=&\alpha (v-u)
 * We discretize it with (linear) finite elements and obtain the ODE
 *  M \dot{U} = Lu + B(t)*u+b(t)
 */
#define WITH_SLOW_SRC
class HeatCoupled
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
  typedef std::array<BlockVector<FieldVector<double,1> >, 2> VectorType;
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
  GridView_MV gridView_mv;
  Basis basis;
  Basis_MV basis_mv;

  std::array<MatrixType, 2> mass, lapl, lapl0;
  MatrixType constRobM;
  std::array<MatrixType, 2> mMaJ;
  VectorType constRobB;

  typedef UMFPack<MatrixType > MSolver;
  std::array<std::shared_ptr<MSolver>, 2> mSolver;

  typedef VTKSequenceWriter< GridView > PvdWriter;
  typedef VTKSequenceWriter< GridView_MV > PvdWriter_MV;
  std::shared_ptr<PvdWriter > pvdWriter;
  std::shared_ptr<PvdWriter_MV > pvdWriter_mv;

  double nu, alpha, v0, sourceVal;
  double bAlph, bVal;
  unsigned nInter;
  bool useLapl0, addConstRobin;
    
  void fillMatrices();
  void fillMatricesZeroLapl();
  void setLaplZero();
  void buildMatrices();
  void setConstRobin();

  public:
  HeatCoupled(int nInter, double nu=1.0e-3, double alpha=1.0e-4, double v0=5.0, double source=100, 
      bool useLapl0=false, bool addConstRobin=false);

  //sets nu and alpha to the new values and updates matrices
  void setParam(double nu, double alpha);

  void setbAlph(double balph)
  { bAlph = balph; if(addConstRobin) setConstRobin(); }
  double getBVal() const 
  { return bVal; }

  //updates M-a*J(t) with given a
  void updateMatrix(double t, double a)
  { 
    for(unsigned i(0); i < 2; ++i) {
      mMaJ[i] = mass[i] ; mMaJ[i].axpy( -a, lapl[i]);
    }
    if(addConstRobin) 
      mMaJ[0].axpy(-a, constRobM); 
  }

  //solves x, such that (M-aJ)x=rhs
  //TODO: move construction to updateMatrix and reuse
  template<typename V >
  void solveMaJ(V& rhs, V& x)
  {
    for(unsigned i(0); i <2 ;++i) {
    UMFPack<MatrixType > cg(mMaJ[i]);
    InverseOperatorResult statistics;
    cg.apply(x[i], rhs[i], statistics);
    }
  }

  //computes out = J*yIn+B(t)*yIn+b(t)
  void operator()(double t, const VectorType& yIn, VectorType& out) const
  {
    for(unsigned i(0); i < 2; ++i) {
      lapl[i].mv(yIn[i], out[i]);
      fastAdd(t, yIn, out);
    }
    if(addConstRobin) {
      out[0] += constRobB[0];
      constRobM.umv(yIn[0], out[0]);
    }
  }

  //adds the fast term, i.e. out += B(t)*yIn+b(t)
  void fastGrid(double t, const VectorType& yIn, VectorType& out) const;
#if 0
  void fastRect(double t, const VectorType& yIn, VectorType& out) const;
  void fastFull(double t, const VectorType& yIn, VectorType& out) const;
#endif
  void fastAdd(double t, const VectorType& yIn, VectorType& out) const
  //{ fastRect(t, yIn, out); }
  { fastGrid(t, yIn, out); }
  //{ fastGridZeroLapl(t, yIn, out); }
  //{ fastFull(t, yIn, out); }

  template<typename F >
  void fastBoundary(const VectorType& yIn, const F& flux, VectorType& out) const;

  //void fastGridZeroLapl(double t, const VectorType& yIn, VectorType& out) const;

  //sets the fast term, i.e. out = B(t)*yIn+b(t)
  void fast(double t, const VectorType& yIn, VectorType& out) const
  { 
    out[0] = 0.0, out[1] = 0.0; fastAdd(t, yIn, out); 
    if(addConstRobin) {
      out[0] += constRobB[0]; 
      //constRobM.umv(yIn, out);
    }
  }

  //compute the slow term, i.e. out = L*yIn
  void slow(double t, const VectorType& yIn, VectorType& out) const
  { 
    for(unsigned i(0); i < 2; ++i)
      lapl[i].mv(yIn[i], out[i]); 
    if(addConstRobin) { 
      constRobM.umv(yIn[0], out[0]); 
    } 
  }

  void slowSrc(double , VectorType& out) const
  { 
#ifdef WITH_SLOW_SRC
    if(addConstRobin)  {
      out[0] = constRobB[0]; 
      out[1] = constRobB[1];
    }
    else {
      out[0] = 0.0;
      out[1] = 0.0;
    }
#else
    out = 0.0;
#endif
  }

  //computes out = M*in
  void Mv(const VectorType& in, VectorType& out) const
  { for(unsigned i(0); i < 2; ++i) mass[i].mv(in[i], out[i]); }

  //computes out = M^{-1}in
  void MinvV(VectorType& in, VectorType& out) const
  { 
    for(unsigned i(0); i < 2; ++i) {
      InverseOperatorResult statistics;
      mSolver[i]->apply(out[i], in[i], statistics);
    }
  }

  void writeResult(std::string fname, std::string fname_mv, const VectorType& sol);
  void startFile(std::string fname, std::string fname_mv, VectorType& data)
  {
    shared_ptr< VTKWriter<GridView>> vtkWriter(new VTKWriter<GridView>(gridView));
    vtkWriter->addVertexData(data[0], "T");
    pvdWriter=std::shared_ptr<PvdWriter >(new PvdWriter(vtkWriter, fname));

    shared_ptr<VTKWriter<GridView_MV>> vtkWriter_mv(new VTKWriter<GridView_MV>(grid_mv->leafGridView()));
    vtkWriter_mv->addVertexData(data[1], "T");
    pvdWriter_mv=std::shared_ptr<PvdWriter_MV>(new PvdWriter_MV(vtkWriter_mv, fname_mv));
  }

  void writeResult(double t) const
  { pvdWriter->write(t); pvdWriter_mv->write(t); }

  //initialize dest with func(x,y)
  template<typename F >
  inline void init(VectorType& dest, const F& func) const
  { Functions::interpolate(basis, dest, func); }

  inline void init(VectorType& dest) const
  { 
    dest[0].resize(basis.size()); 
    dest[1].resize(basis_mv.size());
  }

  const MatrixType& getLaplacian(unsigned pos) const
  { return lapl[pos]; }
  const MatrixType& getMass(unsigned pos) const
  { return mass[pos]; }
  MatrixType getLaplaceWithMove(double t);
  void addFastMatrix(double t, MatrixType& dest) const;
};
#endif
