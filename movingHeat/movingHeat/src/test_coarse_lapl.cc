#include <config.h>
#include <dune/grid/yaspgrid.hh>
#include <dune/common/function.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include "Heat.h"
#include "CoarseLapl.h"
#include "CoarseLapl.hh"

using namespace Dune;

int main(int argc, char* argv[])
{
  MPIHelper::instance(argc, argv);

  double nu(1.0);
  int nInter(5);
  unsigned  globRefine(1);
  static const int dim = 2;
  typedef YaspGrid<dim > GridType;
  typedef typename GridType::LeafGridView LeafView;
  typedef Functions::PQkNodalBasis<LeafView,1> LeafBasis;

  auto func = [](const auto& x) { return 0.1*cos(x[0]*2*M_PI)*cos(2*x[1]*0.5*M_PI); };
  auto laplFunc = [](const auto& x) { return -0.1*M_PI*M_PI*(4+1)*cos(x[0]*2*M_PI)*cos(2*x[1]*0.5*M_PI); };
  //auto func = [](const auto& x) { return 0.1*(x[0]*x[0]+x[1]*x[1]); };
  std::array<int, dim> resolution({nInter,4*nInter});
  GridType grid({1.0,4.0}, resolution); 
  grid.globalRefine(globRefine);

  Heat heat(nInter, globRefine, nu);

  CoarseLapl coarseLapl(grid, nu);

  LeafView leafView(grid.leafGridView());
  LeafBasis leafBasis(leafView);
  CoarseLapl::VectorType fineIn(leafBasis.size()), fineOut(leafBasis.size()), fineLapl(leafBasis.size()), swap(leafBasis.size());
  CoarseLapl::VectorType eLapl(leafBasis.size());
  Functions::interpolate(leafBasis, fineIn, func);
  Functions::interpolate(leafBasis, eLapl, laplFunc);

  fineOut = 0.0;
  coarseLapl.apply(fineIn, swap);

  auto fineLaplMat(coarseLapl.getLapl(grid, nu, grid.maxLevel()));
  fineLaplMat.mv(fineIn, swap);
  heat.MinvV(swap, fineLapl);
  VTKWriter<LeafView> vtkWriter(leafView);
  //vtkWriter.addVertexData(fineOut, "T");
  vtkWriter.addVertexData(fineLapl, "fineLapl");
  vtkWriter.addVertexData(fineIn, "fineIn");
  vtkWriter.write("test_coarse_lapl");

  eLapl -= fineLapl;
  std::cout << "max-norm-diff: " << eLapl.infinity_norm() << std::endl;

  /*auto heatLapl(heat.getLaplacian());
  Heat::VectorType heatVec, heatIn2;
  heat.init(heatVec);
  heat.init(heatIn2, func);
  heatLapl.mv(heatIn2, heatVec);
  heat.MinvV(heatVec, heatIn2);
  heat.writeResult("heatVec", heatIn2);*/
  return 0;
}
