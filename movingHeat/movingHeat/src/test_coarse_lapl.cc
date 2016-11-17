#include <config.h>
#include <dune/grid/yaspgrid.hh>
#include <dune/common/function.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include "CoarseLapl.h"
#include "CoarseLapl.hh"

using namespace Dune;

int main(int argc, char* argv[])
{
  MPIHelper::instance(argc, argv);

  static const int dim = 2;
  typedef YaspGrid<dim > GridType;
  typedef typename GridType::LeafGridView LeafView;
  typedef Functions::PQkNodalBasis<LeafView,1> LeafBasis;

  auto func = [](const auto& x) { return 0.1*(x[0]*x[0]+x[1]*x[1]); };
  int nInter(40);
  std::array<int, dim> resolution({nInter,4*nInter});
  GridType grid({1.0,4.0}, resolution); 
  grid.globalRefine(2);

  double nu(1.0);
  CoarseLapl coarseLapl(grid, nu);

  LeafView leafView(grid.leafGridView());
  LeafBasis leafBasis(leafView);
  CoarseLapl::VectorType fineIn(leafBasis.size()), fineOut(leafBasis.size());
  Functions::interpolate(leafBasis, fineIn, func);

  fineOut = 0.0;
  coarseLapl.apply(fineIn, fineOut);
  VTKWriter<LeafView> vtkWriter(leafView);
  vtkWriter.addVertexData(fineOut, "T");
  vtkWriter.write("test_coarse_lapl");
  return 0;
}
