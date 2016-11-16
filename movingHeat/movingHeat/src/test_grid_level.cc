#include <config.h>
#include <dune/grid/yaspgrid.hh>
#include <dune/common/function.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <array>
#include <string>
#include <sstream>

using namespace Dune;

template<typename GridView >
void write(const GridView& gridView, std::string fname)
{
  VTKWriter<GridView> vtkWriter(gridView);
  vtkWriter.write(fname);
}

int main(int argc, char* argv[])
{
  MPIHelper::instance(argc, argv);
  static const int dim = 2;
  typedef YaspGrid<dim > GridType;

  int nInter(5);
  std::array<int, dim> resolution({nInter,4*nInter});
  GridType grid({1.0,4.0}, resolution); 

  grid.globalRefine(2);
  std::cout << "max Level: " << grid.maxLevel() << std::endl;
  for(unsigned i(0); i < grid.maxLevel() ; ++i) {
    std::stringstream ss; ss << "grid-" << i ;
    write(grid.levelGridView(i), ss.str());
  }
  write(grid.leafGridView(), "leafGrid");
  return 0;
}
