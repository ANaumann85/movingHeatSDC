#include <config.h>
#include <dune/grid/yaspgrid.hh>
#include <dune/common/function.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/fufem/assemblers/transferoperatorassembler.hh>

#include <array>
#include <string>
#include <sstream>
#include <vector>

using namespace Dune;
typedef BlockVector<FieldVector<double,1> > VectorType;
typedef BCRSMatrix<FieldMatrix<double,1,1> > MatrixType;

//Fine-Coarse (FC) and Coarse-Fine (CF) operator, such that:
// FC: from fine mesh to this (selection)
// CF: from coarse mesh to this (interpolation)
struct FCCF
{
  //takes the coarse-to-fine interpolation from fufem
  FCCF(std::shared_ptr<MatrixType > CF = nullptr):
    CF(CF)
  {}

  //uses the coarse-fine from fufem to create the selection
  void createFC(const MatrixType& cfm)
  {
    std::cout << "create FC: " << cfm.N() << " " << cfm.M() << std::endl;
    FC.resize(cfm.M());
    auto rbeg(cfm.begin());
    auto rend(cfm.end());
    for( ; rbeg != rend ; ++rbeg)
    {
      auto cbeg(rbeg->begin());
#if 0
      std::cout << "cur row: " << rbeg.index() << " " << rbeg->getsize();
      auto cend(rbeg->end());
      for( ; cbeg != cend; ++cbeg)
        std::cout << " " << *cbeg ;
      std::cout << std::endl;
#endif
      if(rbeg->getsize() == 1 && *cbeg == 1.0) {
        //std::cout << "adding " << rbeg.index() << " " << cbeg.index() << std::endl;
        FC[cbeg.index()] = rbeg.index();
      }
    }
  }

  void applyFC(const VectorType& in, VectorType& out) const
  { 
    if(FC.size() == out.size()) {
      for(unsigned i(0); i < FC.size(); ++i) {
        out[i] = in[FC[i]];
      }
    }
  }

  void applyCF(const VectorType& in, VectorType& out) const
  { if(CF) CF->mv(in, out); }

  private:
  std::shared_ptr<MatrixType > CF;
  std::vector<unsigned > FC;
};

template<typename Grid >
std::vector<FCCF > getTransferPairs(Grid& grid)
{
  typedef TransferOperatorAssembler<Grid> TransferOperator;
  std::vector<std::shared_ptr<MatrixType> > matrices;
  TransferOperator transOp(grid);
  transOp.assembleMatrixHierarchy(matrices);
  //one pair for each level
  std::vector<FCCF > ret(matrices.size()+1);
  ret[0].createFC(*matrices[0]);
  for(unsigned i(1); i < matrices.size(); ++i) {
    ret[i] = FCCF(matrices[i-1]);
    ret[i].createFC(*matrices[i]);
  }
  ret[matrices.size()] = FCCF(matrices.back());
  return ret;
}

int main(int argc, char* argv[])
{
  MPIHelper::instance(argc, argv);
  static const int dim = 2;
  typedef YaspGrid<dim > GridType;
  typedef typename GridType::LeafGridView LeafView;
  typedef typename GridType::LevelGridView LevelView;
  typedef Functions::PQkNodalBasis<LeafView,1> LeafBasis;
  typedef Functions::PQkNodalBasis<LevelView,1> LevelBasis;


  auto func = [](const auto& x) { return 0.1*(x[0]*x[0]+x[1]*x[1]); };
  int nInter(5);
  std::array<int, dim> resolution({nInter,4*nInter});
  GridType grid({1.0,4.0}, resolution); 

  grid.globalRefine(2);
  auto matrices(getTransferPairs(grid));
  std::cout << "max Level: " << grid.maxLevel() << std::endl;
  std::cout << "number of matrices: " << matrices.size() << std::endl;

  std::vector<VectorType> data(grid.maxLevel()+1), dataTransferCF(grid.maxLevel()+1), dataTransferFC(grid.maxLevel());
  for(unsigned i(0); i <= grid.maxLevel() ; ++i) {
    //std::cout << "matrix size: " << matrices[i]->N() << "x" << matrices[i]->M() << std::endl;
    LevelBasis basis(grid.levelGridView(i));
    data[i].resize(basis.size());

    Functions::interpolate(basis, data[i], func);
  }
  
  /*LeafBasis basis(grid.leafGridView());
  data[grid.maxLevel()].resize(basis.size());
  Functions::interpolate(basis, data[grid.maxLevel()], func);*/
  //write(grid.leafGridView(), "leafGrid", data);

  for(unsigned i(0); i <= grid.maxLevel() ; ++i) {
    if(i != grid.maxLevel()) {
      dataTransferCF[i+1].resize(data[i+1].size());
      matrices[i+1].applyCF(data[i], dataTransferCF[i+1]);
      dataTransferFC[i].resize(data[i].size());
      matrices[i].applyFC(data[i+1], dataTransferFC[i]);
    }
    std::stringstream ss; ss << "grid-" << i ;
    VTKWriter<LevelView> vtkWriter(grid.levelGridView(i));
    vtkWriter.addVertexData(data[i], "T");
    if(i == 0) {
      vtkWriter.addVertexData(dataTransferFC[i], "TFC");
    }
    if(i==1) {
      vtkWriter.addVertexData(dataTransferFC[i], "TFC");
      vtkWriter.addVertexData(dataTransferCF[i], "TCF");
    }
    if(i==2) {
      vtkWriter.addVertexData(dataTransferCF[i], "TCF");
    }
    vtkWriter.write(ss.str());
  }
  return 0;
}
