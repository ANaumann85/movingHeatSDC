
#include <dune/geometry/quadraturerules.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>
#include <dune/fufem/assemblers/transferoperatorassembler.hh>

using namespace Dune;

template<typename Grid >
CoarseLapl::CoarseLapl(Grid& grid, double nu)
{
  typedef typename Grid::LevelGridView LevelView;
  typedef typename Functions::PQkNodalBasis<LevelView,1> LevelBasis;
  transfers = getTransferPairs(grid);
  LevelView view(grid.levelGridView(0));
  LevelBasis basis(view);
  fillLapl(basis, view, nu);
}

template<typename Grid >
std::vector<CoarseLapl::FCCF > CoarseLapl::getTransferPairs(Grid& grid)
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

template<typename Basis, typename View>
void CoarseLapl::fillLapl(Basis& basis, View& gridView, double nu)
{
  {
  MatrixIndexSet occupationPattern;
  getOccupationPattern(basis, occupationPattern);
  occupationPattern.exportIdx(lapl);
  }
  lapl = 0.0;
  auto localView = basis.localView();
  auto localIndexSet = basis.localIndexSet();

  for (const auto& element : elements(gridView))             /*@\label{li:poissonequation_elementloop}@*/
  {
    // { assembler_element_loop_end }

    // Now let's get the element stiffness matrix
    // A dense matrix is used for the element stiffness matrix
    // { assembler_assemble_element_matrix_begin }
    localView.bind(element);
    localIndexSet.bind(localView);

    Matrix<FieldMatrix<double,1,1> > elemLapl;
    //build local mass and laplacian
    using Element = typename Basis::LocalView::Element;
    const int dim = Element::dimension;
    auto geometry = element.geometry();
    const auto& localFiniteElement = localView.tree().finiteElement();
    elemLapl.setSize(localFiniteElement.size(),localFiniteElement.size());
    elemLapl= 0;      // fills the entire matrix with zeros
    int order = 2*(dim*localFiniteElement.localBasis().order());       
    const auto& quadRule = QuadratureRules<double, dim>::rule(element.type(), order); 
    for(auto& qP : quadRule) {
      const auto qPos = qP.position();
      // The transposed inverse Jacobian of the map from the reference element to the element
      const auto jacobian = geometry.jacobianInverseTransposed(qPos);

      // The multiplicative factor in the integral transformation formula
      const auto integrationElement = geometry.integrationElement(qPos);
      std::vector<FieldMatrix<double,1,dim> > referenceGradients;
      localFiniteElement.localBasis().evaluateJacobian(qPos, referenceGradients);

      // Compute the shape function gradients on the real element
      std::vector<FieldVector<double,dim> > gradients(referenceGradients.size());
      for (size_t i=0; i<gradients.size(); i++)
        jacobian.mv(referenceGradients[i][0], gradients[i]);

      for (size_t i=0; i<elemLapl.N(); i++)
        for (size_t j=0; j<elemLapl.M(); j++ ) {
          elemLapl[localView.tree().localIndex(i)][localView.tree().localIndex(j)]
            += ( gradients[i] * gradients[j] ) * qP.weight() * integrationElement;
        }
    }

    //move insert entries to global sparse matrix
    for(size_t i=0; i<elemLapl.N(); i++)
    {
      // The global index of the i-th degree of freedom of the element
      auto row = localIndexSet.index(i);

      for (size_t j=0; j<elemLapl.M(); j++ )
      {
        // The global index of the j-th degree of freedom of the element
        auto col = localIndexSet.index(j);
        lapl[row][col] += (-nu)*elemLapl[i][j];
      }
    }
  }
}
template <class Basis>
void CoarseLapl::getOccupationPattern(const Basis& basis, MatrixIndexSet& nb)
{
  nb.resize(basis.size(), basis.size());

  auto gridView = basis.gridView();

  // A loop over all elements of the grid
  auto localView = basis.localView();
  auto localIndexSet = basis.localIndexSet();

  for (const auto& element : elements(gridView))
  {
    localView.bind(element);
    localIndexSet.bind(localView);

    for(size_t i=0; i<localIndexSet.size(); i++)
    {
      // The global index of the i-th vertex of the element
      auto row = localIndexSet.index(i);

      for (size_t j=0; j<localIndexSet.size(); j++ )
      {
        // The global index of the j-th vertex of the element
        auto col = localIndexSet.index(j);
        nb.add(row,col);
      }
    }
  }
}
