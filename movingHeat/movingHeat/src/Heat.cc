#include "Heat.h"
#include <dune/istl/matrixmarket.hh>
//#include <dune/istl/vector.hh>

#include <dune/grid-glue/extractors/extractorpredicate.hh>
#include <dune/grid-glue/extractors/codim1extractor.hh>

#include <dune/grid-glue/merging/contactmerge.hh>
#include <dune/grid-glue/gridglue.hh>

namespace Helper 
{
  template <class Basis>
  void getOccupationPattern(const Basis& basis, MatrixIndexSet& nb)
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

  template <class GridView>
  struct VerticalFaceDescriptor
  : public GridGlue::ExtractorPredicate<GridView,1>
  {
    virtual bool contains(const typename GridView::template Codim<0>::Entity& element,
        unsigned int face) const
    {
      const int dim = GridView::dimension;
      const auto& refElement = Dune::ReferenceElements<double, dim>::general(element.type());

      // Number of corners of the element face
      int numVertices = refElement.size(face, 1, dim);

      for (int i=0; i<numVertices; i++)
        if ( std::abs(element.geometry().corner(refElement.subEntity(face,1,i,dim))[0] ) > 1e-6 )
          return false;

      return true;
    }
  };
}

Heat::Heat(int nInter, double nu, double alpha, double v0, double source, bool useLapl0):
  L({1.0, 4.0}), lower_mv({-0.5, 2.0}), upper_mv({0.0, 2.5}),
  grid(new GridType(L, std::array<int, dim>({nInter,4*nInter}))),
  hgtmv(lower_mv, upper_mv, std::array<int, 2>({nInter, nInter})),
  grid_mv(new GridType_MV(hgtmv, mf)),
  gridView(grid->leafGridView()), basis(gridView),
  nu(nu), alpha(alpha), v0(v0), sourceVal(source), nInter(nInter), useLapl0(useLapl0)
{ 
  buildMatrices(); 
}

void Heat::setParam(double nu, double alpha)
{
  this->nu=nu;
  this->alpha=alpha;
  buildMatrices(); 
}

void Heat::buildMatrices()
{
  //set nnz structure
  MatrixIndexSet occupationPattern;
  Helper::getOccupationPattern(basis, occupationPattern);
  occupationPattern.exportIdx(mass);
  occupationPattern.exportIdx(lapl);
  if(useLapl0) 
    fillMatricesZeroLapl();
  else
    fillMatrices();
  /*storeMatrixMarket(lapl, "lapl-matrix.mm");
  storeMatrixMarket(mass, "mass-matrix.mm");*/
  mSolver.reset(new MSolver(mass));
}

void Heat::fillMatrices()
{
  lapl = 0.0, mass = 0.0;
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

    Matrix<FieldMatrix<double,1,1> > elemMass, elemLapl;
    //build local mass and laplacian
    using Element = typename Basis::LocalView::Element;
    const int dim = Element::dimension;
    auto geometry = element.geometry();
    const auto& localFiniteElement = localView.tree().finiteElement();
    elemMass.setSize(localFiniteElement.size(),localFiniteElement.size());
    elemMass= 0;      // fills the entire matrix with zeros
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

      //compute the values of the basis functions
      std::vector<FieldVector<double, 1>> basValues;
      localFiniteElement.localBasis().evaluateFunction(qPos, basValues);

      for (size_t i=0; i<elemMass.N(); i++)
        for (size_t j=0; j<elemMass.M(); j++ ) {
          elemLapl[localView.tree().localIndex(i)][localView.tree().localIndex(j)]
            += ( gradients[i] * gradients[j] ) * qP.weight() * integrationElement;
          elemMass[localView.tree().localIndex(i)][localView.tree().localIndex(j)] +=
            basValues[i]*basValues[j]*qP.weight()*integrationElement;
        }
    }

    //move insert entries to global sparse matrix
    for(size_t i=0; i<elemMass.N(); i++)
    {
      // The global index of the i-th degree of freedom of the element
      auto row = localIndexSet.index(i);

      for (size_t j=0; j<elemMass.M(); j++ )
      {
        // The global index of the j-th degree of freedom of the element
        auto col = localIndexSet.index(j);
        mass[row][col] += elemMass[i][j];
        lapl[row][col] += (-nu)*elemLapl[i][j];
      }
    }
  }
}

void Heat::fillMatricesZeroLapl()
{
  lapl = 0.0, mass = 0.0;
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

    Matrix<FieldMatrix<double,1,1> > elemMass, elemLapl;
    //build local mass and laplacian
    using Element = typename Basis::LocalView::Element;
    const int dim = Element::dimension;
    auto geometry = element.geometry();
    const auto& localFiniteElement = localView.tree().finiteElement();
    elemMass.setSize(localFiniteElement.size(),localFiniteElement.size());
    elemMass= 0;      // fills the entire matrix with zeros
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

      //compute the values of the basis functions
      std::vector<FieldVector<double, 1>> basValues;
      localFiniteElement.localBasis().evaluateFunction(qPos, basValues);

      for (size_t i=0; i<elemMass.N(); i++)
        for (size_t j=0; j<elemMass.M(); j++ ) {
#if 1
          if(std::abs(geometry.center()[0]-0.5/nInter) >1e-6) {
            elemLapl[localView.tree().localIndex(i)][localView.tree().localIndex(j)]
              += ( gradients[i] * gradients[j] ) * qP.weight() * integrationElement;
          }
#endif
          elemMass[localView.tree().localIndex(i)][localView.tree().localIndex(j)] +=
            basValues[i]*basValues[j]*qP.weight()*integrationElement;
        }
    }

    //move insert entries to global sparse matrix
    for(size_t i=0; i<elemMass.N(); i++)
    {
      // The global index of the i-th degree of freedom of the element
      auto row = localIndexSet.index(i);

      for (size_t j=0; j<elemMass.M(); j++ )
      {
        // The global index of the j-th degree of freedom of the element
        auto col = localIndexSet.index(j);
        mass[row][col] += elemMass[i][j];
        lapl[row][col] += (-nu)*elemLapl[i][j];
      }
    }
  }
}

void Heat::fastRect(double t, const VectorType& yIn, VectorType& out) const
{
  const double center = 2.25-0.1*t;
  auto rect = [center](const auto& x) { return std::abs(center-x[1]) <= 0.25 ? 1.0 : 0.0 ; };
  //add the neumann part
  fastBoundary(yIn, [&](double uh, const auto& posGlobal) { return rect(posGlobal)*alpha*(v0-uh); }, out);
}

void Heat::fastFull(double t, const VectorType& yIn, VectorType& out) const
{
  fastBoundary(yIn, [&](double uh, const auto& ) { return alpha*(v0-uh); }, out); 
}

template<typename F >
void Heat::fastBoundary(const VectorType& yIn, const F& flux, VectorType& out) const
{
  auto localView = basis.localView();
  auto localIndexSet = basis.localIndexSet();
  for (const auto& bdEl : elements(gridView))
  {
    localView.bind(bdEl);
    localIndexSet.bind(localView);
    for(const auto& inter : intersections(gridView, bdEl))
    {
      //std::cout << "bd center:" << bdEl.geometry().center() << " " << inter.geometry().center() << " " << inter.boundary() << " " << inter.type() <<  std::endl;
      if(inter.boundary() && (std::abs(inter.geometry().center()[0]) < 1.0e-10)) {
        //std::cout << "bd center:" << bdEl.geometry().center() << " " << inter.geometry().center() << " " << inter.boundary() << " " << inter.type() <<  std::endl;
        const int qOrder = 2;
        auto quadRule = QuadratureRules<double, dim-1>::rule(inter.type(), qOrder);
        const auto& localFiniteElement = localView.tree().finiteElement();
        for(const auto& qPos : quadRule) {
          //std::cout << "qPos:" << qPos.position(); // <<std::endl;
          const double det = inter.geometry().integrationElement(qPos.position());
          std::vector<FieldVector<double,1> > shapeFunctionValues;
          //localFiniteElement.localBasis().evaluateFunction(inter.geometryInInside().global(qPos.position()), shapeFunctionValues);
          //std::cout << "\ninInside.type():" <<inter.geometryInInside().type() ;
          //std::cout << "\nqGlobal:" << inter.geometryInInside().global(qPos.position()) ;
          auto posInRef = inter.geometryInInside().global(qPos.position());
          localFiniteElement.localBasis().evaluateFunction(posInRef, shapeFunctionValues); //ist das richtig?? ist das global in Ref, also aus sicht der Linie global
          auto posGlobal = bdEl.geometry().global(posInRef);
          //std::cout << "globalPosBoundary:" << posGlobal << std::endl;
          double uh = 0.0;
          for(unsigned i(0); i < localFiniteElement.size(); ++i) {
            const auto row = localIndexSet.index(i);
            uh += shapeFunctionValues[i]*yIn[row];
          }
          const double fVal = flux(uh, posGlobal); 
          //std::cout << "shapeVals:[";
          for(unsigned i(0); i < localFiniteElement.size(); ++i) {
            //std::cout << shapeFunctionValues[i] << " ";
            const auto row = localIndexSet.index(i);
            out[row] += qPos.weight()*fVal*shapeFunctionValues[i]*det;
          }
          //std::cout << std::endl;
        }

      }
    }
  }
}

void Heat::fastGrid(double t, const VectorType& yIn, VectorType& out) const
{
  if(useLapl0) {
    fastGridZeroLapl(t, yIn, out);
    return; 
  }

  //const double center = 2.25-0.1*t;
  const double dY = -0.1*t;
  //move grid_mv by dY (indirectly through geometrygrid)
  mf.dx[1] = dY;

  typedef GridType_MV::LeafGridView GridView_MV;

  Helper::VerticalFaceDescriptor<GridView> facePredicate0;
  Helper::VerticalFaceDescriptor<GridView_MV> facePredicate1;

  typedef GridGlue::Codim1Extractor<GridView> Extractor0;
  typedef GridGlue::Codim1Extractor<GridView_MV> Extractor1;

  GridGlue::Codim1Extractor<GridView> domEx(grid->leafGridView(), facePredicate0);
  GridGlue::Codim1Extractor<GridView_MV> tarEx(grid_mv->leafGridView(), facePredicate1);

  typedef GridGlue::GridGlue<Extractor0,Extractor1> GlueType;

  // Backend for the computation of the remote intersections
  GridGlue::ContactMerge<dim,double> merger;
  GlueType glue(domEx, tarEx, &merger);

  glue.build();

  const GridView::IndexSet& indexSet0 = grid->leafGridView().indexSet();
  const GridView_MV::IndexSet& indexSet1 = grid_mv->leafGridView().indexSet();

  typedef PQkLocalFiniteElementCache<GridType::ctype, double, dim, 1> TestFECache;
  typedef PQkLocalFiniteElementCache<GridType::ctype, double, dim, 1> FiniteElementCache0;
  FiniteElementCache0 cache0, cache1;
  TestFECache testCache;

  for (const auto& intersection : intersections(glue))
  {
    const FiniteElementCache0::FiniteElementType& nonmortarFiniteElement = cache0.get(intersection.inside().type());
    const FiniteElementCache0::FiniteElementType& mortarFiniteElement    = cache1.get(intersection.outside().type());
    const TestFECache::FiniteElementType&         testFiniteElement      = testCache.get(intersection.inside().type());

    // Select a quadrature rule:  Use order = 2 just for simplicity
    int quadOrder = 2;
    const auto& quad = QuadratureRules<double, dim-1>::rule(intersection.type(), quadOrder);

    // Loop over all quadrature points
    for (size_t l=0; l<quad.size(); l++)
    {
      // compute integration element of overlap
      double integrationElement = intersection.geometry().integrationElement(quad[l].position());

      // quadrature point positions on the reference element
      FieldVector<double,dim> nonmortarQuadPos = intersection.geometryInInside().global(quad[l].position());
      //FieldVector<double,dim> mortarQuadPos    = intersection.geometryInOutside().global(quad[l].position());

      //evaluate all shapefunctions at the quadrature point
      std::vector<FieldVector<double,1> > nonmortarValues,testValues; // mortarValues

      nonmortarFiniteElement.localBasis().evaluateFunction(nonmortarQuadPos,nonmortarValues);
      //mortarFiniteElement   .localBasis().evaluateFunction(mortarQuadPos,mortarValues);
      testFiniteElement     .localBasis().evaluateFunction(nonmortarQuadPos,testValues);

      double uh(0.0);
      for (size_t j=0; j<nonmortarValues.size(); j++) { 
        auto r = indexSet0.subIndex(intersection.inside(), j, dim);
        uh += nonmortarValues[j]*yIn[r];
      }
      const double fVal = (v0-uh)*alpha+sourceVal;
      // Loop over all shape functions of the test space
      for (size_t j=0; j<testFiniteElement.size(); j++)
      {
        int testIdx = indexSet0.subIndex(intersection.inside(),j,dim);
        out[testIdx] += integrationElement*quad[l].weight()*testValues[j]*fVal;

      }

    }

  }
}

void Heat::fastGridZeroLapl(double t, const VectorType& yIn, VectorType& out) const
{
  //const double center = 2.25-0.1*t;
  const double dY = -0.1*t;
  //move grid_mv by dY (indirectly through geometrygrid)
  mf.dx[1] = dY;

  typedef GridType_MV::LeafGridView GridView_MV;

  Helper::VerticalFaceDescriptor<GridView> facePredicate0;
  Helper::VerticalFaceDescriptor<GridView_MV> facePredicate1;

  typedef GridGlue::Codim1Extractor<GridView> Extractor0;
  typedef GridGlue::Codim1Extractor<GridView_MV> Extractor1;

  GridGlue::Codim1Extractor<GridView> domEx(grid->leafGridView(), facePredicate0);
  GridGlue::Codim1Extractor<GridView_MV> tarEx(grid_mv->leafGridView(), facePredicate1);

  typedef GridGlue::GridGlue<Extractor0,Extractor1> GlueType;

  // Backend for the computation of the remote intersections
  GridGlue::ContactMerge<dim,double> merger;
  GlueType glue(domEx, tarEx, &merger);

  glue.build();

  const GridView::IndexSet& indexSet0 = grid->leafGridView().indexSet();
  const GridView_MV::IndexSet& indexSet1 = grid_mv->leafGridView().indexSet();

  typedef PQkLocalFiniteElementCache<GridType::ctype, double, dim, 1> TestFECache;
  typedef PQkLocalFiniteElementCache<GridType::ctype, double, dim, 1> FiniteElementCache0;
  FiniteElementCache0 cache0, cache1;
  TestFECache testCache;

  for (const auto& intersection : intersections(glue))
  {
    const FiniteElementCache0::FiniteElementType& nonmortarFiniteElement = cache0.get(intersection.inside().type());
    const FiniteElementCache0::FiniteElementType& mortarFiniteElement    = cache1.get(intersection.outside().type());
    const TestFECache::FiniteElementType&         testFiniteElement      = testCache.get(intersection.inside().type());

    // Select a quadrature rule:  Use order = 2 just for simplicity
    int quadOrder = 2;
    const auto& quad = QuadratureRules<double, dim-1>::rule(intersection.type(), quadOrder);

    // Loop over all quadrature points
    for (size_t l=0; l<quad.size(); l++)
    {
      // compute integration element of overlap
      double integrationElement = intersection.geometry().integrationElement(quad[l].position());

      // quadrature point positions on the reference element
      FieldVector<double,dim> nonmortarQuadPos = intersection.geometryInInside().global(quad[l].position());
      //FieldVector<double,dim> mortarQuadPos    = intersection.geometryInOutside().global(quad[l].position());

      //evaluate all shapefunctions at the quadrature point
      std::vector<FieldVector<double,1> > nonmortarValues,testValues; // mortarValues

      nonmortarFiniteElement.localBasis().evaluateFunction(nonmortarQuadPos,nonmortarValues);
      //mortarFiniteElement   .localBasis().evaluateFunction(mortarQuadPos,mortarValues);
      testFiniteElement     .localBasis().evaluateFunction(nonmortarQuadPos,testValues);

      double uh(0.0);
      for (size_t j=0; j<nonmortarValues.size(); j++) { 
        auto r = indexSet0.subIndex(intersection.inside(), j, dim);
        uh += nonmortarValues[j]*yIn[r];
      }
      const double fVal = (v0-uh)*alpha;
      // Loop over all shape functions of the test space
      for (size_t j=0; j<testFiniteElement.size(); j++)
      {
        int testIdx = indexSet0.subIndex(intersection.inside(),j,dim);
        out[testIdx] += integrationElement*quad[l].weight()*testValues[j]*fVal;

      }

    }

  }

  addLaplZero(yIn, out);
}

void Heat::addLaplZero(const VectorType& yIn, VectorType& out) const
{
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

    std::vector<double> elemLapl;
    //build local laplacian
    using Element = typename Basis::LocalView::Element;
    const int dim = Element::dimension;
    auto geometry = element.geometry();
    if(std::abs(geometry.center()[0]- 0.5/nInter) < 1e-6) {
      const auto& localFiniteElement = localView.tree().finiteElement();
      elemLapl.resize(localFiniteElement.size());
      for(auto& d : elemLapl) d = 0;      // fills the entire vector with zeros
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

        //compute the values of the basis functions
        //std::vector<FieldVector<double, 1>> basValues;
        //localFiniteElement.localBasis().evaluateFunction(qPos, basValues);

        for (size_t i=0; i<elemLapl.size(); i++)
          for (size_t j=0; j<elemLapl.size(); j++ ) {
              elemLapl[localView.tree().localIndex(i)]
                += yIn[localIndexSet.index(localView.tree().localIndex(j))]*( gradients[i] * gradients[j] ) * qP.weight() * integrationElement;
          }
      }

      //move insert entries to global sparse matrix
      for(size_t i=0; i<elemLapl.size(); i++)
      {
        // The global index of the i-th degree of freedom of the element
        auto row = localIndexSet.index(i);
        // The global index of the j-th degree of freedom of the element
        out[row] += (-nu)*elemLapl[i];
      }
    }
  }
}

void Heat::addFastMatrix(double t, MatrixType& dest) const
{
  //const double center = 2.25-0.1*t;
  const double dY = -0.1*t;
  //move grid_mv by dY (indirectly through geometrygrid)
  mf.dx[1] = dY;

  typedef GridType_MV::LeafGridView GridView_MV;

  Helper::VerticalFaceDescriptor<GridView> facePredicate0;
  Helper::VerticalFaceDescriptor<GridView_MV> facePredicate1;

  typedef GridGlue::Codim1Extractor<GridView> Extractor0;
  typedef GridGlue::Codim1Extractor<GridView_MV> Extractor1;

  GridGlue::Codim1Extractor<GridView> domEx(grid->leafGridView(), facePredicate0);
  GridGlue::Codim1Extractor<GridView_MV> tarEx(grid_mv->leafGridView(), facePredicate1);

  typedef GridGlue::GridGlue<Extractor0,Extractor1> GlueType;

  // Backend for the computation of the remote intersections
  GridGlue::ContactMerge<dim,double> merger;
  GlueType glue(domEx, tarEx, &merger);

  glue.build();

  const GridView::IndexSet& indexSet0 = grid->leafGridView().indexSet();
  const GridView_MV::IndexSet& indexSet1 = grid_mv->leafGridView().indexSet();

  typedef PQkLocalFiniteElementCache<GridType::ctype, double, dim, 1> TestFECache;
  typedef PQkLocalFiniteElementCache<GridType::ctype, double, dim, 1> FiniteElementCache0;
  FiniteElementCache0 cache0, cache1;
  TestFECache testCache;

  for (const auto& intersection : intersections(glue))
  {
    const FiniteElementCache0::FiniteElementType& nonmortarFiniteElement = cache0.get(intersection.inside().type());
    const FiniteElementCache0::FiniteElementType& mortarFiniteElement    = cache1.get(intersection.outside().type());
    const TestFECache::FiniteElementType&         testFiniteElement      = testCache.get(intersection.inside().type());

    // Select a quadrature rule:  Use order = 2 just for simplicity
    int quadOrder = 2;
    const auto& quad = QuadratureRules<double, dim-1>::rule(intersection.type(), quadOrder);

    // Loop over all quadrature points
    for (size_t l=0; l<quad.size(); l++)
    {
      // compute integration element of overlap
      double integrationElement = intersection.geometry().integrationElement(quad[l].position());

      // quadrature point positions on the reference element
      FieldVector<double,dim> nonmortarQuadPos = intersection.geometryInInside().global(quad[l].position());
      //FieldVector<double,dim> mortarQuadPos    = intersection.geometryInOutside().global(quad[l].position());

      //evaluate all shapefunctions at the quadrature point
      std::vector<FieldVector<double,1> > nonmortarValues,testValues; // mortarValues

      nonmortarFiniteElement.localBasis().evaluateFunction(nonmortarQuadPos,nonmortarValues);
      //mortarFiniteElement   .localBasis().evaluateFunction(mortarQuadPos,mortarValues);
      testFiniteElement     .localBasis().evaluateFunction(nonmortarQuadPos,testValues);

      /*double uh(0.0);
      for (size_t j=0; j<nonmortarValues.size(); j++) { 
        auto r = indexSet0.subIndex(intersection.inside(), j, dim);
        uh += nonmortarValues[j]*yIn[r];
      }
      const double fVal = (5.0-uh)*alpha;*/
      // Loop over all shape functions of the test space
      for (size_t j=0; j<testFiniteElement.size(); j++)
      {
        int testR = indexSet0.subIndex(intersection.inside(),j,dim);
        for (size_t c=0; c<testFiniteElement.size(); c++) {
          int basC = indexSet0.subIndex(intersection.inside(),c,dim);
          dest[testR][basC] -= alpha*integrationElement*quad[l].weight()*testValues[j]*nonmortarValues[c];
        }
      }

    }

  }
}

void Heat::writeResult(std::string fname, const VectorType& res)
{
  VTKWriter<GridView> vtkWriter(gridView);
  vtkWriter.addVertexData(res, "solution");
  vtkWriter.write(fname);
}

Heat::MatrixType Heat::getLaplaceWithMove(double t)
{
  MatrixType lwm;
  lwm = lapl;
  addFastMatrix(t, lwm);
  return lwm;
}
