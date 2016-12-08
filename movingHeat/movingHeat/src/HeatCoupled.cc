#include "HeatCoupled.h"
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

HeatCoupled::HeatCoupled(int nInter, double nu, double alpha, double v0, double source, bool useLapl0, bool addConstRobin, int useLaplTilde):
  L({1.0, 4.0}), lower_mv({-0.5, 1.75}), upper_mv({0.0, 2.25}),
  grid(new GridType(L, std::array<int, dim>({nInter,4*nInter}))),
  hgtmv(lower_mv, upper_mv, std::array<int, 2>({nInter, nInter})),
  grid_mv(new GridType_MV(hgtmv, mf)),
  gridView(grid->leafGridView()), gridView_mv(grid_mv->leafGridView()),
  basis(gridView), basis_mv(gridView_mv),
  nu(nu), alpha(alpha), v0(v0), sourceVal(source), 
  bAlph(10.0), bVal(1.0), nInter(nInter), 
  useLapl0(useLapl0), addConstRobin(addConstRobin),useSlowExpl(false)
{ 
  double h(1.0/nInter);
  switch(useLaplTilde) {
    case 0:
      h=0.0;
      break;
    case 1:
      break;
    case 2:
      h*=h;
      break;
    default:
      throw std::runtime_error("only 0(off), 1(h), 2(h^2) allowed");
  }
  this->useLaplTilde = useLaplTilde > 0;
  buildMatrices(h); 
}

void HeatCoupled::setParam(double nu, double alpha)
{
  this->nu=nu;
  this->alpha=alpha;
  if(useLaplTilde)
    throw std::runtime_error(" missing stepsize at this point ");
  buildMatrices(0.0); 
}

void HeatCoupled::buildMatrices(double h)
{
  //set nnz structure
  MatrixIndexSet occupationPattern, occupationPattern_mv;
  Helper::getOccupationPattern(basis, occupationPattern);
  Helper::getOccupationPattern(basis_mv, occupationPattern_mv);
  occupationPattern.exportIdx(mass[0]);
  occupationPattern.exportIdx(lapl[0]);
  occupationPattern_mv.exportIdx(mass[1]);
  occupationPattern_mv.exportIdx(lapl[1]);
  if(useLapl0) {
    //throw std::runtime_error("not yet supported");
    occupationPattern.exportIdx(lapl0[0]);
    fillMatricesZeroLapl();
  }
  else
    fillMatrices();
  /*storeMatrixMarket(lapl, "lapl-matrix.mm");
  storeMatrixMarket(mass, "mass-matrix.mm");*/
  for(unsigned i(0); i < 2; ++i)
    mSolver[i].reset(new MSolver(mass[i]));
  if(addConstRobin)
    setConstRobin();
  if(useLaplTilde) {
    laplTilde0 = lapl[0];
    laplTilde0 *= h;
    laplTilde1 = lapl[1];
    //laplTilde1 *= h*0.25;
    laplTilde1 *= 0.0;
  }
  //fillMovingExchangeMass();
}

void HeatCoupled::setConstRobin()
{
  constRobM = mass[0]; constRobM = 0.0;
  constRobB[0].resize(constRobM.N()); constRobB[0]=0.0;

  auto localView = basis.localView();
  auto localIndexSet = basis.localIndexSet();
  const double bAlph = this->bAlph;
  const double bVal = this->bVal; 
  for (const auto& bdEl : elements(gridView))
  {
    localView.bind(bdEl);
    localIndexSet.bind(localView);
    for(const auto& inter : intersections(gridView, bdEl))
    {
      //std::cout << "bd center:" << bdEl.geometry().center() << " " << inter.geometry().center() << " " << inter.boundary() << " " << inter.type() <<  std::endl;
      if(inter.boundary() && (std::abs(inter.geometry().center()[0]-1.0) < 1.0e-8)) {
        //std::cout << "bd center:" << bdEl.geometry().center() << " " << inter.geometry().center() << " " << inter.boundary() << " " << inter.type() <<  std::endl;
        const int qOrder = 3;
        auto quadRule = QuadratureRules<double, dim-1>::rule(inter.type(), qOrder);
        const auto& localFiniteElement = localView.tree().finiteElement();
        for(const auto& qPos : quadRule) {
          //std::cout << "qPos:" << qPos.position(); // <<std::endl;
          const double det = inter.geometry().integrationElement(qPos.position());
          std::vector<FieldVector<double,1> > shapeFunctionValues;
          localFiniteElement.localBasis().evaluateFunction(inter.geometryInInside().global(qPos.position()), shapeFunctionValues);
          //std::cout << "\ninInside.type():" <<inter.geometryInInside().type() ;
          //std::cout << "\nqGlobal:" << inter.geometryInInside().global(qPos.position()) ;
          //std::cout << "shapeVals:[";
          for(unsigned i(0); i < localFiniteElement.size(); ++i) {
            //std::cout << shapeFunctionValues[i] << " ";
            const auto row = localIndexSet.index(i);
            constRobB[0][row] += qPos.weight()*bAlph*bVal*shapeFunctionValues[i]*det;
            for(unsigned j(0); j < localFiniteElement.size(); ++j) {
              const auto col = localIndexSet.index(j);
              constRobM[row][col] += -qPos.weight()*bAlph*shapeFunctionValues[j]*shapeFunctionValues[i]*det;
            }
          }
          //std::cout << std::endl;
        }

      }
    }
  }
}

void HeatCoupled::fillMovingExchangeMass()
{
  movingExchangeMass = mass[1]; movingExchangeMass = 0.0;

  auto localView = basis_mv.localView();
  auto localIndexSet = basis_mv.localIndexSet();
  const double alpha = -this->alpha;
  for (const auto& bdEl : elements(gridView_mv))
  {
    localView.bind(bdEl);
    localIndexSet.bind(localView);
    for(const auto& inter : intersections(gridView_mv, bdEl))
    {
      //std::cout << "bd center:" << bdEl.geometry().center() << " " << inter.geometry().center() << " " << inter.boundary() << " " << inter.type() <<  std::endl;
      if(inter.boundary() && (std::abs(inter.geometry().center()[0]) < 1.0e-6)) {
        //std::cout << "bd center:" << bdEl.geometry().center() << " " << inter.geometry().center() << " " << inter.boundary() << " " << inter.type() <<  std::endl;
        const int qOrder = 2; //just for consistency...
        auto quadRule = QuadratureRules<double, dim-1>::rule(inter.type(), qOrder);
        const auto& localFiniteElement = localView.tree().finiteElement();
        for(const auto& qPos : quadRule) {
          //std::cout << "qPos:" << qPos.position(); // <<std::endl;
          const double det = inter.geometry().integrationElement(qPos.position());
          std::vector<FieldVector<double,1> > shapeFunctionValues;
          localFiniteElement.localBasis().evaluateFunction(inter.geometryInInside().global(qPos.position()), shapeFunctionValues);
          //std::cout << "\ninInside.type():" <<inter.geometryInInside().type() ;
          //std::cout << "\nqGlobal:" << inter.geometryInInside().global(qPos.position()) ;
          //std::cout << "shapeVals:[";
          for(unsigned i(0); i < localFiniteElement.size(); ++i) {
            //std::cout << shapeFunctionValues[i] << " ";
            const auto row = localIndexSet.index(i);
            for(unsigned j(0); j < localFiniteElement.size(); ++j) {
              const auto col = localIndexSet.index(j);
              movingExchangeMass[row][col] += qPos.weight()*alpha*shapeFunctionValues[j]*shapeFunctionValues[i]*det;
            }
          }
          //std::cout << std::endl;
        }

      }
    }
  }
}

namespace Helper {
using MT = HeatCoupled::MatrixType;
template<typename Basis, typename View>
void fillMatrices(MT& mass, MT& lapl, Basis& basis, View& gridView)
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
        lapl[row][col] += elemLapl[i][j];
      }
    }
  }
} 

}
void HeatCoupled::fillMatrices()
{
  using Helper::fillMatrices;
  for(unsigned i(0); i < 2; ++i) {
    lapl[i] = 0.0, mass[i] = 0.0;
  }
  fillMatrices(mass[0], lapl[0], basis, gridView);
  fillMatrices(mass[1], lapl[1], basis_mv, gridView_mv);
  for(auto& l : lapl)  l*= -nu;
}

void HeatCoupled::fillMatricesZeroLapl()
{
  using Helper::fillMatrices;
#if 1
  lapl[0] = 0.0, mass[0] = 0.0;
  lapl[1] = 0.0, mass[1] = 0.0;
  //classical filling of moving part
  fillMatrices(mass[1], lapl[1], basis_mv, gridView_mv);

  //zero-lapl-filling of non-moving part
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
    if(geometry.center()[0] > 0.5/nInter) {
      for(size_t i=0; i<elemMass.N(); i++)
      {
        // The global index of the i-th degree of freedom of the element
        auto row = localIndexSet.index(i);

        for (size_t j=0; j<elemMass.M(); j++ )
        {
          // The global index of the j-th degree of freedom of the element
          auto col = localIndexSet.index(j);
          mass[0][row][col] += elemMass[i][j];
          lapl[0][row][col] += (-nu)*elemLapl[i][j];
        }
      }
    } else {
      for(size_t i=0; i<elemMass.N(); i++)
      {
        // The global index of the i-th degree of freedom of the element
        auto row = localIndexSet.index(i);

        for (size_t j=0; j<elemMass.M(); j++ )
        {
          // The global index of the j-th degree of freedom of the element
          auto col = localIndexSet.index(j);
          mass[0][row][col] += elemMass[i][j];
          lapl0[0][row][col] += (-nu)*elemLapl[i][j];
        }
      }
    }
  }
#else
  throw std::runtime_error("not yet supported");
#endif
}

#if 0
void HeatCoupled::fastRect(double t, const VectorType& yIn, VectorType& out) const
{
  const double center = 2.25-0.1*t;
  auto rect = [center](const auto& x) { return std::abs(center-x[1]) <= 0.25 ? 1.0 : 0.0 ; };
  //add the neumann part
  fastBoundary(yIn, [&](double uh, const auto& posGlobal) { return rect(posGlobal)*alpha*(v0-uh); }, out);
}

void HeatCoupled::fastFull(double t, const VectorType& yIn, VectorType& out) const
{
  fastBoundary(yIn, [&](double uh, const auto& ) { return alpha*(v0-uh); }, out); 
}

template<typename F >
void HeatCoupled::fastBoundary(const VectorType& yIn, const F& flux, VectorType& out) const
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
#endif

void HeatCoupled::fastGrid(double t, const VectorType& yIn, VectorType& out, unsigned tag) const
{
#if 0
  switch(tag) {
    case 1:
      std::cout << "only 0\n";
      break;
    case 2:
      std::cout << "only 1\n";
      break;
    case 3:
      std::cout << "both\n";
      break;
  }
#endif
  //const double center = 2.25-0.1*t;
  //const double dY = -0.1*t;
  const double dY = -sin(M_PI/10.0*t)*1.75;
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
      FieldVector<double,dim> mortarQuadPos    = intersection.geometryInOutside().global(quad[l].position());

      //evaluate all shapefunctions at the quadrature point
      std::vector<FieldVector<double,1> > nonmortarValues,testValues, mortarValues;

      nonmortarFiniteElement.localBasis().evaluateFunction(nonmortarQuadPos,nonmortarValues);
      mortarFiniteElement   .localBasis().evaluateFunction(mortarQuadPos,mortarValues);
      testFiniteElement     .localBasis().evaluateFunction(nonmortarQuadPos,testValues);

      double uh(0.0);
      for (size_t j=0; j<nonmortarValues.size(); j++) { 
        auto r = indexSet0.subIndex(intersection.inside(), j, dim);
        uh += nonmortarValues[j]*yIn[0][r];
      }
      double vh(0.0);
      for (size_t j=0; j<mortarValues.size(); j++) { 
        auto r = indexSet1.subIndex(intersection.outside(), j, dim);
        vh += mortarValues[j]*yIn[1][r];
      }
      const double fVal0 = (vh-uh)*alpha+sourceVal;
      // Loop over all shape functions of the test space
      if((tag & 1) > 0) {
        for (size_t j=0; j<testFiniteElement.size(); j++)
        {
          int testIdx = indexSet0.subIndex(intersection.inside(),j,dim);
          out[0][testIdx] += integrationElement*quad[l].weight()*testValues[j]*fVal0;
        }
      } 
      double fVal1 ; 
      if(movingExchangeMass.N() > 0) {
        fVal1 = uh*alpha+sourceVal;
      } else { 
        fVal1 = -(vh-uh)*alpha+sourceVal;
      }
      if((tag & 2) > 0) {
        for (size_t j=0; j<mortarFiniteElement.size(); j++)
        {
          int testIdx = indexSet1.subIndex(intersection.outside(),j,dim);
          out[1][testIdx] += integrationElement*quad[l].weight()*mortarValues[j]*fVal1;
        }
      }

    }

  }
  if(useLapl0) {
    lapl0[0].umv(yIn[0], out[0]);
    //throw std::runtime_error("not yet supported");
    //fastGridZeroLapl(t, yIn, out);
  }
}

#if 0
void HeatCoupled::fastGridZeroLapl(double t, const VectorType& yIn, VectorType& out) const
{
#if 0
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

  lapl0.umv(yIn, out);
#else
  throw std::runtime_error("not supported");
#endif
}
#endif

void HeatCoupled::addFastMatrix(double t, MatrixType& dest) const
{
#if 1
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
  if(useLapl0) {
    dest += lapl0[0];

  }
#else
  throw std::runtime_error("not yet supported");
#endif
}

void HeatCoupled::writeResult(std::string fname, std::string fname_mv, const VectorType& res)
{
  VTKWriter<GridView> vtkWriter(gridView);
  vtkWriter.addVertexData(res[0], "solution");
  vtkWriter.write(fname);

  VTKWriter<GridView_MV> vtkWriter_mv(gridView_mv);
  vtkWriter_mv.addVertexData(res[1], "solution");
  vtkWriter_mv.write(fname_mv);
}

HeatCoupled::MatrixType HeatCoupled::getLaplaceWithMove(double t)
{
  MatrixType lwm;
  lwm = lapl[0];
  addFastMatrix(t, lwm);
  return lwm;
}
