#ifndef COARSE_LAPL_H
#define COARSE_LAPL_H

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrixindexset.hh>
#include <vector>
//represents the application of a coarse laplacian on finer mesh
struct CoarseLapl
{
  typedef Dune::BlockVector<Dune::FieldVector<double,1> > VectorType;
  //constructs the coarse laplacian for the level-0-grid
  template<typename Grid >
  CoarseLapl(Grid& grid, double nu);

  //applys the coarse lapl on the fine data
  void apply(const VectorType& in, VectorType& out) const;

  private:
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > MatrixType;
    //Fine-Coarse (FC) and Coarse-Fine (CF) operator, such that:
    // FC: from fine mesh to this (selection)
    // CF: from coarse mesh to this (interpolation, from fufem)
    struct FCCF
    {
      //takes the coarse-to-fine interpolation from fufem
      FCCF(std::shared_ptr<MatrixType > CF = nullptr):
        CF(CF), size(0)
      {
        if(CF)
          size = CF->N(); 
      }

      //uses the coarse-fine from fufem to create the selection
      void createFC(const MatrixType& cfm)
      {
        if(size > 0 && size != cfm.M())
          throw std::runtime_error("the sizes do not match");
        if(size == 0)
          size = cfm.M();
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

      unsigned getsize() const
      { return size; }

      private:
      std::shared_ptr<MatrixType > CF;
      std::vector<unsigned > FC;
      unsigned size;
    };

    template<typename Grid >
    std::vector<FCCF > getTransferPairs(Grid& grid);
    template <class Basis>
    void getOccupationPattern(const Basis& basis, Dune::MatrixIndexSet& nb);
    template<typename Basis, typename View>
    void fillLapl(Basis& basis, View& gridView, double nu);

    MatrixType lapl;
    std::vector<FCCF > transfers;
};

#endif

