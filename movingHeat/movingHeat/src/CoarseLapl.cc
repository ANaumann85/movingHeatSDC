#include <config.h>
#include "CoarseLapl.h"

#include <iostream>
void CoarseLapl::apply(const VectorType& in, VectorType& out) const
{
  //coarsen fine, coarsend input in c1
  VectorType c1(in.size()), c2(in.size());
  c2=in;
  {
    auto curLevel(transfers.rbegin()); 
    curLevel++; // pointer to finest level -> nothing to do
    for( ; curLevel != transfers.rend(); ++curLevel) {
      c1.resize(curLevel->getsize());
      curLevel->applyFC(c2,c1);
      c2.resize(c1.size());
      c2=c1;
    }
  }

  //apply lapl
  lapl.mv(c1,c2);

  //refine lapl
  {
    auto curfLevel(transfers.begin()); //pointer to lowest level -> nothing todo
    curfLevel++;
    for( ; curfLevel != transfers.end(); ++curfLevel) {
      std::cout << "fine size: " << curfLevel->getsize() << std::endl;
      c1.resize(curfLevel->getsize());
      curfLevel->applyCF(c2, c1);
      c2.resize(c1.size());
      c2 = c1;
    }
  }
  out = c1;
}
