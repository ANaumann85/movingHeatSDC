#include "Heat.h"
#include <sstream>
#include <fstream>
#include <string>

int main(int argc, char* argv[])
{
  if(argc < 2) {
    std::cerr << "usage: " << argv[0] << " <nStep>\n";
    return 1;
  }
  //mpi-helper from dune
  MPIHelper::instance(argc, argv);
  Heat heat(10);
  Heat::VectorType dest, yIn;
  heat.init(dest);
  heat.init(yIn); yIn = 0.0;

  double t0=0.0;
  double te=20.0;
  unsigned nStep(40);
  {std::stringstream ss; ss << argv[1]; ss >> nStep; }
  const unsigned int id = 99;

  double dt=(te-t0)/nStep;

  std::string fname;
  {std::stringstream ss; ss << "pointSrc_" << nStep << ".dat"; fname = ss.str();}
  std::fstream file(fname, std::ios_base::trunc | std::ios_base::out);
  for(unsigned i(0); i < nStep; ++i, t0+=dt) {
    heat.fast(t0, yIn, dest);
    file << t0 << " " << dest[id] << std::endl;
  }
  file.close();
  return 0;
}

