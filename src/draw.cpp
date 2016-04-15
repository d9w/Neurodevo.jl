#include <random>
#include "external/cxxopts/src/cxxopts.hpp"
#include "core/config.hpp"
#include "problems/forage.hpp"
#include "core/evaluator.hpp"
#include "core/types.hpp"
#include "core/Viewer.h"
//#include "core/environment_viewer.hpp"

using std::vector;
using std::cout;
using std::endl;
using std::map;

void stepfun() {
  std::cout << "step" << std::endl;
}

int main(int argc, char** argv) {
  using dna_t = Types::DNAType;
  Types t;

  dna_t dna = t.random_dna();

  Viewer v;
  v.run();

  /*
  Forage forager(0);
  Evaluator<Forage> eval;
  eval.evaluate(forager, dna, 0);
  for (auto& fit : *eval.getFitnesses()) cout << fit.first << " : " << fit.second << endl;
  */
}
