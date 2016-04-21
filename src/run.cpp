#include "external/cxxopts/src/cxxopts.hpp"
#include "core/config.hpp"
#include "core/DNA.h"
#include "core/ANN.h"
#include "problems/Forage.h"

int main(int, char**) {
  DNA dna;
  Forage forage(0);
  ANN ann(dna, 0);
  vector<vector<double>> outputs;

  while (!forage.stop()) {
    auto ins = forage.getInputs();
    ann.step(ins, forage.getReward());
    ann.set_outputs(&outputs);
    forage.step(outputs);
  }
 return 0;
}
