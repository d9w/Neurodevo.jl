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

  std::cout << "fitnesses: ";
  auto fitnesses = forage.getFitness();
  for (auto& fit : fitnesses) {
    std::cout << fit.first << ":" << fit.second << " ";
  }
  std::cout << std::endl << "footprint: ";
  auto footprint = forage.getFootprint();
  for (auto& foot : footprint) {
    for (auto& i : foot) {
      std::cout << i << " ";
    }
  }
  return 0;
}
