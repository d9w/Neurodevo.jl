#include "external/cxxopts/src/cxxopts.hpp"
#include "core/Viewer.h"
#include "core/config.hpp"
#include "core/DNA.h"
#include "core/ANN.h"
#include "problems/Forage.h"

using std::vector;
using std::cout;
using std::endl;
using std::map;

DNA dna;
Forage forage{0};
ANN ann{dna, 0};
Viewer v;
vector<vector<double>> outputs;
vector<double> max_morphs = {0.01, 0.01, 0.01};

void display() {
  v.DrawAxes();
  v.DrawANN(ann);
}

void step() {
  if (!forage.stop()) {
    auto ins = forage.getInputs();
    ann.step(ins, forage.getReward());
    ann.set_outputs(&outputs);
    forage.step(outputs);
  }

  display();
}

int main(int, char**) {
  v.run(step);
}
