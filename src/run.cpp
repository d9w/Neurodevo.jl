#include "external/cxxopts/src/cxxopts.hpp"
#include "core/config.hpp"
#include "core/evaluator.hpp"
#include "core/types.hpp"

using std::vector;
using std::cout;
using std::endl;
using std::map;

class Sample {
public:
  using d2 = vector<vector<double>>;

protected:
  int iter = 0;

public:
  Sample() {}

  bool stop() {
    return (iter > 10);
  }

  void setInputs(d2* inputs) {
    inputs->clear();
    inputs->push_back(vector<double> {0});
  }

  void step(const d2* outputs, d2* history) {
    cout << "Outputs: ";
    for (auto& v : *outputs)
      for (auto& o : v) cout << o << " ";
    cout << endl;
    history->push_back(vector<double> {1});
    iter++;
  }

  double getReward() {
    return 0.0;
  }

  void setFitness(map<string,double>* fitnesses) {
    fitnesses->clear();
    fitnesses->operator[]("example") = 0.0;
  }
};

int main(int argc, char** argv) {
  using dna_t = Types::DNAType;

  dna_t t;
  Evaluator<Sample> eval;
  eval.evaluate(t);
  for (auto& fit : *eval.getFitnesses()) cout << fit.first << " : " << fit.second << endl;
}
