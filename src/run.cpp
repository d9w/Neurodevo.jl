#include "external/cxxopts/src/cxxopts.hpp"
#include "core/config.hpp"
#include "problems/forage.hpp"
#include "core/evaluator.hpp"
#include "core/types.hpp"
#include <random>

using std::vector;
using std::cout;
using std::endl;
using std::map;

class Sample {
public:
  using d2 = vector<vector<double>>;

protected:
  int iter;

public:
  std::mt19937 rd;

  Sample(int seed) {
    rd = mt19937(seed);
    iter = 0;
  }

  bool stop() {
    return (iter > 10);
  }

  void setInputs(d2* inputs) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    inputs->clear();
    for (unsigned int x = 0; x < Config::X_SIZE; x++) {
      vector<double> i;
      for (unsigned int y = 0; y < Config::Y_SIZE; y++) {
        i.push_back(dist(rd));
      }
      inputs->push_back(i);
    }
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
  t.addRandomProtein(ProteinType::input, "x");
  t.addRandomProtein(ProteinType::input, "y");
  t.addRandomProtein(ProteinType::input, "z");
  t.addRandomProtein(ProteinType::input, "nt");
  t.addRandomProtein(ProteinType::input, "comm");
  t.addRandomProtein(ProteinType::input, "div");
  t.addRandomProtein(ProteinType::input, "reward");

  t.addRandomProtein(ProteinType::output, "nt");
  t.addRandomProtein(ProteinType::output, "nt_t");
  t.addRandomProtein(ProteinType::output, "f");
  t.addRandomProtein(ProteinType::output, "f_t");
  t.addRandomProtein(ProteinType::output, "comm");

  t.randomReguls(1);
  t.randomParams();

  Forage forager(0);
  Evaluator<Forage> eval;
  eval.evaluate(forager, t, 0);
  for (auto& fit : *eval.getFitnesses()) cout << fit.first << " : " << fit.second << endl;
}
