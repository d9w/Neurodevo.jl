#include "DNA.h"
#include "../external/grgen/classic.hpp"
#include "config.hpp"

DNA::DNA() {
  inputs = {"x","y","z","nt","threshold","weight","comm","div","reward"};
  for (unsigned int i = 0; i < Config::N_M; i++) {
    inputs.push_back("m"+std::to_string(i));
  }

  for (auto& i : inputs) grn.addRandomProtein(ProteinType::input, i);

  outputs = {"nt","nt_t","f","f_t","comm","axon_div","axon_die","axon_none"};
  for (unsigned int i = 0; i < Config::N_M; i++) {
    outputs.push_back("m"+std::to_string(i));
  }
  outputs.push_back("m_thresh");

  for (auto& i : outputs) grn.addRandomProtein(ProteinType::output, i);

  //grn.randomReguls(Config::GRN_REGULS);
  grn.randomReguls(0);
  grn.randomParams();
}

const DNA DNA::random() {
  DNA dna;
  return dna;
}

std::string DNA::concString() const {
  std::string out;
  for (auto& i : inputs) {
    out += std::to_string(grn.getProteinConcentration(i, ProteinType::input)) + " ";
  }
  for (auto& i : outputs) {
    out += std::to_string(grn.getProteinConcentration(i, ProteinType::output)) + " ";
  }
  return out;
}
