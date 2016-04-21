#include "DNA.h"
#include "../external/grgen/classic.hpp"
#include "../external/grgen/common.h"
#include "../external/grgen/grn.hpp"
#include "config.hpp"
#include <sstream>

DNA::DNA() {
  std::vector<std::string> inputs = {"x","y","z","nt","threshold","weight","comm","div","reward"};
  for (unsigned int i = 0; i < Config::N_M; i++) {
    inputs.push_back("m"+std::to_string(i));
  }

  std::vector<std::string> outputs = {"nt","nt_t","f","f_t","comm","axon_div","axon_die","axon_none"};
  for (unsigned int i = 0; i < Config::N_M; i++) {
    outputs.push_back("m"+std::to_string(i));
  }
  outputs.push_back("m_thresh");

  for (auto& i : inputs) grn.addRandomProtein(ProteinType::input, i);

  for (auto& i : outputs) grn.addRandomProtein(ProteinType::output, i);

  grn.randomReguls(Config::GRN_REGULS);
  grn.randomParams();
}

DNA::DNA(const GRN_T &g) : grn(g) {
  reset();
}

DNA::DNA(const std::string &s) : grn(s) {
  reset();
}

DNA DNA::random() {
  DNA dna;
  return dna;
}

DNA DNA::crossover(const DNA &other) {
  GRN_T g = grn.crossover(other.grn);
  DNA res(g);
  return res;
}

void DNA::update() {
  grn.step(Config::GRN_EVO_STEPS);
}

void DNA::mutate() {
  grn.mutate();
}

void DNA::reset() {
  grn.reset();
}

void DNA::setInput(const std::string &input, double val) {
  grn.setProteinConcentration(input, ProteinType::input, val);
}

double DNA::getOutput(const std::string &output) {
  return grn.getProteinConcentration(output, ProteinType::output);
}

std::string DNA::concString() const {
  std::ostringstream stream;
  size_t i=0;
  for (auto& name : grn.getProteinNames(ProteinType::input)) {
    double conc = grn.getProteinConcentration(name, ProteinType::input);
    stream << name << ":" << std::round(conc*1000.0)/1000.0;
    i++;
    if (i < (grn.getProteinSize(ProteinType::input))) {
      stream << ", ";
    }
  }
  i=0;
  stream << " ";
  for (auto& name : grn.getProteinNames(ProteinType::output)) {
    double conc = grn.getProteinConcentration(name, ProteinType::output);
    stream << name << ":" << std::round(conc*1000.0)/1000.0;
    i++;
    if (i < (grn.getProteinSize(ProteinType::output))) {
      stream << ", ";
    }
  }
  return stream.str();
}

std::string DNA::toJSON() const {
  return grn.toJSON();
}
