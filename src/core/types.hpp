#ifndef TYPES_HPP
#define TYPES_HPP
#include "../external/grgen/common.h"
#include "../external/grgen/grn.hpp"
#include "../external/grgen/classic.hpp"
#include "config.hpp"
/**
#include "dna.hpp"
#include "neuron.hpp"
#include "environment.hpp"
*/

struct Types {

  Types() {}

  using DNAType = GRN<Classic>;
  //using NeuronType = Neuron<GRNType>;
  //using EnvironmentType = Environment<NeuronType>;

  DNAType random_dna() {
    DNAType t;
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

    t.randomReguls(Config::GRN_REGULS);
    t.randomParams();
    return t;
  }
};
#endif
