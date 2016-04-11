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
    for (unsigned int i = 0; i < Config::N_M; i++) {
      t.addRandomProtein(ProteinType::input, "m" + std::to_string(i));
    }
    t.addRandomProtein(ProteinType::input, "x");
    t.addRandomProtein(ProteinType::input, "y");
    t.addRandomProtein(ProteinType::input, "z");
    t.addRandomProtein(ProteinType::input, "nt");
    t.addRandomProtein(ProteinType::input, "threshold");
    t.addRandomProtein(ProteinType::input, "weight");
    t.addRandomProtein(ProteinType::input, "comm");
    t.addRandomProtein(ProteinType::input, "div");
    t.addRandomProtein(ProteinType::input, "reward");

    for (unsigned int i = 0; i < Config::N_M; i++) {
      t.addRandomProtein(ProteinType::output, "m" + std::to_string(i));
    }
    t.addRandomProtein(ProteinType::output, "nt");
    t.addRandomProtein(ProteinType::output, "nt_t");
    t.addRandomProtein(ProteinType::output, "f");
    t.addRandomProtein(ProteinType::output, "f_t");
    t.addRandomProtein(ProteinType::output, "comm");
    t.addRandomProtein(ProteinType::output, "axon_div");
    t.addRandomProtein(ProteinType::output, "axon_die");
    t.addRandomProtein(ProteinType::output, "axon_none");

    //t.randomReguls(Config::GRN_REGULS);
    t.randomReguls(10);
    t.randomParams();
    return t;
  }
};
#endif
