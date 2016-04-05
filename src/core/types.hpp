#ifndef TYPES_HPP
#define TYPES_HPP
#include "../external/grgen/common.h"
#include "../external/grgen/grn.hpp"
#include "../external/grgen/classic.hpp"
/**
#include "dna.hpp"
#include "neuron.hpp"
#include "environment.hpp"
*/

struct Types {
  using DNAType = GRN<Classic>;
  //using NeuronType = Neuron<GRNType>;
  //using EnvironmentType = Environment<NeuronType>;
};
#endif
