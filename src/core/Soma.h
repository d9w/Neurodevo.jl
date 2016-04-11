#ifndef SOMA_H
#define SOMA_H
#include "Axon.h"
#include "types.hpp"

class Soma {
 public:
  using GRN = Types::DNAType;
  GRN grn;
  //  Environment env;
  vector<int> position;
  vector<Axon> axons;
  double nt_concentration;
  double next_concentration;
  double threshold;
  double weight;
  bool fired;
  int id;

  Soma();
  Soma(GRN grn, vector<int> position, int id);

  double emission(int);
  void evolve(vector<double> morphogens, double reward);
  bool fire();
  friend std::ostream& operator<<(std::ostream& out, const Soma& s);
};
#endif
