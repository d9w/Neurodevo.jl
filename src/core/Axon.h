#ifndef AXON_H
#define AXON_H
#include "types.hpp"

class Axon {
 public:
  using GRN = Types::DNAType;
  GRN grn;
  //Environment env;
  vector<int> position;
  double division_conc;
  bool marked_for_deletion;
  bool marked_for_branch;
  double weight;
  int age; //action counter

  Axon();
  Axon(GRN grn, vector<int> position);

  double fire();
  void evolve(vector<double> morphogens, double nt_concentration, double soma_concentration, double soma_threshold, double reward);
  int act();
  friend std::ostream& operator<<(std::ostream& out, const Axon& a);
};
#endif
