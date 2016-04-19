#ifndef SOMA_H
#define SOMA_H
#include "Axon.h"
#include "DNA.h"

class Soma {
 public:
  DNA dna;
  vector<int> position;
  vector<Axon> axons;
  double nt_concentration;
  double next_concentration;
  double threshold;
  double weight;
  bool fired;
  int id;

  Soma();
  Soma(DNA dna, vector<int> position, int id);

  double emission(int);
  void evolve(vector<double> morphogens, double reward);
  bool fire();
  friend std::ostream& operator<<(std::ostream& out, const Soma& s);
};
#endif
