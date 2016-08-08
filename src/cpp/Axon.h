#include "GRN.h"

class Axon {
 public:
  GRN grn;
  //Environment env;
  vector<int> position;
  double division_conc;
  bool marked_for_deletion;
  bool marked_for_branch;
  int age; //action counter

  Axon();
  Axon(GRN grn, vector<int> position);

  double fire();
  void evolve(vector<double> morphogens, double nt_concentration, double soma_concentration);
  int act();
  friend std::ostream& operator<<(std::ostream& out, const Axon& a);
};
