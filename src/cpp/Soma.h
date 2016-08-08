#include "GRN.h"
#include "Axon.h"

class Soma {
 public:
  GRN grn;
  GRN axon_grn;
  //  Environment env;
  vector<int> position;
  vector<Axon> axons;
  double nt_concentration;
  double next_concentration;
  bool fired;
  int id;

  Soma();
  Soma(GRN grn, GRN axon_grn, vector<int> position, int id);

  double emission(int);
  void evolve(vector<double> morphogens, double reward);
  bool fire();
  friend std::ostream& operator<<(std::ostream& out, const Soma& s);
};
