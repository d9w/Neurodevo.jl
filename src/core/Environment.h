#include <map>
#include <random>
#include "Soma.h"
#include "types.hpp"
//#include "Snap.h"

using std::vector;
using std::map;

class Environment {
 public:
  using GRN = Types::DNAType;
  int age;
  vector<int> lengths;
  vector<Soma> somas;
  vector<vector<vector<vector<double> > > > morphogens;
  vector<double> max_morphogens;
  map<vector<int>, int> soma_map;
  std::mt19937 rand_engine;
  //PNEGraph ann_graph;
  //PUNGraph ann_pun_graph;

  Environment();
  Environment(vector<int> lengths, GRN grn, int seed);

  void set_random_connectivity();
  void set_morphogens();
  vector<int> move_position(vector<int> position, int morph);
  Soma* soma_at(vector<int> position);
  void develop_grns(const double reward);
  void set_nt_concentration(const vector<vector<double> > inputs);
  void set_outputs(vector<vector<double>> *outputs);
  void fire_ann();
  void axon_actions();
  void set_weights();
  void step(const vector<vector<double> > inputs, double reward);
  /*
  void populate_graph();
  double input_output_diameter();
  double layer_fit(double layer_length);
  double scale_fit(double fit);
  double modularity_fit();
  double clustering_fit();
  bool check_for_edge(vector<int> pos1, vector<int> pos2);
  double symmetry_fit();
  */
};
