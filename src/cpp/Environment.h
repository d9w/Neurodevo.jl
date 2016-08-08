#include <map>
#include "GRN.h"
#include "Soma.h"
//#include "Snap.h"

class Environment {
 public:
  vector<int> lengths;
  vector<Soma> somas;
  vector<vector<vector<vector<double> > > > morphogens;
  vector<double> max_morphogens;
  map<vector<int>, int> soma_map;
  //PNEGraph ann_graph;
  //PUNGraph ann_pun_graph;

  Environment();
  Environment(vector<int> lengths, GRN soma_grn, GRN axon_grn);

  void set_random_connectivity(int seed);
  void set_morphogens();
  vector<int> move_position(vector<int> position, int morph);
  Soma* soma_at(vector<int> position);
  void develop_grns(double reward);
  void set_nt_concentration(vector<vector<double> > inputs);
  void fire_ann();
  void axon_actions();
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
