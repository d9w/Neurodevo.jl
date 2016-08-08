#include <iostream>
#include <cstdio>
#include <time.h>
#include "easylogging++.h"
#include "color.h"
#include "GRN.h"
#include "Environment.h"
#include "Snap.h"

INITIALIZE_EASYLOGGINGPP

int main(int argc, char** argv)
{
  srand(time(NULL));

  el::Configurations conf("log.conf");
  el::Loggers::reconfigureAllLoggers(conf);

  LOG(DEBUG) << BOLD << "============================================" << RESET;
  LOG(DEBUG) << BOLD << "==========TESTING THE ENVIRONMENT===========" << RESET;
  LOG(DEBUG) << BOLD << "============================================" << RESET;
  int x = 10;
  int y = 10;
  int z = 10;
  int num_axons = 0;
  vector<int> lengths = {x, y, z};
  GRN soma_grn = GRN(argv[1]);
  GRN axon_grn = GRN(argv[2]);
  Environment env (lengths, soma_grn, axon_grn);
  vector<vector<double> > inputs;
  int t_action = (int)std::round((soma_grn.t_action+axon_grn.t_action)/2.0);

  for (int step=0; step<GA_EVAL_STEPS; step++) {
    inputs.clear();
    for (int i=0; i<lengths[0]; i++) {
      vector<double> inp;
      for (int j=0; j<lengths[1]; j++) {
        inp.push_back(0.5*sin(i*step+j)+0.5);
        //inp.push_back((double)rand()/RAND_MAX);
      }
      inputs.push_back(inp);
    }
    env.develop_grns();

    env.set_nt_concentration(inputs);
    env.set_morphogens();

    num_axons = 0;
    if (step % t_action == 0) {
      env.axon_actions();
      for (auto& soma : env.somas) {
        for (auto& axon : soma.axons) {
          num_axons++;
        }
      }
      LOG(DEBUG) << BOLD << num_axons << " axons" << RESET;
    }

    env.fire_ann();
    if (step > (GA_EVAL_STEPS-t_action)) {
      int num_fired = 0;
      for(auto& soma : env.somas) {
        if (soma.position[2] != 0 && soma.fired) {
          num_fired++;
        }
      }
      LOG(DEBUG) << GREEN << num_fired << " fired " << RESET;
    }
  }

  LOG(DEBUG) << MAGENTA << "==============random neuron sample===================" << RESET;
  for (auto& soma : env.somas) {
    if ((double)rand()/RAND_MAX < 0.05) {
      std::cout << soma << std::endl;
      for (auto& axon : soma.axons) {
        std::cout << "\t" << CYAN << axon << RESET << std::endl;
      }
    }
  }
  env.populate_graph();
  LOG(DEBUG) << MAGENTA << "===================fitnesses========================" << RESET;
  LOG(DEBUG) << BOLD << "t_action: " << t_action << RESET;
  LOG(DEBUG) << BOLD << "avg layer number: " << env.input_output_diameter() << RESET;
  LOG(DEBUG) << BOLD << "shallow fit: " << env.layer_fit(3.5) << RESET;
  LOG(DEBUG) << BOLD << "deep fit: " << env.layer_fit(10.5) << RESET;
  LOG(DEBUG) << BOLD << "modularity fit: " << env.modularity_fit() << RESET;
  LOG(DEBUG) << BOLD << "symmetry fit: " << env.symmetry_fit() << RESET;
  TVec<TPair<TInt, TInt> > CntV; // vector of pairs of integers (size, count)
  TSnap::GetWccSzCnt(env.ann_graph, CntV);
  // get degree distribution pairs (degree, count)
  TSnap::GetOutDegCnt(env.ann_graph, CntV);
  // get first eigenvector of graph adjacency matrix
  TFltV EigV; // vector of floats
  TSnap::GetEigVec(env.ann_pun_graph, EigV);
  // get diameter of G
  LOG(DEBUG) << BOLD << "diameter: " << TSnap::GetBfsFullDiam(env.ann_graph, 5, true) << RESET;
  // count the number of triads in G, get the clustering coefficient of G
  LOG(DEBUG) << BOLD << "triads: " << TSnap::GetTriads(env.ann_graph) << RESET;
  LOG(DEBUG) << BOLD <<  "clustering coef:" << env.clustering_fit() << RESET;
  return 0;
}
