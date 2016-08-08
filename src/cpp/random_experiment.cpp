#include <iostream>
#include <cstdio>
#include <time.h>
#include <iostream>
#include <sstream>
#include <iterator>
#include <numeric>
#include "easylogging++.h"
#include "color.h"
#include "GRN.h"
#include "Environment.h"
#include "Snap.h"

INITIALIZE_EASYLOGGINGPP

int main(int argc, char** argv)
{
  srand(time(NULL));

  el::Configurations conf("random_log.conf");
  el::Loggers::reconfigureAllLoggers(conf);
  el::Logger* defaultLogger = el::Loggers::getLogger("default");


  int x = 10;
  int y = 10;
  int z = 10;
  int num_axons = 0;
  vector<int> lengths = {x, y, z};

  LOG(INFO) << "iter\taxon\tsoma\tlayers\tshallow\tdeep\tmodular\tsym\tdiam\ttriads\tcluster";

  int max_n = atoi(argv[1]);

  for (int n=0; n<max_n; n++) {
    int soma_reg_length = (int) std::round(((double)rand()/RAND_MAX)*
                                           (SOMA_GRN_REGULS+MAX_ITER*ADD_MUTATION_RATE));
    int axon_reg_length = (int) std::round(((double)rand()/RAND_MAX)*
                                           (AXON_GRN_REGULS+MAX_ITER*ADD_MUTATION_RATE));
    GRN soma_grn = GRN(AXON_GRN_INPUTS, AXON_GRN_OUTPUTS, axon_reg_length);
    GRN axon_grn = GRN(SOMA_GRN_INPUTS, SOMA_GRN_OUTPUTS, soma_reg_length);
    Environment env (lengths, soma_grn, axon_grn);
    //env.set_random_connectivity();
    //env.develop_grns();
    //env.populate_graph();
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

    LOG(DEBUG) << MAGENTA << "==========random neuron sample=============" << RESET;
    for (auto& soma : env.somas) {
      if ((double)rand()/RAND_MAX < 0.01) {
        std::cout << soma << std::endl;
        for (auto& axon : soma.axons) {
          std::cout << "\t" << CYAN << axon << RESET << std::endl;
        }
      }
    }
    env.populate_graph();
    LOG(DEBUG) << MAGENTA << "===============fitnesses===================" << RESET;
    double z_avg_diameter = env.input_output_diameter();
    double shallow_fit = env.layer_fit(3.5);
    double deep_fit = env.layer_fit(10.5);
    double modularity_fit = env.modularity_fit();
    double symmetry_fit = env.symmetry_fit();
    double diameter = TSnap::GetBfsFullDiam(env.ann_graph, 2, true);
    double triads = TSnap::GetTriads(env.ann_graph);
    double clustering_coeff = TSnap::GetClustCf(env.ann_graph);
    int grn_length = soma_reg_length + axon_reg_length;

    const char *fmt = "%03d\t%03d\t%03d\t%05.1f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%05.1f\t%05.1f\t%0.4f\t";
    int sz = std::snprintf(nullptr, 0, fmt, n, axon_reg_length, soma_reg_length, z_avg_diameter,
                           shallow_fit, deep_fit, modularity_fit, symmetry_fit, diameter, triads,
                           clustering_coeff);
    std::vector<char> buf(sz + 1); // note +1 for null terminator
    std::snprintf(&buf[0], buf.size(), fmt, n, axon_reg_length, soma_reg_length, z_avg_diameter,
                  shallow_fit, deep_fit, modularity_fit, symmetry_fit, diameter, triads,
                  clustering_coeff);
    string out = "";
    for (auto& c : buf) {
      out += c;
    }
    LOG(DEBUG) << BOLD << "t_action: " << t_action << RESET;
    LOG(DEBUG) << BOLD << "avg layer number: " << z_avg_diameter << RESET;
    LOG(DEBUG) << BOLD << "shallow fit: " << shallow_fit << RESET;
    LOG(DEBUG) << BOLD << "deep fit: " << deep_fit << RESET;
    LOG(DEBUG) << BOLD << "modularity fit: " << modularity_fit << RESET;
    LOG(DEBUG) << BOLD << "symmetry fit: " << symmetry_fit << RESET;
    LOG(DEBUG) << BOLD << "diameter: " << diameter << RESET;
    LOG(DEBUG) << BOLD << "triads: " << triads << RESET;
    LOG(DEBUG) << BOLD <<  "clustering coeff:" << clustering_coeff << RESET;
    LOG(INFO) << out;
  }
  return 0;
}
