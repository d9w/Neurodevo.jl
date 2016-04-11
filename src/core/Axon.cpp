#include <iostream>
#include "config.hpp"
#include "../external/grgen/common.h"
#include "Axon.h"

Axon::Axon(GRN inp_grn, vector<int> inp_position) {
  grn = GRN(inp_grn);
  for (unsigned int i=0; i<Config::N_D; i++) {
    position.push_back(inp_position[i]);
  }
  marked_for_deletion = false;
  marked_for_branch = false;
  division_conc = 0.2;
  age = 0;
  weight = 0.0;
}

double Axon::fire() {
  /*
  double nt_amount = grn.getProteinConcentration("nt", ProteinType::output);
  double nt_thresh = grn.getProteinConcentration("nt_t", ProteinType::output);

  double fire_conc = 0.0;
  if (nt_amount > 0.0 || nt_thresh > 0.0) {
    fire_conc =  (nt_amount - nt_thresh)/(nt_amount + nt_thresh);
  }
  */
  double fire_conc = weight;
  //std::cout << "Transmitting: " << fire_conc << ", " << nt_amount << ", " << nt_thresh;
  return fire_conc;
}

void Axon::evolve(vector<double> morphogens, double soma_concentration, double soma_threshold, double soma_signal, double reward) {
  for (unsigned int i=0; i<Config::N_M; i++) {
    grn.setProteinConcentration("m"+std::to_string(i), ProteinType::input, morphogens[i]);
  }
  grn.setProteinConcentration("x", ProteinType::input, (double)position[0]/Config::X_SIZE);
  grn.setProteinConcentration("y", ProteinType::input, (double)position[1]/Config::Y_SIZE);
  grn.setProteinConcentration("z", ProteinType::input, (double)position[2]/Config::Z_SIZE);
  grn.setProteinConcentration("nt", ProteinType::input,
                              std::exp(soma_concentration/((double)Config::AXON_MAX_NUMBER+2.0)));
  grn.setProteinConcentration("weight", ProteinType::input,
                              std::exp(weight));
  grn.setProteinConcentration("threshold", ProteinType::input,
                              std::exp(soma_threshold));
  grn.setProteinConcentration("comm", ProteinType::input, soma_signal);
  grn.setProteinConcentration("div", ProteinType::input, division_conc);
  grn.setProteinConcentration("reward", ProteinType::input, reward);

  grn.step(Config::GRN_EVO_STEPS);
}

int Axon::act() {
  marked_for_branch = false;
  marked_for_deletion = false;

  int action = 0;
  double max_action_conc = 0.0;
  vector<string> outputs = {"axon_div", "axon_die", "axon_none"};
  for (unsigned int i=0; i<outputs.size()+Config::N_M; i++) {
    string out = "";
    if (i<Config::N_M) {
      out = "m" + std::to_string(i);
    } else {
      out = outputs[i-Config::N_M];
    }
    double conc = grn.getProteinConcentration(out, ProteinType::output);
    if (conc >= max_action_conc) {
      action = i;
      max_action_conc = conc;
    }
  }

  switch(action) {
  case Config::N_M: marked_for_branch=true; break;
  case Config::N_M+1: marked_for_deletion=true; break;
  case Config::N_M+2: break;
  default: break;
  }
  return action;
}

std::ostream& operator<<(std::ostream& out, const Axon& a) {
  string grnstring = "";
  return out << "Axon at (" << a.position[0] << "," << a.position[1] << "," << a.position[2]
             << ")\tage: " << a.age << "weight:" << a.weight << "\tgrn: [" << grnstring << "]";
}
