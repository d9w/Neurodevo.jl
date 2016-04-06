#include <iostream>
#include "config.hpp"
#include "../external/grgen/common.h"
#include "Axon.h"

Axon::Axon(GRN inp_grn, vector<int> inp_position) {
  grn = GRN(inp_grn);
  for (int i=0; i<Config::N_D; i++) {
    position.push_back(inp_position[i]);
  }
  marked_for_deletion = false;
  marked_for_branch = false;
  division_conc = 0.2;
  age = 0;
}

double Axon::fire() {
  double nt_amount = grn.getProteinConcentration("nt", ProteinType::output);
  double nt_thresh = grn.getProteinConcentration("nt_t", ProteinType::output);

  double fire_conc = 0.0;
  if (nt_amount > 0.0 || nt_thresh > 0.0) {
    fire_conc =  (nt_amount - nt_thresh)/(nt_amount + nt_thresh);
  }
  std::cout << "Transmitting: " << fire_conc << ", " << nt_amount << ", " << nt_thresh;
  return fire_conc;
}

void Axon::evolve(vector<double> morphogens, double soma_concentration, double soma_signal) {
  /*
  for (int i=0; i<N_M; i++) {
    grn.proteins[i].concentration = morphogens[i];
  }
  */
  grn.setProteinConcentration("x", ProteinType::input, (double)position[0]/Config::X_SIZE);
  grn.setProteinConcentration("y", ProteinType::input, (double)position[1]/Config::Y_SIZE);
  grn.setProteinConcentration("z", ProteinType::input, (double)position[2]/Config::Z_SIZE);
  grn.setProteinConcentration("nt", ProteinType::input,
                              std::exp(soma_concentration/((double)Config::AXON_MAX_NUMBER+2.0)));
  grn.setProteinConcentration("comm", ProteinType::input, soma_signal);
  grn.setProteinConcentration("div", ProteinType::input, division_conc);

  grn.step(Config::GRN_EVO_STEPS);
}

int Axon::act() {
  /*
  marked_for_branch = false;
  marked_for_deletion = false;

  int action = 0;
  double max_action_conc = 0.0;
  for (int i=AXON_GRN_INPUTS; i<=AXON_GRN_OUTPUT_NOTHING; i++) {
    if (grn.proteins[i].concentration > max_action_conc) {
      action = i-AXON_GRN_INPUTS;
      max_action_conc = grn.proteins[i].concentration;
    }
  }

  switch(action) {
  case N_M: marked_for_branch=true; break;
  case N_M+1: marked_for_deletion=true; break;
  case N_M+2: break;
  default: break;
  }
  return action;
  */
  return 0;
}

std::ostream& operator<<(std::ostream& out, const Axon& a) {
  string grn_outputs = "";
  /*
  for (int c=0; c<(N_G_I+GRN_OUTPUTS); c++) {
    grn_outputs += std::to_string((int)std::round(a.grn.proteins[c].concentration*100));
    if (c<(N_G_I+GRN_OUTPUTS-1)) {
      grn_outputs += ", ";
    }
  }
  */
  return out << "Axon at (" << a.position[0] << "," << a.position[1] << "," << a.position[2]
             << ")\tage: " << a.age << "\tgrn concs: [" << grn_outputs << "]";
}
