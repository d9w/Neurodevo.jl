#include <iostream>
#include "easylogging++.h"
#include "Constants.h"
#include "GRN.h"
#include "Axon.h"

Axon::Axon(GRN inp_grn, vector<int> inp_position) {
  grn = inp_grn.copy();
  for (int i=0; i<N_D; i++) {
    position.push_back(inp_position[i]);
  }
  marked_for_deletion = false;
  marked_for_branch = false;
  division_conc = 0.2;
  age = 0;
}

double Axon::fire() {
  double nt_amount = grn.proteins[GRN_OUTPUT_NT].concentration;
  double nt_thresh = grn.proteins[GRN_OUTPUT_NT_T].concentration;

  double fire_conc = 0.0;
  if (nt_amount > 0.0 || nt_thresh > 0.0) {
    fire_conc =  (nt_amount - nt_thresh)/(nt_amount + nt_thresh);
  }
  //LOG(DEBUG) << "Transmitting: " << fire_conc << ", " << nt_amount << ", " << nt_thresh;
  return fire_conc;
}

void Axon::evolve(vector<double> morphogens, double soma_concentration, double soma_signal) {
  for (int i=0; i<N_M; i++) {
    grn.proteins[i].concentration = morphogens[i];
  }
  grn.proteins[GRN_INPUT_X].concentration = (double)position[0]/X_SIZE;
  grn.proteins[GRN_INPUT_Y].concentration = (double)position[1]/Y_SIZE;
  grn.proteins[GRN_INPUT_Z].concentration = (double)position[2]/Z_SIZE;
  grn.proteins[GRN_INPUT_NEUROTRANSMITTER].concentration =
    std::exp(soma_concentration/((double)AXON_MAX_NUMBER+2.0));
  grn.proteins[GRN_INPUT_COMM].concentration = soma_signal;
  grn.proteins[GRN_INPUT_DIVISION].concentration = division_conc;

  grn.evolve(GRN_EVO_STEPS);
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
  for (int c=0; c<(N_G_I+GRN_OUTPUTS); c++) {
    grn_outputs += std::to_string((int)std::round(a.grn.proteins[c].concentration*100));
    if (c<(N_G_I+GRN_OUTPUTS-1)) {
      grn_outputs += ", ";
    }
  }
  return out << "Axon at (" << a.position[0] << "," << a.position[1] << "," << a.position[2]
             << ")\tage: " << a.age << "\tgrn concs: [" << grn_outputs << "]";
}
