#include <iostream>
#include "easylogging++.h"
#include "Constants.h"
#include "GRN.h"
#include "Soma.h"

Soma::Soma(GRN inp_grn, GRN inp_axon_grn, vector<int> inp_position, int sid) {
  id = sid;
  grn = inp_grn.copy();
  axon_grn = inp_axon_grn.copy();
  for (int i=0; i<N_D; i++) {
    position.push_back(inp_position[i]);
  }
  //position.swap(inp_position);

  nt_concentration = 0.0;
  next_concentration = 0.0;
  fired = false;

  if (position[2] < X_SIZE-1) { // don't give output neurons an axon for now
    axons.push_back(Axon(axon_grn, position));
  }
}

double Soma::emission(int m) {
  /*
  double morph = grn.proteins[SOMA_GRN_INPUTS+m].concentration;
  double thresh = grn.proteins[SOMA_GRN_OUTPUT_MORPHOGEN_THRESH].concentration;
  if (morph == 0) {
    return morph;
  }
  return morph / (thresh + morph);
  */
  return 0;
}

void Soma::evolve(vector<double> morphogens, double reward) {
  // decide GRN action output and call it
  for (int i=0; i<N_M; i++) {
    grn.proteins[i].concentration = morphogens[i];
  }
  grn.proteins[GRN_INPUT_X].concentration = (double)position[0]/X_SIZE;
  grn.proteins[GRN_INPUT_Y].concentration = (double)position[1]/Y_SIZE;
  grn.proteins[GRN_INPUT_Z].concentration = (double)position[2]/Z_SIZE;
  grn.proteins[GRN_INPUT_NEUROTRANSMITTER].concentration =
    std::exp(nt_concentration/((double)AXON_MAX_NUMBER+2.0));
  double axon_input = 0.0;
  for (auto axon: axons) {
    axon_input += axon.grn.proteins[GRN_OUTPUT_COMM].concentration;
  }
  if (axon_input > SOMA_AXON_INPUT_THRESH) {
    axon_input = SOMA_AXON_INPUT_THRESH;
  }
  grn.proteins[GRN_INPUT_COMM].concentration = axon_input;
  grn.proteins[GRN_INPUT_DIVISION].concentration = 1.0;
  grn.proteins[GRN_INPUT_REWARD].concentration = reward;

  grn.evolve(GRN_EVO_STEPS);
}

bool Soma::fire() {
  double p1 = grn.proteins[GRN_OUTPUT_F_NULL].concentration;
  double p2 = grn.proteins[GRN_OUTPUT_FT_NULL].concentration;
  double vt = 1.0 + (double)(p1-p2)/(p1+p2);
  double vr = 0.0;
  fired = false;
  if (nt_concentration > vt) {
    nt_concentration = vr;
    fired = true;
  }
  return fired;
}

std::ostream& operator<<(std::ostream& out, const Soma& s) {
  string grn_outputs = "";
  for (int c=0; c<(N_G_I+GRN_OUTPUTS); c++) {
    grn_outputs += std::to_string((int)std::round(s.grn.proteins[c].concentration*100));
    if (c<(N_G_I+GRN_OUTPUTS-1)) {
      grn_outputs += ", ";
    }
  }
  return out << "Soma " << s.id << " at (" << s.position[0] << "," << s.position[1] << ","
             << s.position[2] << ")\tnt: (" << s.nt_concentration
             << "," << s.next_concentration << ")\tgrn concs: [" << grn_outputs << "]";
}
