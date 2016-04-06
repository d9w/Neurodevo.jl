#include <iostream>
#include "Soma.h"
#include "../external/grgen/common.h"
#include "config.hpp"

Soma::Soma(GRN inp_grn, vector<int> inp_position, int sid) {
  id = sid;
  grn = GRN(inp_grn);
  for (unsigned int i=0; i<Config::N_D; i++) {
    position.push_back(inp_position[i]);
  }

  nt_concentration = 0.0;
  next_concentration = 0.0;
  fired = false;

  if (position[2] < static_cast<int>(Config::Z_SIZE)-1) { // don't give output neurons an axon for now
    axons.push_back(Axon(grn, position));
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
  /*
  for (int i=0; i<N_M; i++) {
    grn.proteins[i].concentration = morphogens[i];
  }
  */
  grn.setProteinConcentration("x", ProteinType::input, (double)position[0]/Config::X_SIZE);
  grn.setProteinConcentration("y", ProteinType::input, (double)position[1]/Config::Y_SIZE);
  grn.setProteinConcentration("z", ProteinType::input, (double)position[2]/Config::Z_SIZE);
  grn.setProteinConcentration("nt", ProteinType::input,
                              std::exp(nt_concentration/((double)Config::AXON_MAX_NUMBER+2.0)));
  double axon_input = 0.0;
  for (auto axon: axons) {
    axon_input += axon.grn.getProteinConcentration("comm", ProteinType::output);
  }
  if (axon_input > Config::SOMA_AXON_INPUT_THRESH) {
    axon_input = Config::SOMA_AXON_INPUT_THRESH;
  }
  grn.setProteinConcentration("comm", ProteinType::input, axon_input);
  grn.setProteinConcentration("div", ProteinType::input, 1.0);
  grn.setProteinConcentration("reward", ProteinType::input, reward);

  grn.step(Config::GRN_EVO_STEPS);
}

bool Soma::fire() {
  double p1 = grn.getProteinConcentration("f", ProteinType::output);
  double p2 = grn.getProteinConcentration("f_t", ProteinType::output);
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
  /*
  for (int c=0; c<(N_G_I+GRN_OUTPUTS); c++) {
    grn_outputs += std::to_string((int)std::round(s.grn.proteins[c].concentration*100));
    if (c<(N_G_I+GRN_OUTPUTS-1)) {
      grn_outputs += ", ";
    }
  }
  */
  return out << "Soma " << s.id << " at (" << s.position[0] << "," << s.position[1] << ","
             << s.position[2] << ")\tnt: (" << s.nt_concentration
             << "," << s.next_concentration << ")\tgrn concs: [" << grn_outputs << "]";
}
