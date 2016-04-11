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
  threshold = 0.0;
  weight = 0.0;
  fired = false;

  if (position[2] < static_cast<int>(Config::Z_SIZE)-1) { // don't give output neurons an axon for now
    axons.push_back(Axon(grn, position));
  }

}

double Soma::emission(int m) {
  double morph = grn.getProteinConcentration("m"+std::to_string(m), ProteinType::output);
  double thresh = grn.getProteinConcentration("m_thresh", ProteinType::output);
  if (morph == 0) {
    return morph;
  }
  return morph / (thresh + morph);
}

void Soma::evolve(vector<double> morphogens, double reward) {
  for (unsigned int i=0; i<Config::N_M; i++) {
    grn.setProteinConcentration("m"+std::to_string(i), ProteinType::input, morphogens[i]);
  }
  grn.setProteinConcentration("x", ProteinType::input, (double)position[0]/Config::X_SIZE);
  grn.setProteinConcentration("y", ProteinType::input, (double)position[1]/Config::Y_SIZE);
  grn.setProteinConcentration("z", ProteinType::input, (double)position[2]/Config::Z_SIZE);
  grn.setProteinConcentration("nt", ProteinType::input,
                              std::exp(nt_concentration/((double)Config::AXON_MAX_NUMBER+2.0)));
  grn.setProteinConcentration("weight", ProteinType::input,
                              std::exp(weight));
  grn.setProteinConcentration("threshold", ProteinType::input,
                              std::exp(threshold));
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
  /*
  double p1 = grn.getProteinConcentration("f", ProteinType::output);
  double p2 = grn.getProteinConcentration("f_t", ProteinType::output);
  double vt = 1.0 + (double)(p1-p2)/(p1+p2);
  */
  double vt = threshold;
  double vr = 0.0;
  fired = false;
  if (nt_concentration > vt) {
    nt_concentration = vr;
    fired = true;
  }
  return fired;
}

std::ostream& operator<<(std::ostream& out, const Soma& s) {
  string grnstring = "";
  return out << "Soma " << s.id << " at (" << s.position[0] << "," << s.position[1] << ","
             << s.position[2] << ")\tnt: (" << s.nt_concentration
             << "," << s.threshold << "," << s.weight << ")\tgrn: [" << grnstring << "]";
}
