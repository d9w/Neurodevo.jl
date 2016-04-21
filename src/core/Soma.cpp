#include <iostream>
#include "Soma.h"
#include "config.hpp"

Soma::Soma(const DNA& inp_dna, vector<int> inp_position, int sid) {
  id = sid;
  dna = DNA(inp_dna);
  for (unsigned int i=0; i<Config::N_D; i++) {
    position.push_back(inp_position[i]);
  }

  nt_concentration = 0.0;
  next_concentration = 0.0;
  threshold = 0.0;
  weight = 0.0;
  fired = false;

  if (position[2] < static_cast<int>(Config::Z_SIZE)-1) { // don't give output neurons an axon for now
    axons.push_back(Axon(dna, position));
  }
}

double Soma::emission(int m) {
  double morph = dna.getOutput("m"+std::to_string(m));
  double thresh = dna.getOutput("m_thresh");
  if (morph == 0) {
    return morph;
  }
  return morph / (thresh + morph);
}

void Soma::evolve(vector<double> morphogens, double reward) {
  for (unsigned int i=0; i<Config::N_M; i++) {
    dna.setInput("m"+std::to_string(i), morphogens[i]);
  }
  dna.setInput("x", (double)position[0]/Config::X_SIZE);
  dna.setInput("y", (double)position[1]/Config::Y_SIZE);
  dna.setInput("z", (double)position[2]/Config::Z_SIZE);
  dna.setInput("nt", std::exp(nt_concentration/((double)Config::AXON_MAX_NUMBER+2.0)));
  dna.setInput("weight", std::exp(weight));
  dna.setInput("threshold", std::exp(threshold));
  double axon_input = 0.0;
  for (auto axon: axons) {
    axon_input += axon.dna.getOutput("comm");
  }
  if (axon_input > Config::SOMA_AXON_INPUT_THRESH) {
    axon_input = Config::SOMA_AXON_INPUT_THRESH;
  }
  dna.setInput("comm", axon_input);
  dna.setInput("div", 1.0);
  dna.setInput("reward", reward);

  dna.update();
}

bool Soma::fire() {
  /*
  double p1 = grn.getProteinConcentration("f", ProteinType::output);
  double p2 = grn.getProteinConcentration("f_t", ProteinType::output);
  double vt = 1.0 + (double)(p1-p2)/(p1+p2);
  */
  if (!(position[2] == 0) && nt_concentration > 0.0)
    std::cout << "[" << id << "," << nt_concentration << "," << threshold << "] ";
  double vt = threshold;
  double vr = 0.0;
  fired = false;
  if (nt_concentration > vt) {
    //if (!(position[2] == 0)) std::cout << "[" << id << "," << nt_concentration << "," << threshold << "] ";
    nt_concentration = vr;
    fired = true;
  }
  return fired;
}

std::ostream& operator<<(std::ostream& out, const Soma& s) {
  auto dnastring = s.dna.concString();
  return out << "Soma " << s.id << " at (" << s.position[0] << "," << s.position[1] << ","
             << s.position[2] << ")\tnt: (" << s.nt_concentration
             << "," << s.threshold << "," << s.weight << ")\tdna: [" << dnastring << "]";
}
