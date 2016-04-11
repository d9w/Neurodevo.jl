#include <iostream>
#include <cmath>
#include <random>
#include "config.hpp"
#include "Environment.h"

using namespace std;

Environment::Environment() {}

Environment::Environment(vector<int> inp_lengths, GRN grn, int seed) {

  age = 0;

  rand_engine = std::mt19937(seed);

  for (unsigned int i=0; i<Config::N_D; i++) {
    lengths.push_back(inp_lengths[i]);
  }

  for (int x_i=0; x_i < lengths[0]; x_i++) {
    vector<vector<vector<double> > > ys;
    for (int y_i=0; y_i < lengths[1]; y_i++) {
      vector<vector<double> > zs;
      for (int z_i=0; z_i < lengths[2]; z_i++) {
        vector<double> concs;
        for (unsigned int c_i=0; c_i < Config::N_M; c_i++) {
          concs.push_back(0.0);
        }
        zs.push_back(concs);
      }
      ys.push_back(zs);
    }
    morphogens.push_back(ys);
  }

  int id =0;
  for (int x=0; x < lengths[0]; x++) {
    for (int y=0; y < lengths[1]; y++) {
      for (int z=0; z < lengths[2]; z++) {
        //if ((z == 4 || (x == 0 || x == lengths[0]-1 || y == 0 || y == lengths[1]-1))
        //    || z == 0 || z == lengths[2]-1) {
          //if (z == 0 || z == lengths[2]-1 || z == 4) {
        id = somas.size();
        soma_map[{x, y, z}] = id;
        somas.push_back(Soma(grn, {x, y, z}, id));
          //}
      }
    }
  }

  //ann_graph = TNEGraph::New();
  //ann_pun_graph = TUNGraph::New();

  for (unsigned int m=0; m< Config::N_M; m++) {
    max_morphogens.push_back(0.0);
  }
}

void Environment::set_random_connectivity() {
  // add axons
  for (auto& soma : somas) {
    if (soma.position[2] != lengths[2]-1) {
      for (unsigned int i=1; i < Config::AXON_MAX_NUMBER; i++) {
        soma.axons.push_back(Axon(soma.axons[0].grn, soma.axons[0].position));
      }
    }
  }

  std::uniform_int_distribution<int> xdist(0, Config::X_SIZE-1);
  std::uniform_int_distribution<int> ydist(0, Config::Y_SIZE-1);
  std::uniform_int_distribution<int> zdist(0, Config::Z_SIZE-1);
  std::uniform_real_distribution<double> dist(0.0, 2.0);

  // reposition axons
  for (auto& soma : somas) {
    soma.threshold = dist(rand_engine);
    soma.weight = dist(rand_engine)-1.0;
    for (auto& axon : soma.axons) {
      axon.weight = dist(rand_engine)-1.0;
    }
    if (soma.position[2] != lengths[2]-1) {
      for (unsigned int a=0; a<soma.axons.size(); a++) {
        int x = xdist(rand_engine);
        int y = ydist(rand_engine);
        if (a < 2) {
          soma.axons[a].position = {x, y, Config::X_SIZE-1};
        } else {
          /*
          if (!(x == 0 || x == Config::X_SIZE-1)) {
            y = dist(rand_engine) > 0.5 ? 0 : Config::Y_SIZE-1;
          }
          */
          soma.axons[a].position = {x, y, zdist(rand_engine)};
        }
      }
    }
  }
}

void Environment::set_morphogens() {
  double emission = 0.0;
  double distance = 0.0;
  double morph = 0.0;
  vector<int> positions = {0, 0, 0};
  for (int x=0; x<lengths[0]; x++) {
    for (int y=0; y<lengths[1]; y++) {
      for (int z=0; z<lengths[2]; z++) {
        for (unsigned int m=0; m < Config::N_M; m++) {
          morphogens[x][y][z][m] = 0.0;
          for (auto& soma : somas) {
            // might want to have discretized distance while in a discretized space
            for (int i=0; i<3; i++) {
              positions[i] = soma.position[i];
            }
            distance = std::max(0.5, pow((x-positions[0]),2)+pow((y-positions[1]),2)+
                                            pow((z-positions[2]),2));
            emission = soma.emission(m);
            morphogens[x][y][z][m] += emission/distance;
          }
          morph = morphogens[x][y][z][m];
          if (morph > max_morphogens[m]) {
            max_morphogens[m] = morph;
          }
        }
      }
    }
  }

  // normalize morphogens based on max for each
  for (unsigned int m=0; m<Config::N_M; m++) {
    if (max_morphogens[m] > 0.0) {
      for (int x=0; x<lengths[0]; x++) {
        for (int y=0; y<lengths[1]; y++) {
          for (int z=0; z<lengths[2]; z++) {
            morphogens[x][y][z][m] /= max_morphogens[m];
          }
        }
      }
    }
  }
}

vector<int> Environment::move_position(vector<int> position, int morph) {
  vector<vector<int>> pos_mods = {{0, 0, 1}, {0, 0, -1}, {0, 1, 0}, {0, -1, 0}, {1, 0, 0}, {-1, 0, 0}};
  int x=position[0];
  int y=position[1];
  int z=position[2];
  double max_morph = morphogens[x][y][z][morph];
  vector<int> pos_mod = {0, 0, 0};
  double conc = 0.0;
  for (auto& mod : pos_mods) {
    if (((x+mod[0])>=0) && ((x+mod[0])<lengths[0])) {
      if (((y+mod[1])>=0) && ((y+mod[1])<lengths[1])) {
        if (((z+mod[2])>=0) && ((z+mod[2])<lengths[2])) {
          conc = morphogens[x+mod[0]][y+mod[1]][z+mod[2]][morph];
          if (conc > max_morph) {
            max_morph = conc;
            pos_mod[0] = mod[0];
            pos_mod[1] = mod[1];
            pos_mod[2] = mod[2];
          }
        }
      }
    }
  }
  return {x+pos_mod[0], y+pos_mod[1], z+pos_mod[2]};
}

Soma* Environment::soma_at(vector<int> position) {
  auto it = soma_map.find(position);
  if (it != soma_map.end()) {
    return &(somas[it->second]);
  }
  return NULL;
}

void Environment::develop_grns(const double reward) {
  double soma_signal;
  double soma_concentration;
  double soma_threshold;
  for (auto& soma : somas) {
    soma_signal = soma.grn.getProteinConcentration("comm", ProteinType::output);
    soma.evolve(morphogens[soma.position[0]][soma.position[1]][soma.position[2]], reward);
    for (auto& axon : soma.axons) {
      soma_concentration = 0.0;
      soma_threshold = 0.0;
      auto rec_soma = soma_at(axon.position);
      if (rec_soma != NULL && rec_soma != &soma) {
        soma_concentration = rec_soma->nt_concentration;
        soma_threshold = rec_soma->threshold;
      }
      axon.evolve(morphogens[soma.position[0]][soma.position[1]][soma.position[2]],
                  soma_concentration, soma_threshold, soma_signal, reward);
    }
  }
}

void Environment::set_nt_concentration(const vector<vector<double> > inputs) {
  for (int x=0; x < lengths[0]; x++) {
    for (int y=0; y < lengths[1]; y++) {
      soma_at({x, y, 0})->next_concentration += inputs[y][x];
    }
  }

  for (auto& soma: somas) {
    double delta_conc = 0;
    //double delta_conc = (soma.nt_concentration - vr) * (soma.nt_concentration - vt); //QIF
    //double delta_conc = 0.1*soma.nt_concentration; //LIF
    soma.nt_concentration = soma.nt_concentration + delta_conc + soma.next_concentration;
    soma.next_concentration = 0.0;
  }
}

void Environment::set_outputs(vector<vector<double>> *outputs) {
  outputs->clear();
  for (int x=0; x < lengths[0]; x++) {
    vector<double> o;
    for (int y=0; y < lengths[1]; y++) {
      if (soma_at({x, y, 0})->fired) {
        o.push_back(1.0);
      } else {
        o.push_back(0.0);
      }
    }
    outputs->push_back(o);
  }
}

void Environment::fire_ann() {
  //std::cout << "Firing ANN ";
  for (auto& soma : somas) {
    if (soma.fire()) {
      //std::cout << soma.id << " ";
      for (auto& axon : soma.axons) {
        auto rec_soma = soma_at(axon.position);
        if (rec_soma != NULL && rec_soma != &soma && rec_soma->position[2] != 0) {
            //&& (soma.position[2] != 0 || rec_soma->position[2] != lengths[2]-1)) {
          rec_soma->next_concentration += (axon.weight + soma.weight)/2.0;
        }
      }
    }
  }
  //std::cout << std::endl;
}

void Environment::axon_actions() {
  int action = 0;
  for (auto& soma: somas) {
    for (auto& axon : soma.axons) {
      if (axon.age>0) {
        action = axon.act();
        if (action < static_cast<int>(Config::N_M)) {
          // move the axon along the morphogen gradient
          auto new_pos = move_position(axon.position, action);
          //std::cout << " move (" << axon.position[0] << "," << axon.position[1] << "," << axon.position[2] << ") ";
          //std::cout << " to (" << new_pos[0] << "," << new_pos[1] << "," << new_pos[2] << ") " << std::endl;
          for (unsigned int i=0; i<Config::N_D; i++) {
            axon.position[i] = new_pos[i];
          }
        }
      }
      axon.age++;
    }
    auto size = soma.axons.size();
    for (unsigned int i=0; i<size; i++) {
      if (soma.axons[i].marked_for_branch) {
        if (soma.axons.size() < Config::AXON_MAX_NUMBER) {
          soma.axons.push_back(Axon(soma.axons[i].grn, soma.axons[i].position));
          if (soma.axons[i].division_conc >= Config::AXON_DIVISION_REDUCTION) {
            soma.axons[i].division_conc -= Config::AXON_DIVISION_REDUCTION;
          }
        }
      }
    }
    for (unsigned int i=0; i<soma.axons.size(); ) {
      if (soma.axons[i].marked_for_deletion) {
        soma.axons.erase(soma.axons.begin()+i);
      } else {
        i++;
      }
    }
  }
}

void Environment::set_weights() {
  for (auto& soma : somas) {
    double f = soma.grn.getProteinConcentration("f", ProteinType::output);
    double f_t = soma.grn.getProteinConcentration("f_t", ProteinType::output);
    if (f > 0 || f_t > 0) {
      soma.threshold += (f - f_t)/(f + f_t);
      if (soma.threshold < 0.0) soma.threshold = 0.0;
    }
    double nt = soma.grn.getProteinConcentration("nt", ProteinType::output);
    double nt_t = soma.grn.getProteinConcentration("nt_t", ProteinType::output);
    if (nt > 0 || nt_t > 0) {
      soma.weight += (nt - nt_t)/(nt + nt_t);
    }
    for (auto& axon : soma.axons) {
      nt = axon.grn.getProteinConcentration("nt", ProteinType::output);
      nt_t = axon.grn.getProteinConcentration("nt_t", ProteinType::output);
      if (nt > 0 || nt_t > 0) {
        axon.weight += (nt - nt_t)/(nt + nt_t);
      }
    }
  }
}

void Environment::step(const vector<vector<double> > inputs, double reward) {
  // evolve GRN
  /*
  int axon_count = 0;
  for (auto& soma : somas) axon_count += soma.axons.size();
  //std::cout << "============================= DEVELOPMENT ===========================================" << std::endl;
  std::cout << somas.size() << " somata and " << axon_count << " axons ";
  std::uniform_int_distribution<int> somadist(0, somas.size()-1);
  std::cout << somas[somadist(rand_engine)] << std::endl;
  */
  develop_grns(reward);

  // act on GRN decision
  //std::cout << "=============================== GRN DECISIONS =======================================" << std::endl;
  set_weights();
  set_nt_concentration(inputs);
  //std::cout << "actions: ";
  if (age % 10 == 0) axon_actions();
  //std::cout << std::endl;

  // fire resultant ANN
  //std::cout << "================================= FIRE ANN ==========================================" << std::endl;
  fire_ann();
  age += 1;
}

/*
void Environment::populate_graph() {
  ann_graph->Clr();
  ann_pun_graph->Clr();
  for (auto& soma : somas) {
    ann_graph->AddNode(soma.id);
    ann_pun_graph->AddNode(soma.id);
  }

  for (auto& soma : somas) {
    for (auto& axon : soma.axons) {
      auto rec_soma = soma_at(axon.position);
      if (rec_soma != NULL && rec_soma != &soma && rec_soma->position[2] != 0) {
        ann_graph->AddEdge(soma.id, rec_soma->id);
        ann_pun_graph->AddEdge(soma.id, rec_soma->id);
      }
    }
  }
}

double Environment::input_output_diameter() {
  int diameter = 0;
  int count = 0;
  int path = -1;
  for (auto& soma : somas) {
    if (soma.position[2] == 0) {
      for (auto& rec_soma : somas) {
        if (rec_soma.position[2] == lengths[2]-1) {
          path = TSnap::GetShortPath(ann_graph, soma.id, rec_soma.id, true);
          if (path > 0) {
            count++;
            diameter += path;
          }
        }
      }
    }
  }
  if (count == 0) {
    return 0.0;
  }
  return (double) diameter/(1.0*count);
}

double Environment::layer_fit(double layer_length) {
  double d = input_output_diameter();
  if (d == 0.0) {
    return d;
  }
  if (d < layer_length-0.1 || d > layer_length+0.1) {
    return 0.1/std::abs(d-layer_length);
  }
  return 1.0;
}

double Environment::scale_fit(double fit) {
  double d = input_output_diameter();
  if (d == 0.0) {
    fit *= 0.1;
  }
  return fit;
}

double Environment::modularity_fit() {
  TCnComV CmtyV;
  return scale_fit(TSnap::CommunityCNM(ann_pun_graph, CmtyV));
}

double Environment::clustering_fit() {
  return scale_fit(TSnap::GetClustCf(ann_graph));
}

bool Environment::check_for_edge(vector<int> pos1, vector<int> pos2) {
  int source = -1;
  int dest = -1;
  auto it = soma_map.find(pos1);
  if (it != soma_map.end()) {
    source = it->second;
  }
  it = soma_map.find(pos2);
  if (it != soma_map.end()) {
    dest = it->second;
  }
  if (source < 0 || dest < 0) {
    return false;
  }
  return ann_pun_graph->IsEdge(source, dest);
}

double Environment::symmetry_fit() {
  int n_edges = ann_pun_graph->GetEdges();
  if (n_edges == 0) {
    return 0.0;
  }
  int x_sym = 0;
  int y_sym = 0;
  for (TUNGraph::TEdgeI EI = ann_pun_graph->BegEI(); EI < ann_pun_graph->EndEI(); EI++) {
    auto source_pos = somas[EI.GetSrcNId()].position;
    auto dest_pos = somas[EI.GetDstNId()].position;
    source_pos[0] = lengths[0] - source_pos[0];
    dest_pos[0] = lengths[0] - dest_pos[0];
    if (check_for_edge(source_pos, dest_pos)) {
      x_sym++;
    }
    source_pos[0] = lengths[0] - source_pos[0];
    dest_pos[0] = lengths[0] - dest_pos[0];
    source_pos[1] = lengths[1] - source_pos[1];
    dest_pos[1] = lengths[1] - dest_pos[1];
    if (check_for_edge(source_pos, dest_pos)) {
      y_sym++;
    }
  }
  double fit = 1.0/(2.0*n_edges) * (x_sym*1.0 + y_sym*1.0);
  return scale_fit(fit);
}
*/
