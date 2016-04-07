#include <iostream>
#include <cmath>
#include "config.hpp"
#include "Environment.h"

using namespace std;

Environment::Environment() {}

Environment::Environment(vector<int> inp_lengths, GRN grn) {

  for (unsigned int i=0; i<Config::N_D; i++) {
    lengths.push_back(inp_lengths[i]);
  }
  //lengths.swap(inp_lengths);

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
        if ((z == 4 && (x == 0 || x == lengths[0]-1 || y == 0 || y == lengths[1]-1))
            || z == 0 || z == lengths[2]-1) {
          //if (z == 0 || z == lengths[2]-1 || z == 4) {
          id = somas.size();
          soma_map[{x, y, z}] = id;
          somas.push_back(Soma(grn, {x, y, z}, id));
        }
      }
    }
  }

  //ann_graph = TNEGraph::New();
  //ann_pun_graph = TUNGraph::New();

  for (unsigned int m=0; m< Config::N_M; m++) {
    max_morphogens.push_back(0.0);
  }
}

void Environment::set_random_connectivity(int seed) {
  srand(123456+100*seed);
  // add axons
  for (auto& soma : somas) {
    if (soma.position[2] != lengths[2]-1) {
      for (unsigned int i=1; i < Config::AXON_MAX_NUMBER; i++) {
        soma.axons.push_back(Axon(soma.axons[0].grn, soma.axons[0].position));
      }
    }
  }

  // reposition axons
  for (auto& soma : somas) {
    if (soma.position[2] != lengths[2]-1) {
      soma.axons[0].position = {rand() % Config::X_SIZE, rand() % Config::Y_SIZE, Config::Z_SIZE-1};
      soma.axons[1].position = {rand() % Config::X_SIZE, rand() % Config::Y_SIZE, Config::Z_SIZE-1};
      for (unsigned int a=2; a<soma.axons.size(); a++) {
        int x = rand() % Config::X_SIZE;
        int y = rand() % Config::Y_SIZE;
        if (!(x == 0 || x == Config::X_SIZE-1)) {
          y = rand()>RAND_MAX/2.0 ? 0 : Config::Y_SIZE-1;
        }
        soma.axons[a].position = {x, y, 4};
        //soma.axons[a].position = somas[rand() % somas.size()].position;
      }
    }
  }
  srand(time(NULL));
}

void Environment::set_morphogens() {
  double emission = 0.0;
  double distance = 0.0;
  double morph = 0.0;
  int position;
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

void Environment::develop_grns(double reward) {
  double soma_signal;
  double soma_concentration;
  for (auto& soma : somas) {
    soma_signal = soma.grn.getProteinConcentration("comm", ProteinType::output);
    soma.evolve(morphogens[soma.position[0]][soma.position[1]][soma.position[2]], reward);
    for (auto& axon : soma.axons) {
      soma_concentration = 0.0;
      auto rec_soma = soma_at(axon.position);
      if (rec_soma != NULL && rec_soma != &soma) {
        soma_concentration = rec_soma->nt_concentration;
      }
      axon.evolve(morphogens[soma.position[0]][soma.position[1]][soma.position[2]],
                  soma_concentration, soma_signal, reward);
    }
  }
}

void Environment::set_nt_concentration(vector<vector<double> > inputs) {
  for (int x=0; x < lengths[0]; x++) {
    for (int y=0; y < lengths[1]; y++) {
      soma_at({x, y, 0})->next_concentration += inputs[y][x];
    }
  }

  for (auto& soma: somas) {
    /*
    double p1 = soma.grn.proteins[SOMA_GRN_OUTPUT_FIRING].concentration;
    double p2 = soma.grn.proteins[SOMA_GRN_OUTPUT_FIRING_THRESH].concentration;
    double vt = 1.0 + std::abs((double)(p1-p2)/(p1+p2));
    double vr = 2.0 - vt;
    */
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
      //std::cout << "fired soma " << soma.id << " ";
      for (auto& axon : soma.axons) {
        auto rec_soma = soma_at(axon.position);
        if (rec_soma != NULL && rec_soma != &soma && rec_soma->position[2] != 0) {
            //&& (soma.position[2] != 0 || rec_soma->position[2] != lengths[2]-1)) {
          double nt_amount = soma.grn.getProteinConcentration("nt", ProteinType::output);
          double nt_thresh = soma.grn.getProteinConcentration("nt_t", ProteinType::output);

          double fire_conc = 0.0;
          if (nt_amount > 0.0 || nt_thresh > 0.0) {
            fire_conc = (nt_amount - nt_thresh)/(nt_amount + nt_thresh);
          }
          rec_soma->next_concentration += axon.fire()/2.0;
          rec_soma->next_concentration += fire_conc/2.0;
        }
      }
    }
  }
  //std::cout << std::endl;
}

void Environment::axon_actions() {
  /*
  int action = 0;
  int size = 0;
  for (auto& soma: somas) {
    for (auto& axon : soma.axons) {
      //cout << axon.grn.toString();
      if (axon.age>0) {
        action = axon.act();
        if (action < N_M) {
          // move the axon along the morphogen gradient
          auto new_pos = move_position(axon.position, action);
          //cout << " move (" << axon.position[0] << "," << axon.position[1] << "," << axon.position[2] << ") ";
          //cout << " to (" << new_pos[0] << "," << new_pos[1] << "," << new_pos[2] << ") ";
          for (int i=0; i<N_D; i++) {
            axon.position[i] = new_pos[i];
          }
        }
      }
      axon.age++;
    }
    size = soma.axons.size();
    for (int i=0; i<size; i++) {
      if (soma.axons[i].marked_for_branch) {
        if (soma.axons.size() < AXON_MAX_NUMBER) {
          soma.axons.push_back(Axon(soma.axons[i].grn, soma.axons[i].position));
          if (soma.axons[i].division_conc >= AXON_DIVISION_REDUCTION) {
            soma.axons[i].division_conc -= AXON_DIVISION_REDUCTION;
          }
        }
      }
    }
    for (int i=0; i<soma.axons.size(); ) {
      if (soma.axons[i].marked_for_deletion) {
        soma.axons.erase(soma.axons.begin()+i);
      } else {
        i++;
      }
    }
  }
  */
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
