#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <mpi.h>
#include <limits>
#include <cmath>
#include "easylogging++.h"
#include "color.h"
#include "Constants.h"
#include "GA.h"
#include "GRN.h"
#include "Environment.h"
#include "Spikes.h"

using namespace std;

GA::GA(string s_dir, int w_size, int w_rank) {

  save_dir = s_dir;

  world_size = w_size;
  world_rank = w_rank;

  if (world_rank == 0) {
    for (int i=0; i<NUM_POP; i++) {
      //axons.push_back(GRN(s_dir+"/input_axon_"+std::to_string(i)+".grn"));
      //somas.push_back(GRN(s_dir+"/input_soma_"+std::to_string(i)+".grn"));
      grns.push_back(GRN(GRN_INPUTS, GRN_OUTPUTS, GRN_REGULS));
      fits.push_back(0.0);
    }
  }

  evaluate(0);
}

GA::~GA() {
}

void GA::save_all(int i) {
  // save grns
  for (int a=0; a<NUM_POP; a++) {
    grns[a].saveToFile(save_dir+"/"+std::to_string(a)+"_"+std::to_string(i)+".grn");
  }
}

double GA::env_fit(int iter, int a, int s, int seed) {

  vector<int> lengths = {X_SIZE, Y_SIZE, Z_SIZE};
  vector<vector<double> > inputs;
  double fit = 0.0;
  GRN grn = GRN(save_dir+"/"+std::to_string(a)+"_"+std::to_string(iter)+".grn");
  Environment env(lengths, grn, grn);
  env.set_random_connectivity(seed);
  Spikes spikes(50, 50, seed);

  while (spikes.tsteps < GA_EVAL_STEPS && spikes.robot.life >= 0) {
    int step = spikes.tsteps;
    inputs.clear();
    inputs.push_back({1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
    for (int i=1; i<X_SIZE; i++) {
      inputs.push_back({0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    }
    auto nearby_food = spikes.nearby_food();
    double reward = 0.0;
    for (auto food : nearby_food) {
      double dist = 1.0 - (food[0] / (double) spikes.robot.sight);
      double angle = food[1];
      int sensor = std::floor(angle/(M_PI/8.0));

      reward = std::max(reward, 1.0-dist);
      for (int i=1; i<X_SIZE; i++) {
        if ((double)i/(X_SIZE-1) < dist) {
          inputs[i][sensor] = 1.0;
        }
      }
    }

    env.develop_grns(reward);

    env.set_nt_concentration(inputs);

    env.fire_ann();

    double right = 0;
    double left = 0;
    for (auto& soma : env.somas) {
      if (soma.position[2] == lengths[2]-1) {
        if (soma.fired) {
          if (soma.position[0] < 5) {
            left += 1.0;
          }
          if (soma.position[0] > 2) {
            right += 1.0;
          }
        }
      }
    }

    right /= (8.0*20.0);
    left /= (8.0*20.0);

    spikes.robot.move(right, left);
    spikes.step();
  }

   LOG(DEBUG) << "rank: " << world_rank << " iter: " << iter
              << " ate: " << (HOTDOGS-spikes.food.size()) << " survived: " << spikes.tsteps;

   return (HOTDOGS-spikes.food.size())/(double)spikes.tsteps;
   //return (double)spikes.tsteps/(50*(HOTDOGS+1));
}

void GA::evaluate(int iter_count) {
  if (world_rank == 0) {
    save_all(iter_count);
  }

  int eval_num;
  MPI_Barrier(MPI_COMM_WORLD);

  int a = world_rank;
  int s = world_rank;

  LOG(INFO) << "rank: " << world_rank << " iter: " << iter_count
            << " evaluating [" << a << "][" << s << "]: ";

  double fit = 2.0;

  for (int i=0; i<5; i++) {
    fit = std::min(fit, env_fit(iter_count, a, s, i));
  }

  LOG(INFO) << "rank: " << world_rank << " iter: " << iter_count
            << " fit: " << fit;

  MPI_Barrier(MPI_COMM_WORLD);
  if (world_rank !=0) {
    MPI_Send(&fit, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  } else {
    fits[0] = fit;
    for (int source = 1; source < world_size; ++source) {
      MPI_Recv(&fit, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      fits[source] = fit;
    }
  }

  if (world_rank == 0) {
    //printing
    string fit_str = " fits: ";
    for (auto& f : fits) {
      if (f > 0.2) {
        fit_str += GREEN + std::to_string(f) + RESET + BOLD + " ";
      } else {
        fit_str += std::to_string(f) + " ";
      }
    }
    LOG(INFO) << BOLD << "iter: " << iter_count << fit_str << RESET;
  }
}

void GA::mutate(GRN *child) {
  // add, delete, or mutate
  if (((double)rand()/RAND_MAX) < ADD_MUTATION_RATE) {
    // add a gene
    Protein protein;
    protein.id = (int)rand()%ID_SIZE;
    protein.enhancer = (int)rand()%ID_SIZE;
    protein.inhibiter = (int)rand()%ID_SIZE;
    protein.type = Protein::PROTEIN_REGUL;
    protein.concentration = 0.0;
    child->proteins.push_back(protein);
  }
  if (((double)rand()/RAND_MAX) < DELETE_MUTATION_RATE) {
    // delete a gene
    int reg_genes = 0;
    // this could be a GRN method
    for (auto p : child->proteins) {
      if (p.type == Protein::PROTEIN_REGUL) {
        reg_genes++;
      }
    }
    if (reg_genes > 0) {
      child->proteins.erase(child->proteins.begin()+child->proteins.size()-(1+(rand() % reg_genes)));
    }
  }
  child->mutate(MUTATION_RATE);
}

//assumes parent1 and parent2 have the same number of input and output genes
GRN GA::crossover(GRN& parent1, GRN& parent2) {
  vector<Protein> proteins;
  vector<Protein> regs1;
  vector<Protein> regs2;
  double beta = rand()>RAND_MAX/2.0 ? parent1.beta : parent2.beta;
  double delta = rand()>RAND_MAX/2.0 ? parent1.delta : parent2.delta;
  double t_action = rand()>RAND_MAX/2.0 ? parent1.t_action : parent2.t_action;

  for (auto i=0; i<std::max(parent1.proteins.size(), parent2.proteins.size()); i++) {
    if (i<std::min(parent1.proteins.size(), parent2.proteins.size())) {
      auto prot = parent1.proteins[i];
      if (prot.type == Protein::PROTEIN_INPUT || prot.type == Protein::PROTEIN_OUTPUT) {
        if (rand()>RAND_MAX/2.0) {
          proteins.push_back(parent1.proteins[i].copy());
        } else {
          proteins.push_back(parent2.proteins[i].copy());
        }
      } else {
        regs1.push_back(parent1.proteins[i].copy());
        regs2.push_back(parent2.proteins[i].copy());
      }
    } else {
      if (i < parent1.proteins.size()) {
        regs1.push_back(parent1.proteins[i].copy());
      }
      if (i < parent2.proteins.size()) {
        regs2.push_back(parent2.proteins[i].copy());
      }
    }
  }

  for (int i=0; i<regs1.size(); i++) {
    int index = rand() % (i+1);
    auto temp = regs1[index];
    regs1[index] = regs1[i];
    regs1[i] = temp;

  }

  for (int i=0; i<regs2.size(); i++) {
    int index = rand() % (i+1);
    auto temp = regs2[index];
    regs2[index] = regs2[i];
    regs2[i] = temp;
  }

  for (int i=0; i<regs1.size(); ) {
    bool matched = false;
    for (int j=0; j<regs2.size(); j++) {
      if (regs1[i].getDistance(regs2[j]) < 0.15) {
        if ((double) rand()/RAND_MAX < 0.5) {
          proteins.push_back(regs1[i].copy());
        } else {
          proteins.push_back(regs2[j].copy());
        }
        matched = true;
        regs1.erase(regs1.begin()+i);
        regs2.erase(regs2.begin()+j);
        break;
      }
    }
    if (!matched) {
      i++;
    }
  }

  double add_rand = (double)rand()/RAND_MAX;
  if (add_rand < 1.0/3.0) {
    for (int i=0; i<regs1.size(); i++) {
      proteins.push_back(regs1[i].copy());
    }
  } else if (add_rand < 2.0/3.0) {
    for (int i=0; i<regs2.size(); i++) {
      proteins.push_back(regs2[i].copy());
    }
  }

  return GRN(proteins, beta, delta, t_action);
}

// fits needs to be %tour_size
vector<int> GA::select(vector<double> *pop_fits) {
  int num_winners = pop_fits->size()/TOURNEY_SIZE;
  int competitors [pop_fits->size()];
  vector<int> winners;

  for (int c=0; c<pop_fits->size(); c++) {
    competitors[c] = c;
  }

  for (int c=0; c<pop_fits->size(); c++) {
    int index = rand() % (c+1);
    int temp = competitors[index];
    competitors[index] = competitors[c];
    competitors[c] = temp;
  }


  for (int t=0; t<num_winners; t++) {
    int winner = -1;
    double winner_fit = -std::numeric_limits<double>::max();
    for (int c=0; c<TOURNEY_SIZE; c++) {
      int competitor = competitors[TOURNEY_SIZE*t + c];
      if (pop_fits->at(competitor) > winner_fit) {
        winner = competitor;
        winner_fit = pop_fits->at(winner);
      }
    }
    winners.push_back(winner);
  }
  return winners;
}

void GA::run() {

  int num_winners = NUM_POP/TOURNEY_SIZE;

  for (int i=1; i<MAX_ITER; i++) {

    if (world_rank == 0) {

      auto winner_inds = select(&fits);

      vector<GRN> children;

      for (int c=0; c<(NUM_POP-num_winners); c++) {
        int s1 = rand() % winner_inds.size();
        int s2 = std::max(0, (int)(rand() % winner_inds.size())-1);
        if (s2 == s1) {
          s2++;
        }
        auto child = crossover(grns[winner_inds[s1]], grns[winner_inds[s2]]);
        children.push_back(child);
      }

      // mutate
      for (int c=0; c<children.size(); c++) {
        mutate(&(children[c]));
      }

      // elitism
      for (int c=0; c<num_winners; c++) {
        children.push_back(grns[winner_inds[c]].copy());
      }

      // reset original populations
      grns.clear();
      for (int c=0; c<children.size(); c++) {
        grns.push_back(children[c].copy());
      }
      children.clear();

      for (int a=0; a<NUM_POP; a++) {
        fits[a] = 0.0;
      }
    }
    evaluate(i);
  }
}
