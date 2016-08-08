#include <iostream>
#include <cstdio>
#include <time.h>
#include <cmath>
#include "easylogging++.h"
#include "color.h"
#include "GRN.h"
#include "Environment.h"
#include "Spikes.h"
#include "Constants.h"

INITIALIZE_EASYLOGGINGPP

int main(int argc, char** argv)
{
  srand(0);

  el::Configurations conf("log.conf");
  el::Loggers::reconfigureAllLoggers(conf);

  vector<int> lengths = {X_SIZE, Y_SIZE, Z_SIZE};

  GRN soma_grn = GRN(argv[1]);
  GRN axon_grn = GRN(argv[2]);
  int t_action = (int)std::round((soma_grn.t_action+axon_grn.t_action)/2.0);
  Environment env (lengths, soma_grn, axon_grn);
  env.set_random_connectivity(2);
  Spikes spikes(50, 50, 2);
  int last_life = spikes.robot.life;

  vector<vector<double> > inputs;

  while (spikes.tsteps < GA_EVAL_STEPS && spikes.robot.life >= 0) {
    //for (int step=0; step<GA_EVAL_STEPS; step++) {
    int step = spikes.tsteps;
    inputs.clear();
    inputs.push_back({1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
    for (int i=1; i<X_SIZE; i++) {
      inputs.push_back({0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    }
    auto nearby_food = spikes.nearby_food();
    for (auto food : nearby_food) {
      double dist = 1.0 - (food[0] / (double) spikes.robot.sight);
      double angle = food[1];
      int sensor = std::floor(angle/(M_PI/8.0));

      for (int i=1; i<X_SIZE; i++) {
        if ((double)i/(X_SIZE-1) < dist) {
          inputs[i][sensor] = 1.0;
        }
      }
    }

    env.develop_grns();

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

    right /= (8.0*5.0);
    left /= (8.0*5.0);

    spikes.robot.move(right, left);
    spikes.step();

    LOG(DEBUG) << "=========================================================================";

    if (spikes.robot.life > last_life) {
      LOG(INFO) << spikes.tsteps << " " << GREEN << spikes.robot << RESET;
    } else {
      LOG(INFO) << spikes.tsteps << " " << spikes.robot;
    }

    for (auto& soma : env.somas) {
      LOG(DEBUG) << spikes.tsteps << " " << soma;
      for (auto& axon : soma.axons) {
        LOG(DEBUG) << spikes.tsteps << " " << axon;
      }
    }

    last_life = spikes.robot.life;
  }

  LOG(INFO) << " ate: " << (HOTDOGS-spikes.food.size()) << " survived: " << spikes.tsteps;
}

