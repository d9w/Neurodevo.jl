#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include "Forage.h"
#include "../core/config.hpp"

using std::vector;
using std::map;
using std::string;

Forage::Forage(int seed) {
  rand_engine = std::mt19937(seed);
  x = 50;
  y = 50;
  tsteps = 0;

  std::uniform_real_distribution<double> dist(0.0, 1.0);

  for (unsigned int i=0; i<Config::HOTDOGS; i++) {
    double food_x = dist(rand_engine) * x;
    double food_y = dist(rand_engine) * y;
    food.push_back({food_x, food_y});
  }
  robot = Robot((double)x/2.0, (double)y/2.0, 0.0, Config::ROBOT_SIZE, Config::ROBOT_LIFE, Config::ROBOT_SIGHT);
}

const d2 Forage::getInputs() {
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  d2 inputs;
  inputs.push_back(vector<double> (Config::X_SIZE, 1.0));
  // inputs[0] is an x-length row at y=0
  for (unsigned int i = 1; i < Config::Y_SIZE; i++) inputs.push_back(vector<double> (Config::X_SIZE, 0.0));
  reward = 0.0;
  for (auto foo : nearby_food()) {
    double distance = 1.0 - (foo[0] / (double) robot.sight);
    double angle = foo[1]; //should be from 0 to M_PI
    int sensor = std::floor(angle/(M_PI/Config::X_SIZE));

    reward = std::max(reward, distance);
    for (unsigned int i=1; i<Config::Y_SIZE; i++) {
      if (dist(rand_engine) > 0.5) {
        //if ((double)i/(Config::Y_SIZE-1) < distance) {
        inputs[i][sensor] = 1.0;
      }
    }
  }
  return inputs;
}

const d2 Forage::getFootprint() {
  d2 history;
  history.push_back(vector<double>{robot.x,robot.y});
  return history;
}

void Forage::step(const d2 outputs) {
  assert(static_cast<int>(outputs.size()) == Config::Y_SIZE);
  assert(static_cast<int>(outputs.at(0).size()) == Config::X_SIZE);
  double right = 0;
  double left = 0;
  int right_count = 0;
  int left_count = 0;
  unsigned int right_bound = std::max(0, static_cast<int>(std::floor(Config::X_SIZE/2.0)-1));
  unsigned int left_bound = std::min(right_bound+2, Config::X_SIZE);
  vector<double> fire_signature;
  for (unsigned int i = 0; i < Config::Y_SIZE; i++) {
    for (unsigned int j = 0; j < Config::X_SIZE; j++) {
      double out = outputs[i][j];
      if (j < left_bound) {
        left += out;
        left_count += 1;
      }
      if (j > right_bound) {
        right += out;
        right_count += 1;
      }
      fire_signature.push_back(out);
    }
  }

  right /= right_count;
  left /= left_count;

  robot.move(right, left);

  while (robot.x > x) {
    robot.x -= x;
  }
  while (robot.x < 0) {
    robot.x += x;
  }
  while (robot.y > y) {
    robot.y -= y;
  }
  while (robot.y < 0) {
    robot.y += y;
  }
  while (robot.theta < 0) {
    robot.theta += 2*M_PI;
  }
  while (robot.theta > 2*M_PI) {
    robot.theta -= 2*M_PI;
  }

  //std::cout << tsteps << " " << robot << " food: " << food.size() <<std::endl;

  tsteps += 1;
  robot.life -=1;
}

const map<string, double> Forage::getFitness() {
  //std::cout << "calling set fitnesses with " << food.size() << "/" << Config::HOTDOGS << " food in " << tsteps << std::endl;
  //(*fitnesses)["consumed"] = (Config::HOTDOGS - food.size())/static_cast<double>(tsteps);
  map<string, double> fitnesses;
  fitnesses["consumed"] = 1.0 - food.size()/static_cast<double>(Config::HOTDOGS);
  return fitnesses;
}

bool Forage::stop() {
  return (robot.life <= 0 || tsteps > Config::PROBLEM_EVAL_STEPS);
}

const d2 Forage::nearby_food() {
  d2 nearby;
  bool eat;
  double distance;

  for (unsigned int i=0; i<food.size(); ) {
    eat = false;
    distance = std::sqrt((food[i][0] - robot.x) * (food[i][0] - robot.x) +
                    (food[i][1] - robot.y) * (food[i][1] - robot.y));

    if (distance <= robot.sight) {
      if (distance <= robot.size) {
        eat = true;
      } else {
        double angle = std::atan2(food[i][1] - robot.y, food[i][0] - robot.x);
        if (robot.theta > M_PI) {
          angle -= (robot.theta - 2*M_PI);
        } else {
          angle -= robot.theta;
        }
        if (angle >= 0 && angle < M_PI) {
          nearby.push_back({distance, angle});
        }
      }
    }

    if (eat) {
      food.erase(food.begin() + i);
      robot.life = robot.max_life;
    } else {
      i++;
    }
  }
  return nearby;
}
