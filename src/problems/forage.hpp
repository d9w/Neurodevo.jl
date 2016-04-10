#ifndef FORAGE_H
#define FORAGE_H
#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include "robot.hpp"
#include "../core/config.hpp"

using std::vector;
using std::map;
using std::string;

class Forage {
public:
  using d2 = vector<vector<double>>;

protected:
  int x;
  int y;
  unsigned int tsteps;
  d2 food;
  Robot robot;
  std::mt19937 rand_engine;
  double reward;

public:

  Forage(int seed) {
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

  void setInputs(d2 *inputs) {
    inputs->clear();
    inputs->push_back(vector<double> (Config::X_SIZE, 1.0));
    // inputs[0] is an x-length row at y=0
    for (unsigned int i = 1; i < Config::Y_SIZE; i++) inputs->push_back(vector<double> (Config::X_SIZE, 0.0));
    reward = 0.0;
    for (auto foo : nearby_food()) {
      double dist = 1.0 - (foo[0] / (double) robot.sight);
      double angle = foo[1]; //should be from 0 to M_PI
      int sensor = std::floor(angle/(M_PI/Config::X_SIZE));

      reward = std::max(reward, dist);
      for (unsigned int i=1; i<Config::Y_SIZE; i++) {
        if ((double)i/(Config::Y_SIZE-1) < dist) {
          (*inputs)[i][sensor] = 1.0;
        }
      }
    }
  }

  void setFootprint(vector<double> *history) {
    history->clear();
    history->push_back({robot.x,robot.y});
  }

  void step(const d2 *outputs) {
    assert(static_cast<int>(outputs->size()) == Config::Y_SIZE);
    assert(static_cast<int>(outputs->at(0).size()) == Config::X_SIZE);
    double right = 0;
    double left = 0;
    unsigned int right_bound = std::max(0, static_cast<int>(std::floor(Config::X_SIZE/2.0)-1));
    unsigned int left_bound = std::min(right_bound+2, Config::X_SIZE);
    vector<double> fire_signature;
    for (unsigned int i = 0; i < Config::Y_SIZE; i++) {
      for (unsigned int j = 0; j < Config::X_SIZE; j++) {
        double out = (*outputs)[i][j];
        if (j < left_bound) {
          left += out;
        }
        if (j > right_bound) {
          right += out;
        }
        fire_signature.push_back(out);
      }
    }

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

  void setFitness(map<string, double> *fitnesses) {
    //std::cout << "calling set fitnesses with " << food.size() << "/" << Config::HOTDOGS << " food in " << tsteps << std::endl;
    //(*fitnesses)["consumed"] = (Config::HOTDOGS - food.size())/static_cast<double>(tsteps);
    (*fitnesses)["consumed"] = 1.0 - food.size()/static_cast<double>(Config::HOTDOGS);
  };

  bool stop() {
    return (robot.life <= 0 || tsteps > Config::PROBLEM_EVAL_STEPS);
  }

  double getReward() {
    return reward;
  }

  const d2 nearby_food() {
    d2 nearby;
    bool eat;
    double dist;
    for (unsigned int i=0; i<food.size(); ) {
      eat = false;

      dist = std::sqrt((food[i][0] - robot.x) * (food[i][0] - robot.x) +
                      (food[i][1] - robot.y) * (food[i][1] - robot.y));

      if (dist <= robot.sight) {
        if (dist <= robot.size) {
          eat = true;
        } else {
          double angle = std::atan2(food[i][1] - robot.y, food[i][0] - robot.x);
          if (robot.theta > M_PI) {
            angle -= (robot.theta - 2*M_PI);
          } else {
            angle -= robot.theta;
          }
          //std::cout << "food at " << food[i][0] << "," << food[i][1] << " robot at " << robot.x << " " << robot.y << " theta " << robot.theta << " phi " << angle << std::endl;
          if (angle >= 0 && angle < M_PI) {
            nearby.push_back({dist, angle});
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
};
#endif
