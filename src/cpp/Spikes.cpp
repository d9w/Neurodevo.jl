#include "Spikes.h"
#include "Constants.h"
#include <vector>
#include <cmath>
#include <iostream>

Spikes::Spikes() {};

Spikes::Spikes(int s, int q, int seed) {
  srand(123456+100*seed);
  x = s;
  y = q;
  tsteps = 0;

  for (int i=0; i<HOTDOGS; i++) {
    double food_x = (double)rand()/RAND_MAX * x;
    double food_y = (double)rand()/RAND_MAX * y;
    food.push_back({food_x, food_y});
    //std::cout << "r " << radius << " x " << food_x << " y " << food_y << std::endl;
  }
  robot = Robot((double)x/2.0, (double)y/2.0, 0.0, 2.0, 100.0, 7.5);
  srand(time(NULL));
}

void Spikes::step() {
  // assume robot has already moved
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

  tsteps += 1;
  robot.life -=1;
}

std::vector<std::vector<double> > Spikes::nearby_food() {
  std::vector<std::vector<double> > nearby;
  bool eat;
  double dist;
  for (int i=0; i<food.size(); ) {
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
