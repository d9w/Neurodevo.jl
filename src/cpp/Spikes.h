#ifndef SPIKES_H
#define SPIKES_H

#include <vector>
#include "Robot.h"

class Spikes {
 public:
  int x;
  int y;
  int tsteps;
  std::vector<std::vector<double> > food;
  Robot robot;

  Spikes();
  Spikes(int s, int q, int seed);

  void step();
  std::vector<std::vector<double> > nearby_food();
};

#endif
