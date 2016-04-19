#ifndef FORAGE_H
#define FORAGE_H
#include "Robot.h"
#include <map>

using std::vector;
using std::map;
using std::string;
using d2 = vector<vector<double>>;

class Forage {
 protected:
  int x;
  int y;
  unsigned int tsteps;
  d2 food;
  Robot robot;
  std::mt19937 rand_engine;
  double reward;

public:

  Forage(int seed);

  const d2 getInputs();
  const d2 getFootprint();
  void step(const d2 outputs);
  const map<string, double> getFitness();
  bool stop();
  double getReward() { return reward; }
  const d2 nearby_food();
};
#endif
