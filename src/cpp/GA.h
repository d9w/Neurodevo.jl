#ifndef __GA_H__
#define __GA_H__

#include "GRN.h"
#include "MNIST.h"

class GA {
public:
  std::vector<double> fits;
  std::vector<GRN> grns;
  std::vector<std::vector<double>> mnist_images;
  std::vector<int> mnist_labels;
  string save_dir;
  MNIST mnist;
  int world_rank;
  int world_size;

  GA(string s_dir, int world_size, int world_rank);
  ~GA();

  void save_all(int i);
  double env_fit(int iter, int a, int s, int seed);
  double length_fit(int iter, int a, int s);
  void evaluate(int i);
  void mutate(GRN *child);
  GRN crossover(GRN& parent1, GRN& parent2);
  vector<int> select(std::vector<double> *pop_fits);
  void run();
};

#endif
