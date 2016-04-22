#ifndef DNA_H
#define DNA_H
#include "config.hpp"
#include "../external/grgen/common.h"
#include "../external/grgen/grn.hpp"
#include "../external/grgen/classic.hpp"

class DNA {
 public:
  using GRN_T = GRN<Classic>;
 protected:
  GRN_T grn;

 public:
  int foo = 0;

  DNA();
  DNA(const GRN_T &g);
  DNA(const std::string &s);
  static DNA random();
  DNA crossover(const DNA &other);

  void update();
  void mutate();
  void reset();
  void setInput(const std::string &input, double val);
  double getOutput(const std::string &output);
  std::string concString() const;
  std::string toJSON() const;
};
#endif
