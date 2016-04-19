#ifndef DNA_H
#define DNA_H
#include "config.hpp"
#include "../external/grgen/common.h"
#include "../external/grgen/grn.hpp"
#include "../external/grgen/classic.hpp"

class DNA {
 private:
  vector<string> inputs;
  vector<string> outputs;

 public:
  using GRN_T = GRN<Classic>;

  GRN_T grn;

  DNA();

  DNA(const GRN_T &g) : grn(g) { reset(); }

  DNA(const std::string &s) : grn(s) { reset(); }

  DNA(const DNA &other) : grn(other.grn) {}

	DNA &operator=(const DNA &other) {
		if (this != &other) {
			grn = other.grn;
			reset();
		}
		return *this;
	}

  const DNA random();

	void update() { grn.step(Config::GRN_EVO_STEPS); }

	DNA crossover(const DNA &other) {
		GRN_T g = grn.crossover(other.grn);
		DNA res(g);
		return res;
	}

	void mutate() { grn.mutate(); }

	void reset() { grn.reset(); }

	void setInput(const std::string &input, double val) {
		grn.setProteinConcentration(input, ProteinType::input, val);
	}

	double getOutput(const std::string &output) const {
		auto r = grn.getProteinConcentration(output, ProteinType::output);
		return r;
	}

	std::string toJSON() const { return grn.toJSON(); }

  std::string concString() const;

};
#endif
