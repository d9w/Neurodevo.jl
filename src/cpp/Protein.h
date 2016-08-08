#ifndef PROTEIN_H
#define PROTEIN_H

#include <cmath>
#include <vector>
#include <stdlib.h>

#include "Constants.h"

class Protein {
public:
	static const int PROTEIN_INPUT = 0;
	static const int PROTEIN_OUTPUT = 1;
	static const int PROTEIN_REGUL= 2;

	int id;
	double concentration;
	int enhancer;
	int inhibiter;
	int type;

	Protein();
	Protein(int ID, int typ, double conc, int enh, int inh);
	Protein(const Protein& p);

    inline int getMatchingValue(Protein p) {return abs(p.concentration-concentration);}
    inline double getConcentration() {return concentration;}
    inline int getID() {return id;}
	Protein copy();
  double getDistance(const Protein& p);

};

#endif
