#include "Protein.h"
#include <cmath>

Protein::Protein(int ID, int t, double conc, int enh, int inh) {
	id=ID;
	concentration=conc;
	enhancer=enh;
	inhibiter=inh;
	type=t;
}

Protein::Protein(){
	id=0;
	concentration=0;
	enhancer=0;
	inhibiter=0;
	type=0;
}

Protein::Protein(const Protein& p) {
	id=p.id;
	concentration=p.concentration;
	enhancer=p.enhancer;
	inhibiter=p.inhibiter;
	type=p.type;
}

Protein Protein::copy() {
	return Protein(id, type, concentration, enhancer, inhibiter);
}

double Protein::getDistance(const Protein& p) {
  double a = 0.75;
  double b = 0.125;
  double c = 0.125;
  return (a*std::abs(id-p.id) + b*std::abs(enhancer-p.enhancer) + c*std::abs(inhibiter-p.inhibiter))/
    (double)ID_SIZE;
}
