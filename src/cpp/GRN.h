/*
 *  GRN.h
 *  GRN_perso
 *
 *  Created by Sylvain Cusat-Blanc on 19/10/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef GRN_MODEL_H
#define GRN_MODEL_H

#include <vector>
#include <string>
#include <stdio.h>

#include "Protein.h"

using namespace std;

class GRN {
public:
	vector<Protein> proteins;
	double beta;
	double delta;
  double t_action;

	GRN(vector<Protein> prot, double beta, double delta, double t_action);
	GRN();
	GRN(const GRN& g);
	GRN(string fileName);
  GRN(int inp_length, int output_length, int regul_length);

  inline vector<Protein> getProteins() {return proteins;}
//	inline unsigned int getNbProteins() {return nbProteins;};
	void evolve(unsigned int nbSteps);
	void duplicateProteins(int nbDup, double mutProb);
	void duplicateRegulatoryProteins(int nbDup, double mutProb);
	GRN copy();
	void mutate(double mutRate);
	void reset();
	void saveToFile(string fileName);
	string toString();
	void updateSignatures();

  unsigned int currentStep;
	double maxEnhance;
	double maxInhibit;
	vector < vector < double > > enhanceMatching;
	vector < vector < double > > inhibitMatching;

private:
};

#endif
