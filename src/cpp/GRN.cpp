#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <time.h>

#include "GRN.h"

using namespace std;

GRN::GRN(vector<Protein> p, double b, double d, double t) {
	proteins=p;
	beta=b;
	delta=d;
  t_action=t;

	currentStep=0;

	maxEnhance=0;
	maxInhibit=0;

	updateSignatures();
}

GRN::GRN() {
	proteins=vector<Protein>();
	beta=1;
	delta=1;
  t_action=(T_ACTION_MIN+T_ACTION_MAX)/2.0;
	currentStep=0;
}

GRN::GRN(const GRN& g) {
	proteins=vector<Protein>();
	for (unsigned int i=0; i<g.proteins.size(); i++) {
		proteins.push_back(Protein(g.proteins[i]));
	}

	beta=g.beta;
	delta=g.delta;
  t_action=g.t_action;
	currentStep=g.currentStep;

	maxEnhance=g.maxEnhance;
	maxInhibit=g.maxInhibit;

	updateSignatures();
}

GRN::GRN(int inp_length, int output_length, int regul_length) {
	proteins=vector<Protein>();
  for (int i=0; i<(inp_length + output_length + regul_length); i++) {
    Protein protein;
    protein.id=(int)rand()%ID_SIZE;
    protein.enhancer=(int)rand()%ID_SIZE;
    protein.inhibiter=(int)rand()%ID_SIZE;
    if (i < inp_length) {
      protein.type = Protein::PROTEIN_INPUT;
    } else if (i < (inp_length+output_length)) {
      protein.type = Protein::PROTEIN_OUTPUT;
    } else {
      protein.type = Protein::PROTEIN_REGUL;
    }
    protein.concentration = 0.0;
    proteins.push_back(protein);
  }
  beta=(double)rand()/RAND_MAX*1.5+0.5;
  delta=(double)rand()/RAND_MAX*1.5+0.5;
  t_action=(double)rand()/RAND_MAX*(T_ACTION_MAX-T_ACTION_MIN)+T_ACTION_MIN;
  currentStep = 0;
	maxEnhance=0;
	maxInhibit=0;
  updateSignatures();
}

void GRN::updateSignatures() {
	// calculating signatures
  enhanceMatching.clear();
  inhibitMatching.clear();
  //cout << "Protein size " << proteins.size() << endl;
	for (unsigned int j=0; j<proteins.size(); j++) {
		vector<double> enh;
		vector<double> inh;
		for (unsigned int k=0; k<proteins.size(); k++) {
			enh.push_back(ID_SIZE-abs(proteins[j].enhancer-proteins[k].id));
			maxEnhance=std::max(maxEnhance, enh[k]);
			inh.push_back(ID_SIZE-abs(proteins[j].inhibiter-proteins[k].id));
			maxInhibit=std::max(maxInhibit, inh[k]);
      //cout << "Enhancer: " << proteins[j].enhancer << endl;
      //cout << "ID: " << proteins[k].id << endl;
      //cout << "maxEnhance: " << maxEnhance << endl;
      //cout << "maxInhibit: " << maxInhibit << endl;
		}
    //cout << "enh: " << enh.size() << endl;
    //cout << "inh: " << inh.size() << endl;
		enhanceMatching.push_back(enh);
		inhibitMatching.push_back(inh);
	}
  //cout << "enhance: " << enhanceMatching.size() << " " << enhanceMatching[0].size();
  //cout << "inhibit: " << inhibitMatching.size() << " " << inhibitMatching[0].size();
	for (unsigned int j=0; j<proteins.size(); j++) {
		for (unsigned int k=0; k<proteins.size(); k++) {
			enhanceMatching[j][k]=exp(beta*(enhanceMatching[j][k]-maxEnhance));
			inhibitMatching[j][k]=exp(beta*(inhibitMatching[j][k]-maxInhibit));
		}
	}
}


void GRN::evolve(unsigned int nbSteps) {
	double enhance, inhibit;
	unsigned int i,j,k,nbDiv;
	double sumConcentration;
//	for (i=0; i<proteins.size(); i++) {
//		cerr << proteins[i].concentration << "\t";
//	}
//	cerr << proteins.size() << endl;


	for (unsigned int step=0; step<nbSteps; step++) {
		vector<Protein> nextProteins=vector<Protein>();

		// Calculating next step protein concentrations
		for (j=0; j<proteins.size(); j++) {
			// For an input protein, the concentration is not calculated
			if (proteins[j].type==Protein::PROTEIN_INPUT) {
				nextProteins.push_back(proteins[j]);
			} else {
				enhance=0;
				inhibit=0;
				// Calculating enhancing and inhibiting factor for the current protein
				for (k=0; k<proteins.size(); k++) {
					if (proteins[k].type!=Protein::PROTEIN_OUTPUT) {
						enhance+=proteins[k].concentration*enhanceMatching[j][k];
						inhibit+=proteins[k].concentration*inhibitMatching[j][k];
				    }
				}
				// if (j=5) //cout << enhance << "   " << inhibit << endl;
				// Calculating the next concentration of current protein
				nextProteins.push_back(Protein(proteins[j].id, proteins[j].type, max(0.0,proteins[j].concentration+delta/proteins.size()*(enhance-inhibit)), proteins[j].enhancer, proteins[j].inhibiter));
			}
		}

		// Reajusting proteins concentration so that their sum stay equal to 1
		sumConcentration = 0;
		nbDiv=0;
		for (i=0; i<proteins.size(); i++) {
			if (proteins[i].type!=Protein::PROTEIN_INPUT) {
				sumConcentration+=nextProteins[i].concentration;
				nbDiv++;
			}
		}
		////cout << sumConcentration << "\t" << nbDiv << "\t" << (sumConcentration-1)/nbDiv << "\t";
		for (i=0; i<proteins.size(); i++) {
			if (proteins[i].type!=Protein::PROTEIN_INPUT) {
        if (sumConcentration > 0.0) {
          nextProteins[i].concentration/=sumConcentration;
        } else {
          nextProteins[i].concentration = 0.0;
        }
			}
		}

		// Finalizing the step by switching the vector
		proteins.clear();
		proteins=nextProteins;

		//for (i=0; i<nextProteins.size(); i++) {
		//	//cout << proteins[i].concentration << "\t";
		//}
		////cout << endl;
		currentStep++;
	}
}

void GRN::duplicateProteins(int nbDup, double mutProb) {
	//srand(time(NULL))
	int nbProteins=proteins.size();
	for (int dup=0; dup<nbDup; dup++) {
		for (int i=0; i<nbProteins; i++) {
				proteins.push_back(Protein(
									   (double)rand()/(double)RAND_MAX>mutProb?proteins[i].id:rand()/ID_SIZE,
									   proteins[i].type,
									   proteins[i].concentration,
									   (double)rand()/(double)RAND_MAX>mutProb?proteins[i].enhancer:rand()/ID_SIZE,
									   (double)rand()/(double)RAND_MAX>mutProb?proteins[i].inhibiter:rand()/ID_SIZE));
		}
	}
	updateSignatures();
}

void GRN::duplicateRegulatoryProteins(int nbDup, double mutProb) {
	//srand(time(NULL))
	int nbProteins=proteins.size();
	for (int dup=0; dup<nbDup; dup++) {
		for (int i=0; i<nbProteins; i++) {
			if (proteins[i].type==Protein::PROTEIN_REGUL) {
				proteins.push_back(Protein(
										   (double)rand()/(double)RAND_MAX>mutProb?proteins[i].id:rand()/ID_SIZE,
										   proteins[i].type,
										   proteins[i].concentration,
										   (double)rand()/(double)RAND_MAX>mutProb?proteins[i].enhancer:rand()/ID_SIZE,
										   (double)rand()/(double)RAND_MAX>mutProb?proteins[i].inhibiter:rand()/ID_SIZE));
			}
		}
	}
	updateSignatures();
}

GRN GRN::copy() {
	GRN res = GRN();
	res.proteins=vector<Protein>();
	for (unsigned int i=0; i<this->proteins.size(); i++) {
		res.proteins.push_back(Protein(this->proteins[i]));
	}

	res.beta=this->beta;
	res.delta=this->delta;
  res.t_action=this->t_action;
	res.currentStep=this->currentStep;

	res.maxEnhance=this->maxEnhance;
	res.maxInhibit=this->maxInhibit;

	res.enhanceMatching=vector < vector < double > >();
	for (unsigned int i=0;i<this->enhanceMatching.size(); i++) {
		vector<double> l = vector<double>();
		for (unsigned int j=0; j<this->enhanceMatching[i].size(); j++) {
			l.push_back(this->enhanceMatching[i][j]);
		}
		res.enhanceMatching.push_back(l);
	}
	res.inhibitMatching=vector < vector < double> >();
	for (unsigned int i=0;i<this->inhibitMatching.size(); i++) {
		vector<double> l = vector<double>();
		for (unsigned int j=0; j<this->inhibitMatching[i].size(); j++) {
			l.push_back(this->inhibitMatching[i][j]);
		}
		res.inhibitMatching.push_back(l);
	}

	return res;
}

void GRN::mutate(double mutRate) {
  //cout << "Begin mutation" << endl;
  //cout << proteins.size() << endl;
	for (int i=0; i<proteins.size(); i++) {
    //cout << "protein " << i << " id " << proteins[i].id << " enh " << proteins[i].enhancer << "inh" << proteins[i].inhibiter << endl;
		proteins[i].id=(double)rand()/RAND_MAX<mutRate?rand()%ID_SIZE:proteins[i].id;
		proteins[i].enhancer=(double)rand()/RAND_MAX<mutRate?rand()%ID_SIZE:proteins[i].enhancer;
		proteins[i].inhibiter=(double)rand()/RAND_MAX<mutRate?rand()%ID_SIZE:proteins[i].inhibiter;
    //cout << "protein " << i << " now id " << proteins[i].id << " enh " << proteins[i].enhancer << "inh" << proteins[i].inhibiter << endl;
	}
  //cout << "Mutated the proteins" << endl;
	beta=(double)rand()/RAND_MAX<mutRate?((double)rand()/RAND_MAX*1.5+0.5):beta;
	delta=(double)rand()/RAND_MAX<mutRate?((double)rand()/RAND_MAX*1.5+0.5):delta;
	t_action=(double)rand()/RAND_MAX<mutRate?((double)rand()/RAND_MAX*(T_ACTION_MAX-T_ACTION_MIN)+T_ACTION_MIN):t_action;
  //cout << "Mutated the parameters to " << beta << " and " << delta << endl;
	updateSignatures();
  //cout << "Updated the signatures" << endl;
}

void GRN::reset() {
	for (int i=0; i<proteins.size(); i++) {
		proteins[i].concentration=1.0/(double)proteins.size();
	}
	currentStep=0;
}

void GRN::saveToFile(string fileName) {
	ofstream file;
	file.open(fileName.c_str());

	file.write((char*)&beta, sizeof(double));
	file.write((char*)&delta, sizeof(double));
	file.write((char*)&t_action, sizeof(double));
	file.write((char*)&currentStep, sizeof(int));
	int nbProt = proteins.size();
	file.write((char*)&nbProt, sizeof(int));

	for (int i=0; i<proteins.size(); i++) {
		file.write((char*)&proteins[i].id, sizeof(int));
		file.write((char*)&proteins[i].concentration, sizeof(double));
		file.write((char*)&proteins[i].type, sizeof(int));
		file.write((char*)&proteins[i].enhancer, sizeof(int));
		file.write((char*)&proteins[i].inhibiter, sizeof(int));
	}

	file.close();
}

GRN::GRN(string fileName) {
	ifstream file;
	file.open(fileName.c_str());

	file.read((char*)&beta,sizeof(double));
	file.read((char*)&delta,sizeof(double));
	file.read((char*)&t_action,sizeof(double));
	file.read((char*)&currentStep,sizeof(int));

	int nbProt;
	file.read((char*)&nbProt,sizeof(int));

	int protID;
	double protConc;
	int protType;
	int protEN;
	int protIN;
	for (int i=0; i<nbProt; i++) {
		file.read((char*)&protID,sizeof(int));
		file.read((char*)&protConc,sizeof(double));
		file.read((char*)&protType,sizeof(int));
		file.read((char*)&protEN,sizeof(int));
		file.read((char*)&protIN,sizeof(int));
		proteins.push_back(Protein(protID, protType, protConc, protEN, protIN));
	}
	file.close();
	maxEnhance=0;
	maxInhibit=0;
	updateSignatures();
}


string GRN::toString() {
	stringstream res;
	res << "GRN\tbeta=" << beta << " ; delta=" << delta << " ; t_action=" << t_action << endl;
	res << "\tcurrent step=" << currentStep << " ; nbProt=" << proteins.size() << endl;
	res << "\tProteins:" << endl;
	for (int i=0; i<proteins.size(); i++) {
		res << "\t\tProtein " << i << ":\t[" << proteins[i].id << " ; " << proteins[i].concentration << " ; " << proteins[i].type << " ; " << proteins[i].enhancer << " ; " << proteins[i].inhibiter << "]" << endl;
	}
	return res.str();
}
