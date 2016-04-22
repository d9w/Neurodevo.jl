#define CLUSTER
#include "external/gaga/gaga.hpp"
#include "core/config.hpp"
#include "core/DNA.h"
#include "core/ANN.h"
#include "problems/Forage.h"

int main(int, char**) {

	std::string evaluatorName;
  bool novelty = false;

  GAGA::GA<DNA> ga(0, nullptr);

  ga.setEvaluator([](auto &i) {
      Forage forage(0);
      ANN ann(i.dna, 0);
      vector<vector<double>> outputs;

      while (!forage.stop()) {
        auto ins = forage.getInputs();
        ann.step(ins, forage.getReward());
        ann.set_outputs(&outputs);
        forage.step(outputs);
      }

      for (auto& fit : forage.getFitness()) i.fitnesses[fit.first] = fit.second;
      i.footprint.clear();
      i.footprint = forage.getFootprint();
    });

  if (novelty) {
    ga.enableNovelty();
    ga.setMinNoveltyForArchive(Config::NOVELTY_MIN);
  }
  ga.setSaveFolder("evos");
  ga.setVerbosity(1);
  ga.setPopSize(Config::NUM_POP);
  ga.setMutationProba(Config::MUTATION_RATE);
  ga.setCrossoverProba(Config::CROSSOVER_RATE);
  ga.setNbElites(Config::NUMBER_ELITES);
  ga.setNbSavedElites(0);
  ga.setPopSaveInterval(10);
  ga.setSelectionMethod(GAGA::SelectionMethod::paretoTournament);
  ga.setTournamentSize(Config::TOURNAMENT_SIZE);


  ga.initPopulation(DNA::random);

  /*
  vector<GAGA::Individual<DNA>> pop;
  for (unsigned int i = 0; i < Config::NUM_POP; ++i) {
    pop.push_back(GAGA::Individual<dna_t>(t.random_dna()));
  }
	ga.setPopulation(pop);
  */

  ga.step();
  return 0;
}
