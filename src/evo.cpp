#define CLUSTER
#include "external/gaga/gaga/gaga.hpp"
#include "external/cxxopts/src/cxxopts.hpp"
#include "core/config.hpp"
#include "core/types.hpp"
#include "core/evaluator.hpp"
#include "problems/forage.hpp"

int main(int argc, char** argv) {
	using dna_t = Types::DNAType;

	std::string evaluatorName;
  bool novelty;
	try {
		cxxopts::Options options(argv[0]);
		options.add_options()("e,evaluator", "evaluator name",
		                      cxxopts::value<std::string>(evaluatorName));
		options.add_options()("n,novelty", "enable novelty");
		options.parse(argc, argv);
    novelty = options["novelty"].as<bool>();
	} catch (const cxxopts::OptionException& e) {
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	} catch (const std::bad_cast& e) {
		std::cout << "bad cast: " << e.what() << std::endl;
		exit(1);
	}
	if (evaluatorName == "length" || evaluatorName == "") {
    std::cout << "length evaluator" << std::endl;
  } else {
    std::cerr << "No valid evaluator found, aborting." << std::endl;
    exit(1);
  }

  GAGA::GA<dna_t> ga(0, nullptr);
  ga.setEvaluator([](auto &i) {
      Forage forager(0);
      Evaluator<Forage> eval;
      eval.evaluate(forager, i.dna, 0);
      for (auto& fit : *eval.getFitnesses()) i.fitnesses[fit.first] = fit.second;
      i.footprint.clear();
      i.footprint = (*eval.getHistory());
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
  vector<GAGA::Individual<dna_t>> pop;
  for (unsigned int i = 0; i < Config::NUM_POP; ++i) {
    dna_t t;
    t.addRandomProtein(ProteinType::input, "x");
    t.addRandomProtein(ProteinType::input, "y");
    t.addRandomProtein(ProteinType::input, "z");
    t.addRandomProtein(ProteinType::input, "nt");
    t.addRandomProtein(ProteinType::input, "comm");
    t.addRandomProtein(ProteinType::input, "div");
    t.addRandomProtein(ProteinType::input, "reward");

    t.addRandomProtein(ProteinType::output, "nt");
    t.addRandomProtein(ProteinType::output, "nt_t");
    t.addRandomProtein(ProteinType::output, "f");
    t.addRandomProtein(ProteinType::output, "f_t");
    t.addRandomProtein(ProteinType::output, "comm");

    t.randomReguls(1);
    t.randomParams();

    pop.push_back(GAGA::Individual<dna_t>(t));
  }
	ga.setPopulation(pop);
  ga.step(Config::GENERATIONS);
  return 0;
}
