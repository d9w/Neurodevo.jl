#include "external/gaga/gaga/gaga.hpp"
#include "external/cxxopts/src/cxxopts.hpp"
#include "core/config.hpp"
#include "core/types.hpp"

int main(int argc, char** argv) {
	using environment_t = Types::DNAType;
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
  ga.setEvaluator([](auto &i) { i.fitnesses["length"] = i.dna.getProteinSize(ProteinType::regul); });
  if (novelty) {
    ga.enableNovelty();
    ga.setMinNoveltyForArchive(Config::NOVELTY_MIN);
  }
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
    t.config.ADD_RATE = 0.10;
    t.config.DEL_RATE = 0.10;
    t.config.MODIF_RATE = 0.80;
    t.addRandomProtein(ProteinType::input, "input");
    t.addRandomProtein(ProteinType::output, "output");
    t.randomReguls(1);
    t.randomParams();
    pop.push_back(GAGA::Individual<dna_t>(t));
  }
	ga.setPopulation(pop);
  ga.step(10);
  return 0;
}
