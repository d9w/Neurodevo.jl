#ifndef CONFIG_HPP
#define CONFIG_HPP
struct Config {
  // grn parametets
  static constexpr unsigned int ID_SIZE = 32;
  static constexpr unsigned int N_D = 3;
  static constexpr unsigned int GRN_STEPS_PER_UPDATE = 1;

  // Environment
  static constexpr unsigned int N_M = 0;
  static constexpr unsigned int X_SIZE = 8;
  static constexpr unsigned int Y_SIZE = 8;
  static constexpr unsigned int Z_SIZE = 8;
  static constexpr unsigned int GRN_EVO_STEPS = 1;
  static constexpr double AXON_DIVISION_START = 0.5;
  static constexpr double AXON_DIVISION_REDUCTION = 0.1;
  static constexpr unsigned int AXON_MAX_NUMBER = 5;
  static constexpr double SOMA_AXON_INPUT_THRESH = 1.0;
  static constexpr double T_ACTION_MIN = 5.0; //range for axon actions
  static constexpr double T_ACTION_MAX = 20.0;

  // firing grn inputs
  static constexpr unsigned int GRN_INPUT_X = N_M; // position x
  static constexpr unsigned int GRN_INPUT_Y = (N_M + 1); // position y
  static constexpr unsigned int GRN_INPUT_Z = (N_M + 2); // position z
  static constexpr unsigned int GRN_INPUT_NEUROTRANSMITTER = (N_M + 3); //neurotransmitter concentration
  static constexpr unsigned int GRN_INPUT_COMM = (N_M + 4); // communication input
  static constexpr unsigned int GRN_INPUT_DIVISION = (N_M + 5); // division concentration, 1 in soma
  static constexpr unsigned int GRN_INPUT_REWARD = (N_M + 6); // division concentration, 1 in soma
  static constexpr unsigned int GRN_INPUTS = (N_M + 7); // number of axon grn inputs
  static constexpr unsigned int N_G_I = GRN_INPUTS;
  // firing grn outpus
  static constexpr unsigned int GRN_OUTPUT_COMM = (N_G_I + N_M); // commuincation output
  static constexpr unsigned int GRN_OUTPUT_F_NULL = (N_G_I + N_M + 1); // firing threshold/null
  static constexpr unsigned int GRN_OUTPUT_FT_NULL = (N_G_I + N_M + 2); // firing threshold thresh/null
  static constexpr unsigned int GRN_OUTPUT_NT = (N_G_I + N_M + 3); // neurotransmitter amount
  static constexpr unsigned int GRN_OUTPUT_NT_T = (N_G_I + N_M + 4); // neurotransmitter amount thresh
  static constexpr unsigned int GRN_OUTPUTS = (N_M + 5); // number of axon grn outputs
  static constexpr unsigned int GRN_REGULS = 3; // number of starting regulatory proteins

  // Eval problem parameters
  static constexpr unsigned int HOTDOGS = 40;
  static constexpr double ROBOT_SIZE = 1.0;
  static constexpr double ROBOT_SIGHT = 7.5;
  static constexpr double ROBOT_LIFE = 100.0;
  static constexpr unsigned int IMAGE_PROCESS_STEP = 20;

  // GA parameters
  static constexpr unsigned int GA_COLLABORATORS = 5;
  static constexpr unsigned int AXON_GA_POP = 20;
  static constexpr unsigned int SOMA_GA_POP = 20;
  static constexpr unsigned int NUM_POP = 200;
  static constexpr unsigned int GENERATIONS = 400;
  static constexpr unsigned int TOURNAMENT_SIZE = 5;
  static constexpr double MUTATION_RATE = 0.75;
  static constexpr double MOD_MUTATION_RATE = 0.8;
  static constexpr double ADD_MUTATION_RATE = 0.1;
  static constexpr double DEL_MUTATION_RATE = 0.1;
  static constexpr double CROSSOVER_RATE = 0.40;
  static constexpr unsigned int NUMBER_ELITES = 4;
  static constexpr unsigned int PROBLEM_EVAL_STEPS = 200;
  static constexpr double NOVELTY_MIN = 5.0;
};
#endif
