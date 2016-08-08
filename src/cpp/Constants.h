// GRN parameters
#define ID_SIZE 32  // Protein ID Max value
#define N_M 0 // number of morphogens
#define N_D 3 // number of dimensions TODO: very hardcoded currently

#define GRN_EVO_STEPS 1
#define SOMA_GRN_EVO_STEPS 1
#define AXON_GRN_EVO_STEPS 1
#define AXON_DIVISION_REDUCTION 0.2
#define AXON_MAX_NUMBER 5 // maybe get rid of this, but for testing it is necessary
#define SOMA_AXON_INPUT_THRESH 1.0
#define T_ACTION_MIN 5
#define T_ACTION_MAX 20

// GRN INPUTS
#define GRN_INPUT_X N_M // position x
#define GRN_INPUT_Y (N_M + 1) // position y
#define GRN_INPUT_Z (N_M + 2) // position z
#define GRN_INPUT_NEUROTRANSMITTER (N_M + 3) // neurotransmitter concentration (post-synaptic in axon)
#define GRN_INPUT_COMM (N_M + 4) // communication input
#define GRN_INPUT_DIVISION (N_M + 5) // division concentration, 1 in soma
#define GRN_INPUT_REWARD (N_M + 6) // division concentration, 1 in soma
#define GRN_INPUTS (N_M + 7) // number of axon grn inputs
#define N_G_I GRN_INPUTS
// grn outpus
#define GRN_OUTPUT_COMM (N_G_I + N_M) // commuincation output
#define GRN_OUTPUT_F_NULL (N_G_I + N_M + 1) // firing threshold/null
#define GRN_OUTPUT_FT_NULL (N_G_I + N_M + 2) // firing threshold threshold/null
#define GRN_OUTPUT_NT (N_G_I + N_M + 3) // neurotransmitter amount
#define GRN_OUTPUT_NT_T (N_G_I + N_M + 4) // neurotransmitter amount threshold
#define GRN_OUTPUTS (N_M + 5) // number of axon grn outputs
#define GRN_REGULS 3 // number of starting regulatory proteins

/*
// AXON GRN INPUTS, after all morphogens
#define axon_grn_input_x n_m // position x
#define axon_grn_input_y (n_m + 1) // position y
#define axon_grn_input_z (n_m + 2) // position z
#define axon_grn_input_neurotransmitter (n_m + 3) // post-synaptic neurotransmitter concentration
#define axon_grn_input_soma (n_m + 4) // soma connection input
#define axon_grn_input_division (n_m + 5) // division concentration
#define axon_grn_inputs (n_m + 6) // number of axon grn inputs
#define a_g_i axon_grn_inputs
// axon grn outpus, after all morphogens
#define axon_grn_output_branch (a_g_i + n_m) // branch action
#define axon_grn_output_apoptosis (a_g_i + n_m + 1) // apoptosis
#define axon_grn_output_nothing (a_g_i + n_m + 2) // do nothing action
#define axon_grn_output_neurotransmitter (a_g_i + n_m + 3) // neurotransmitter amount
#define axon_grn_output_neurotransmitter_thresh (a_g_i + n_m + 4) // neurotransmitter amount threshold
#define axon_grn_output_soma (a_g_i + n_m + 5) // soma connection output
#define axon_grn_outputs (n_m + 6) // number of axon grn outputs
#define axon_grn_reguls 3 // number of starting regulatory proteins

// SOMA GRN INPUTS, after all morphogens
#define SOMA_GRN_INPUT_X N_M // position x
#define SOMA_GRN_INPUT_Y (N_M + 1) // position y
#define SOMA_GRN_INPUT_Z (N_M + 2) // position z
#define SOMA_GRN_INPUT_NEUROTRANSMITTER (N_M + 3) // neurotransmitter concentration
#define SOMA_GRN_INPUT_AXON (N_M + 4) // axon connection input (summed and cut)
#define SOMA_GRN_INPUTS (N_M + 5)
#define S_G_I SOMA_GRN_INPUTS
// SOMA GRN OUTPUTS, after all morphogens except the central morphogen
#define SOMA_GRN_OUTPUT_MORPHOGEN_THRESH (S_G_I + N_M) // morphogen output threshold
#define SOMA_GRN_OUTPUT_AXON (S_G_I + N_M + 1) // axon connection output
#define SOMA_GRN_OUTPUT_FIRING (S_G_I + N_M + 2) // firing threshold
#define SOMA_GRN_OUTPUT_FIRING_THRESH (S_G_I + N_M + 3) // firing threshold
#define SOMA_GRN_OUTPUTS (N_M + 4)
#define SOMA_GRN_REGULS 3
*/

// Environment parameters
#define X_SIZE 8
#define Y_SIZE 8
#define Z_SIZE 8

//GA parameters
#define GA_COLLABORATORS 5 // the number of collaborators chosen
#define AXON_GA_POP 20
#define SOMA_GA_POP 20
#define NUM_POP 20
#define MAX_ITER 400
#define TOURNEY_SIZE 5
#define MUTATION_RATE 0.10
#define ADD_MUTATION_RATE 0.10
#define DELETE_MUTATION_RATE 0.05
#define CROSSOVER_RATE 0.40
#define GA_EVAL_STEPS 2500
#define IMAGE_PROCESS_STEP 20

#define HOTDOGS 40
