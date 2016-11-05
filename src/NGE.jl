module NGE

include("types.jl")
include("constants.jl")
include("model.jl")
include("inits.jl")
# include("biocont.jl")
include("astrocyte.jl")
Model() = astrocyte_model()

export Cell, Synapse, Model, Controller, CellInputs, MORPHOGEN_REMAINDER,
  STEP_INIT, STEP_TOTAL, CHANGE_STOP, AXON_MAX, MIN_DIST, N_MORPHS, N_PARAMS,
  N_CTYPES, N_NSC, DIMS, N_D, NT_SPIKE, V_T, add_cell!, add_synapse!,
  apoptosis!, update_morphogens!, cell_division!, synapse_update!,
  cell_movement!, fire!, step!

module graphing
# include("graph.jl")
# include("random.jl")
end

module SNN
include("SNN.jl")
end

end # module NGE
