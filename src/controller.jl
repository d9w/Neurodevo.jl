# default (random) controller and docstrings

include("types.jl")

"""
determines if a cell should divide
applied at each time step for each cell
return 1 bool
"""
function division(morphogens::Vector{Float64}, cell::CellInputs)
  rand(Bool)
end

"""
child cell type, int in 1:4
applied upon positive division decision
return 1 int in 1:N_CTYPES
"""
function child_type(morphogens::Vector{Float64}, cell::CellInputs)
  rand(1:N_CTYPES)
end

"""
child cell params
Applied upon positive division decision
return N_PARAMS int in 1:N_MORPHS
"""
function child_params(morphogens::Vector{Float64}, ccell_type::Int64, pcell::CellInputs)
  rand(1:N_MORPHS, N_PARAMS)
end

"""
child cell position, applied as a diff of the parent cell position
Applied upon positive division decision
return N_D floats
"""
function child_position(morphogens::Vector{Float64}, cell::CellInputs)
  rand(N_D)
end

"""
apoptosis, programmed cell death
applied for each cell every each time step
return 1 bool
"""
function apoptosis(morphogens::Vector{Float64}, cell::CellInputs)
  rand(Bool)
end

"""
morphogen diffusion
applied for each cell at each grid point every time step
return N_M floats
"""
function morphogen_diff(morphogens::Vector{Float64}, dist::Vector{Float64}, cell::CellInputs)
  randn(N_MORPHS)
end

"""
the diff in position of a cell
applied at each time step for each cell
return N_D floats
"""
function cell_movement(morphogens::Vector{Float64}, gradients::Array{Float64}, cell::CellInputs)
  randn(N_D)
end

"""
determines if a synapse is formed, true if formed
applied to each axon, soma pair not in a synapse at each timestep
return 1 bool
"""
function synapse_formation(dist::Float64, acell::CellInputs, bcell::CellInputs)
  rand(Bool)
end

"""
update to synaptic weight
applied to each possible synapse (soma,axon pair) at each timestep
return 1 float
"""
function synapse_weight(a_morphs::Vector{Float64}, b_morphs::Vector{Float64},
                        acell::CellInputs, bcell::CellInputs)
  randn()
end

"""
determines synaptic output based on inputs
applied to each receptor cell of a synapse at each timestep
return 1 float
"""
function nt_output(synapse_input::Float64, cell::CellInputs)
  randn()
end

"""
added to the cell's neurotransmitter value
applied to each cell at each timestep
return 1 float
"""
function nt_update(synapse_input::Float64, synapse_output::Float64, cell::CellInputs)
  randn()
end

function Controller()
  Controller(division, child_type, child_params, child_position, apoptosis, morphogen_diff, cell_movement,
             synapse_formation, synapse_weight, nt_output, nt_update)
end
