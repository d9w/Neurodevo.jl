using Distances

"""
determines if a cell should divide, true for division
applied at each time step for each cell
"""
function division(morphogens::Vector{Float64}, cell_type::Int64, cell_params::Vector{Int64}, bias::Float64)
  support_m = morphogens[cell_params[1]]
  thresh = [5.0, 1.0, 3.0, 0.5]
  support_m > thresh[cell_type]
end

"""
division branch direction, true for left
applied upon positive division decision
"""
function child_branch(morphogens::Vector{Float64}, cell_type::Int64, cell_params::Vector{Int64}, bias::Float64)
  probs = bias .< [0.25, 0.4, 0.9, 0.7]
  probs[cell_type]
end

"""
child cell type, int in 1:4
applied upon positive division decision
"""
function child_type(morphogens::Vector{Float64}, cell_type::Int64, cell_params::Vector{Int64}, branch_left::Bool)
  right = [2, 2, 0, 0]
  left = [3, 0, 4, 4]
  branch_left ? left[cell_type] : right[cell_type]
end

"""
child cell params
Applied upon positive division decision
"""
function child_params(morphogens::Vector{Float64}, pcell_type::Int64, pcell_params::Vector{Int64}, ccell_type::Int64,
                      bias::Float64)
  params = convert(Vector{Int64},ceil(((1:4)*16*bias)%4))
  if pcell_type == 1
    if ccell_type == 2
      params[1] = pcell_params[1]
      params[3] = pcell_params[2]
    elseif ccell_type == 3
      params[1] = pcell_params[2]
      params[2] = pcell_params[1]
    end
  elseif pcell_type == 3
    if ccell_type == 4
      params[4] = pcell_params[3]
    end
  end
  params
end

"""
morphogen diffusion
applied for each cell at each grid point
"""
function morphogen_diff(n_m::Int64, morphogen::Int64, cell_type::Int64, cell_params::Vector{Int64}, dist::Float64)
  diff = 0.0
  factor = 0.0
  if cell_type == 1 || cell_type == 3
    if cell_params[3] == morphogen
      factor = cell_params[4] / n_m
    end
  elseif cell_type == 2
    if cell_params[3] == morphogen
      factor = cell_params[4] / n_m
    elseif cell_params[1] == morphogen
      factor = -cell_params[2] / n_m
    end
  end
  if factor != 0.0
    diff = factor * exp(-10.0*(dist^2))
  end
  diff
end

"""
child cell position, applied as a diff of the parent cell position
Applied upon positive division decision
"""
function child_position(morphogens::Vector{Float64}, parent_cell_type::Int64, parent_cell_params::Vector{Int64},
                        world_dims::Array{Float64}, bias::Float64)
  support_m = morphogens[parent_cell_params[1]]
  dir = convert(Int64,ceil(length(world_dims) * bias))
  pos = zeros(world_dims)
  pos[dir] = 0.001*mean(world_dims)/max(support_m,1.0)
  pos
end

"""
the diff in position of a cell
applied at each time step for each cell
"""
function cell_movement(morphogens::Vector{Float64}, gradients::Array{Float64}, cell_type::Int64,
                       cell_params::Vector{Int64}, world_dims::Array{Float64})
  velocity = [0.05, 0.0, 0.0, 0.1]
  v = velocity[cell_type] * mean(world_dims)
  new_pos = zeros(world_dims)
  if v > 0.0
    # calculate new position
    follow = [1, 0, 0, 3]
    repel = [2, 0, 0, 4]
    follow_morph = cell_params[follow[cell_type]]
    repel_morph = cell_params[repel[cell_type]]
    new_pos += v .* vec(gradients[follow_morph,:] - gradients[repel_morph,:])
  end
  new_pos
end

"""
determines if a synapse is formed, true if formed
applied to each axon, soma pair not in a synapse at each timestep
"""
function synapse_formation(dist::Float64)
  dist < 0.01
end

"""
update to synaptic weight
applied to each possible synapse (soma,axon pair) at each timestep
"""
function synapse_weight(soma_morphs::Vector{Float64}, axon_morphs::Vector{Float64},
                        soma_params::Vector{Int64}, axon_params::Vector{Int64})
  soma_morphs[soma_params[1]]+axon_morphs[axon_params[1]]-soma_morphs[soma_params[2]]-axon_morphs[axon_params[2]]
end

"""
synaptic survival, true if the synapse survives
Applied to each synapse at each timestep
"""
function synapse_survival(synapse_weight::Float64)
  synapse.weight >= 0
end

"""
synaptic weight modification based on reward
applied to each synapse upon reward signal, requires supervisor
"""
function reward(soma_morphs::Vector{Float64}, axon_morphs::Vector{Float64},
                soma_params::Vector{Int64}, axon_params::Vector{Int64})
  # use morphogens to detect if synapse is being inhibited or enhanced
  # reinforce that behavior
  2*synapse_weight(soma_morphs, axon_morphs, soma_params, axon_params)
end
