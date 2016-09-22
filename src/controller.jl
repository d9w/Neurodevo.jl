using Distances

type Controller
  division::Function
  child_type::Function
  child_params::Function
  child_position::Function
  apoptosis::Function
  morphogen_diff::Function
  cell_movement::Function
  synapse_formation::Function
  synapse_weight::Function
  synapse_survival::Function
end

"""
determines if a cell should divide
applied at each time step for each cell
"""
function division(morphogens::Vector{Float64}, cell_type::Int64, cell_params::Vector{Int64}, cell_velocity::Float64)
  stationary = cell_velocity < 0.1
  otherm = [4,2,3,3]
  stationary && morphogens[cell_params[otherm[cell_type]]] > mean(morphogens)
end

"""
child cell type, int in 1:4
applied upon positive division decision
"""
function child_type(morphogens::Vector{Float64}, cell_type::Int64, cell_params::Vector{Int64})
  ctypes = [3,2,4,4]
  ctype = ctypes[cell_type]
  if cell_type == 1
    if morphogens[cell_params[3]] < 1/5*(morphogens[cell_params[1]] + morphogens[cell_params[2]])
      ctype = 2
    end
  end
  ctype
end

"""
child cell params
Applied upon positive division decision
"""
function child_params(morphogens::Vector{Float64}, pcell_type::Int64, pcell_params::Vector{Int64}, ccell_type::Int64)
  params = Array{Int64}(4)
  smorph = +(morphogens...)
  n_m = length(morphogens)
  maxv, maxm = findmax(morphogens)
  minv, minm = findmin(morphogens)
  if ccell_type == 1
    params[1] = ((maxm + minm) % n_m) + 1
    params[2] = maxm
    params[3] = minm
    params[4] = 1+convert(Int64,floor((n_m-1.0)*morphogens[params[3]]/smorph))
  elseif ccell_type == 2
    if pcell_type == 1
      params[1] = pcell_params[morphogens[pcell_params[1]] > morphogens[pcell_params[2]] ? 1 : 2]
      params[2] = 1+convert(Int64,floor((n_m-1.0)*(smorph-morphogens[params[1]])/smorph))
      params[3] = pcell_params[3]
      params[4] = 1+convert(Int64,floor((n_m-1.0)*morphogens[params[3]]/smorph))
    else
      params[1] = maxm
      params[2] = 1+convert(Int64,floor((n_m-1.0)*(smorph-morphogens[params[1]])/smorph))
      params[3] = minm
      params[4] = 1+convert(Int64,floor((n_m-1.0)*morphogens[params[3]]/smorph))
    end
  elseif ccell_type == 3
    params[1] = maxm
    params[2] = minm
    params[3] = pcell_params[3]
    params[4] = pcell_params[4]
  elseif ccell_type == 4
    params[1] = minm
    params[2] = maxm
    params[3] = ((pcell_params[3] + floor(Int64, n_m/2.0)) % n_m) + 1
    params[4] = pcell_params[3]
  end
  params
end

"""
child cell position, applied as a diff of the parent cell position
Applied upon positive division decision
"""
function child_position(morphogens::Vector{Float64}, parent_cell_type::Int64, parent_cell_params::Vector{Int64})
  support = [1,3,1,1]
  ant = [2,1,2,2]
  maxm = maximum(morphogens)
  supportm = morphogens[parent_cell_params[support[parent_cell_type]]]
  antm = morphogens[parent_cell_params[ant[parent_cell_type]]]
  0.1*(supportm*sin(1:N_D) - antm*cos(1:N_D))/maxm
end

"""
apoptosis, programmed cell death
applied for each cell every each time step
"""
function apoptosis(morphogens::Vector{Float64}, cell_type::Int64, cell_params::Vector{Int64}, cell_velocity::Float64)
  stationary = cell_velocity < 0.1
  otherm = [4,2,3,3]
  scale = [0.05,1.0,0.1,0.3]
  stationary && (morphogens[cell_params[otherm[cell_type]]] < scale[cell_type]*mean(morphogens))
end

"""
morphogen diffusion
applied for each cell at each grid point every time step
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
    diff = factor * exp(-(dist^2))
  end
  diff
end

"""
the diff in position of a cell
applied at each time step for each cell
"""
function cell_movement(morphogens::Vector{Float64}, gradients::Array{Float64}, cell_type::Int64,
                       cell_params::Vector{Int64})
  velocity = [0.05, 0.0, 0.0, 0.1]
  v = velocity[cell_type]
  new_pos = zeros(N_D)
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
                        soma_params::Vector{Int64}, axon_params::Vector{Int64}, reward::Float64)
  (1.0 + reward) * (soma_morphs[soma_params[1]]+axon_morphs[axon_params[1]]
                    -soma_morphs[soma_params[2]]-axon_morphs[axon_params[2]])
end

"""
synaptic survival, true if the synapse survives
Applied to each synapse at each timestep
"""
function synapse_survival(synapse_weight::Float64)
  synapse_weight >= 0
end

function Controller()
  Controller(division, child_type, child_params, child_position, apoptosis, morphogen_diff, cell_movement,
             synapse_formation, synapse_weight, synapse_survival)
end
