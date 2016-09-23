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
  spike::Function
end

type CellInputs
  id::Float64
  p_id::Float64
  ctype::Int64
  params::Vector{Int64}
  ntconc::Float64
end

"""
determines if a cell should divide
applied at each time step for each cell
return 1 bool
"""
function division(morphogens::Vector{Float64}, cell::CellInputs)
  if cell.id > 0.5
    otherm = [4,2,3,3]
    lower_scale = [1.6,0.8,0.8,0.8]
    upper_scale = [1.8,1.0,1.1,1.1]
    morph = morphogens[cell.params[otherm[cell.ctype]]]
    mm = mean(morphogens)
    return (morph > lower_scale[cell.ctype] * mm) && (morph < upper_scale[cell.ctype] * mm)
  else
    return false
  end
end

"""
child cell type, int in 1:4
applied upon positive division decision
return 1 int in 1:N_CTYPES
"""
function child_type(morphogens::Vector{Float64}, cell::CellInputs)
  ctypes = [3,2,4,4]
  ctype = ctypes[cell.ctype]
  if cell.ctype == 1
    if (morphogens[cell.params[1]] + morphogens[cell.params[2]]) < 1.5*morphogens[cell.params[3]]
      ctype = 2
    end
  end
  ctype
end

"""
child cell params
Applied upon positive division decision
return N_PARAMS int in 1:N_MORPHS
"""
function child_params(morphogens::Vector{Float64}, cell::CellInputs, ccell.ctype::Int64)
  params = Array{Int64}(4)
  smorph = +(morphogens...)
  n_m = length(morphogens)
  maxv, maxm = findmax(morphogens)
  minv, minm = findmin(morphogens)
  if ccell.ctype == 1
    params[1] = maxm + minm
    params[2] = maxm
    params[3] = minm
    params[4] = 1+convert(Int64,floor((n_m-1.0)*morphogens[params[3]]/smorph))
  elseif ccell.ctype == 2
    if pcell.ctype == 1
      params[1] = pcell.params[morphogens[pcell.params[1]] > morphogens[pcell.params[2]] ? 1 : 2]
      params[2] = 1+convert(Int64,floor((n_m-1.0)*(smorph-morphogens[params[1]])/smorph))
      params[3] = pcell.params[3]
      params[4] = 1+convert(Int64,floor((n_m-1.0)*morphogens[params[3]]/smorph))
    else
      params[1] = maxm
      params[2] = 1+convert(Int64,floor((n_m-1.0)*(smorph-morphogens[params[1]])/smorph))
      params[3] = minm
      params[4] = 1+convert(Int64,floor((n_m-1.0)*morphogens[params[3]]/smorph))
    end
  elseif ccell.ctype == 3
    params[1] = maxm
    params[2] = minm
    params[3] = mod(pcell.params[4] + minm, n_m) + 1
    params[4] = mod(pcell.params[3] + maxm, n_m) + 1
  elseif ccell.ctype == 4
    params[1] = minm
    params[2] = maxm
    params[3] = mod(pcell.params[4] + minm, n_m) + 1
    params[4] = mod(pcell.params[3] + maxm, n_m) + 1
  end
  @assert all(x->(x>=1)&&(x<=4), params)
  #TODO: mod all by N_M, in model
  params
end

"""
child cell position, applied as a diff of the parent cell position
Applied upon positive division decision
return N_D floats
"""
function child_position(morphogens::Vector{Float64}, pcell::CellInputs)
  support = [1,3,1,1]
  ant = [2,1,2,2]
  maxm = maximum(morphogens)
  supportm = morphogens[parent_cell.params[support[parent_cell.ctype]]]
  antm = morphogens[parent_cell.params[ant[parent_cell.ctype]]]
  0.01*(sin(supportm/maxm*(1:N_D)) - cos(antm/maxm*(1:N_D)))
end

"""
apoptosis, programmed cell death
applied for each cell every each time step
return 1 bool
"""
function apoptosis(morphogens::Vector{Float64}, cell::CellInputs)
  if cell.id > 0.5
    otherm = [1,3,1,1]
    lower_scale = [0.05,1.4,0.1,0.5]
    upper_scale = [0.06,5.0,0.3,0.7]
    morph = morphogens[cell.params[otherm[cell.ctype]]]
    mm = mean(morphogens)
    return (morph > lower_scale[cell.ctype] * mm) && (morph < upper_scale[cell.ctype] * mm)
  else
    return false
  end
end

"""
morphogen diffusion
applied for each cell at each grid point every time step
return N_M floats
"""
function morphogen_diff(n_m::Int64, morphogen::Int64, dist::Float64, cell::CellInputs)
  diff = 0.0
  factor = 0.0
  if cell.ctype == 1 || cell.ctype == 3
    if cell.params[3] == morphogen
      factor = cell.params[4] / n_m
    end
  elseif cell.ctype == 2
    if cell.params[3] == morphogen
      factor = cell.params[4] / n_m
    elseif cell.params[1] == morphogen
      factor = -cell.params[2] / n_m
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
return N_D floats
"""
function cell_movement(morphogens::Vector{Float64}, gradients::Array{Float64}, cell::CellInputs)
  velocity = [0.05, 0.005, 0.01, 0.1]
  v = velocity[cell.ctype]
  new_pos = zeros(N_D)
  if v > 0.0
    mg = maximum(gradients)
    follow = [1, 1, 1, 3]
    repel = [2, 3, 2, 4]
    follow_morph = cell.params[follow[cell.ctype]]
    repel_morph = cell.params[repel[cell.ctype]]
    new_pos += v/mg .* vec(gradients[follow_morph,:] - gradients[repel_morph,:])
  end
  new_pos
end

"""
determines if a synapse is formed, true if formed
applied to each axon, soma pair not in a synapse at each timestep
return 1 bool
"""
function synapse_formation(dist::Float64, acell::CellInput, bcell::CellInputs)
  form = false
  if acell.ctype == 4
    if acell.p_id == bcell.id
      form = true
    else
      if dist < 0.1 && bcell.ctype == 3
        form = true
      end
    end
  end
  form
end

"""
update to synaptic weight
applied to each possible synapse (soma,axon pair) at each timestep
return 1 float
"""
function synapse_weight(a_morphs::Vector{Float64}, b_morphs::Vector{Float64}, acell::CellInputs, bcell::CellInputs,
                        reward::Float64)
  weight = 0.0
  if acell.ctype == 3
    if bcell.ctype == 4
      weight = (1.0 + reward) * ((a_morphs[acell.params[1]]+b_morphs[bcell.params[1]])
                                 -(a_morphs[acell.params[2]]+b_morphs[bcell.params[2]]))
    end
  end
  weight
end

"""
determines if a cell spikes
applied to each receptor cell of a synapse at each timestep
return 1 bool
"""
function spike(cell::CellInputs)
  (cell.ctype == 3 && cell.ntconc > 1.0) ||  cell.ctype == 4 && cell.ntconc > 0.0)
end

function Controller()
  Controller(division, child_type, child_params, child_position, apoptosis, morphogen_diff, cell_movement,
             synapse_formation, synapse_weight, spike)
end
