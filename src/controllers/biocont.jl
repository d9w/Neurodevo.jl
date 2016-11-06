using Distances
# Rules for a biology based controller

"""
"""
function division(morphogens::Vector{Float64}, cell::CellInputs)
  otherm = [4,2,3,3]
  lower_scale = [1.6,0.8,0.8,0.8]
  upper_scale = [1.8,1.0,1.1,1.1]
  morph = morphogens[cell.params[otherm[cell.ctype]]]
  mm = mean(morphogens)
  return (morph > lower_scale[cell.ctype] * mm) && (morph < upper_scale[cell.ctype] * mm)
end

"""
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
"""
function child_params(morphogens::Vector{Float64}, ccell_type::Int64, pcell::CellInputs)
  params = Array{Int64}(4)
  smorph = +(morphogens...)
  n_m = length(morphogens)
  maxv, maxm = findmax(morphogens)
  minv, minm = findmin(morphogens)
  if ccell_type == 1
    params[1] = maxm + minm
    params[2] = maxm
    params[3] = minm
    params[4] = 1+convert(Int64,floor((n_m-1.0)*morphogens[params[3]]/smorph))
  elseif ccell_type == 2
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
  elseif ccell_type == 3
    params[1] = maxm
    params[2] = minm
    params[3] = mod(pcell.params[4] + minm, n_m) + 1
    params[4] = mod(pcell.params[3] + maxm, n_m) + 1
  elseif ccell_type == 4
    params[1] = minm
    params[2] = maxm
    params[3] = mod(pcell.params[4] + minm, n_m) + 1
    params[4] = mod(pcell.params[3] + maxm, n_m) + 1
  end
  for p in eachindex(params)
    params[p] = mod(params[p], n_m)
  end
  params[params.==0]=n_m
  params
end

"""
"""
function child_position(morphogens::Vector{Float64}, cell::CellInputs)
  support = [1,3,1,1]
  ant = [2,1,2,2]
  maxm = maximum(morphogens)
  supportm = morphogens[cell.params[support[cell.ctype]]]
  antm = morphogens[cell.params[ant[cell.ctype]]]
  0.01*(sin(supportm/maxm*(1:N_D)) - cos(antm/maxm*(1:N_D)))
end

"""
"""
function apoptosis(morphogens::Vector{Float64}, cell::CellInputs)
  otherm = [1,3,1,1]
  lower_scale = [0.05,1.4,0.1,0.5]
  upper_scale = [0.06,5.0,0.3,0.7]
  morph = morphogens[cell.params[otherm[cell.ctype]]]
  mm = mean(morphogens)
  return (morph > lower_scale[cell.ctype] * mm) && (morph < upper_scale[cell.ctype] * mm)
end

"""
"""
function morphogen_diff(morphogens::Vector{Float64}, dist::Float64, cell::CellInputs)
  diffs = zeros(length(morphogens))
  for morphogen in eachindex(diffs)
    diff = 0.0
    factor = 0.0
    if cell.ctype == 1 || cell.ctype == 3
      if cell.params[3] == morphogen
        factor = cell.params[4] / N_MORPHS
      end
    elseif cell.ctype == 2
      if cell.params[3] == morphogen
        factor = cell.params[4] / N_MORPHS
      elseif cell.params[1] == morphogen
        factor = -cell.params[2] / N_MORPHS
      end
    end
    if factor != 0.0
      diff = factor * exp(-(dist^2))
    end
    diffs[morphogen] = diff
  end
  diffs
end

"""
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
"""
function synapse_formation(dist::Float64, acell::CellInputs, bcell::CellInputs)
  form = false
  if acell.ctype == 4
    if dist < 0.1 && bcell.ctype == 3
      form = true
    end
  end
  form
end

"""
"""
function synapse_weight(a_morphs::Vector{Float64}, b_morphs::Vector{Float64},
                        acell::CellInputs, bcell::CellInputs)
  weight = 0.0
  if acell.ctype == 3
    if bcell.ctype == 4
      weight = 1.0 + ((a_morphs[acell.params[1]]+b_morphs[bcell.params[1]])
                      -(a_morphs[acell.params[2]]+b_morphs[bcell.params[2]]))
    end
  end
  weight
end

"""
"""
function nt_output(synapse_input::Float64, cell::CellInputs)
  out = 0.0
  if (cell.ctype == 3 && cell.ntconc > 2.0) ||  cell.ctype == 4
    out = cell.nt_conc
  end
  out
end

"""
"""
function nt_update(synapse_input::Float64, synapse_output::Float64, cell::CellInputs)
  update = 0.0
  if cell.ctype == 3
    if cell.ntconc > 2.0
      update = -cell.ntconc
    else
      update = synapse_input
    end
  end
  update
end

function Controller()
  Controller(division, child_type, child_params, child_position, apoptosis, morphogen_diff, cell_movement,
             synapse_formation, synapse_weight, nt_output, nt_update)
end
