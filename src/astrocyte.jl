# controller functions adding static astrocyte-like glial cells (ctype 2) to a soma (ctype 1) only network

function astro_morphogen_diff(morphogens::Vector{Float64}, dist::Vector{Float64}, cell::CellInputs)
  diffs = -0.1*morphogens
  euc = sqrt(mapreduce(x->x^2, +, dist))
  if cell.ctype == 1 && euc < 0.05
    # activity information for stdp-like modulation
    diffs[1] += 0.1*cell.ntconc
  elseif cell.ctype == 2
    # inactivity information for stdp-like modulation
    if euc < 0.05
      diffs[2] += 0.1*(1.0-cell.ntconc)
    end
    if (dist[3] < 0.1 && dist[3] > -0.1) && euc > 0.01
      # homeostasis within layers
      diffs[3] += 0.1*cell.ntconc
    end
    if dist[3] < -0.1
      # backprop communication about error and local firing
      diffs[4] += - 0.1 * (morphogens[4] + 0.05*(morphogens[1]-morphogens[2])) / (1-exp(euc))
    end
  end
  diffs
end

function astro_synapse_formation(dist::Float64, acell::CellInputs, bcell::CellInputs)
  (acell.ctype == 1) && ((bcell.ctype == 1) && (acell.cparams[1] == bcell.cparams[1]+1)
                         || (bcell.ctype == 2) && dist < 0.02)
end

function astro_synapse_weight(a_morphs::Vector{Float64}, b_morphs::Vector{Float64},
                              acell::CellInputs, bcell::CellInputs)
  weight = 0.0
  if (acell.ctype == 1) && (bcell.ctype == 1)
    if (a_morphs[1] > b_morphs[2]) && (b_morphs[1] > a_morphs[2])
      weight += 0.1 + 0.1*b_morphs[4]
    elseif (a_morphs[2] > b_morphs[1]) && (b_morphs[2] > a_morphs[1])
      weight -= 0.2
    end
  end
  weight
end

function astro_nt_output(synapse_input::Float64, cell::CellInputs)
  # ReLu
  out = 0.0
  if synapse_input > 0.0
    out = synapse_input
  end
  out
end

function astro_nt_update(synapse_input::Float64, synapse_output::Float64, cell::CellInputs)
  synapse_input - cell.ntconc
end

function astrocyte_controller()
  rdivision = (morphogens::Vector{Float64}, cell::CellInputs)->false
  rchild_type = (morphogens::Vector{Float64}, cell::CellInputs)->1
  rchild_params = (morphogens::Vector{Float64}, ccell_type::Int64, pcell::CellInputs)->ones(N_PARAMS)
  rchild_position = (morphogens::Vector{Float64}, cell::CellInputs)->zeros(N_D)
  rapoptosis = (morphogens::Vector{Float64}, cell::CellInputs)->false
  # astro_morphogen_diff
  rcell_movement = (morphogens::Vector{Float64}, gradients::Array{Float64}, cell::CellInputs)->zeros(N_D)
  # astro_synapse_formation
  # astro_synapse_weight

  Controller(rdivision, rchild_type, rchild_params, rchild_position, rapoptosis, astro_morphogen_diff, rcell_movement,
             astro_synapse_formation, astro_synapse_weight, astro_nt_output, astro_nt_update)
end
