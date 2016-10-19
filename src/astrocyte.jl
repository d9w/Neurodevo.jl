# controller functions adding static astrocyte-like glial cells (ctype 2) to a soma (ctype 1) only network

function astro_morphogen_diff(morphogens::Vector{Float64}, dist::Float64, cell::CellInputs)
  diffs = -0.1*morphogens
  if cell.ctype == 1
    diffs[1] += cell.ntconc
  elseif cell.ctype == 2
    diffs[2] += 0.4*(morphogens[1] - morphogens[3])
    diffs[3] += 0.1*(morphogens[3] - morphogens[1])
  end
  diffs
end

function astro_synapse_formation(dist::Float64, acell::CellInputs, bcell::CellInputs)
  (acell.ctype == 1) && (bcell.ctype == 1) && (acell.cparams[1] == bcell.cparams[1]+1)
end

function astro_synapse_weight(reward::Float64, a_morphs::Vector{Float64}, b_morphs::Vector{Float64},
                              acell::CellInputs, bcell::CellInputs)
  weight = 0.0
  if acell.ctype == 1 && bcell.ctype == 1
    if a_morphs[1] > b_morphs[1]
      weight += 0.1
    elseif a_morphs[1] < b_morphs[1]
      weight -= 0.1
    end
    weight = weight + 0.25*a_morphs[2] + 0.25*b_morphs[2] - 0.25*a_morphs[3] - 0.25*a_morphs[3]
  end
  weight
end

function astrocyte_controller()
  rdivision = (morphogens::Vector{Float64}, cell::CellInputs)->false
  rchild_type = (morphogens::Vector{Float64}, cell::CellInputs)->rand(1:N_CTYPES)
  rchild_params = (morphogens::Vector{Float64}, ccell_type::Int64, pcell::CellInputs)->rand(1:N_MORPHS, N_PARAMS)
  rchild_position = (morphogens::Vector{Float64}, cell::CellInputs)->randn(N_D)
  rapoptosis = (morphogens::Vector{Float64}, cell::CellInputs)->false

  rcell_movement = (morphogens::Vector{Float64}, gradients::Array{Float64}, cell::CellInputs)->zeros(N_D)

  rsynapse_weight = (reward::Float64, a_morphs::Vector{Float64}, b_morphs::Vector{Float64}, acell::CellInputs,
                 bcell::CellInputs)->randn()

  rsynapse_output = (synapse_input::Float64, synapse_weight::Float64, cell::CellInputs)->rand()
  rnt_update = (synapse_input::Float64, synapse_output::Float64, cell::CellInputs)->randn()

  Controller(rdivision, rchild_type, rchild_params, rchild_position, rapoptosis, astro_morphogen_diff, rcell_movement,
             rsynapse_formation, rsynapse_weight, rsynapse_output, rnt_update)
end
