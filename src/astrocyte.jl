# controller functions adding static astrocyte-like glial cells (ctype 2) to a soma (ctype 1) only network

function astro_morphogen_diff(cthresh::Float64, beta::Float64, hstat::Float64, locd::Float64)
  return function morphogen_diff(morphogens::Vector{Float64}, dist::Vector{Float64}, cell::CellInputs)
    # diffs = -0.1*morphogens # TODO: absorption should be global
    diffs = zeros(morphogens)
    euc = sqrt(mapreduce(x->x^2, +, dist))
    if cell.ctype == 2
      # astro ntconc is corresponding activity over time
      diffs[1] += (cell.ntconc-cthresh)*exp(-beta*euc)
      diffs[2] -= (cell.ntconc-cthresh)*exp(-beta*euc)
      if (dist[3] == 0.0) && euc > locd
        # homeostasis within layers
        diffs[3] += hstat*(cell.ntconc-cthresh)
      end
      if dist[3] < 0.0
        # backprop communication about error and local firing
        diffs[4] += (morphogens[4] + (cell.ntconc-cthresh)) * exp(-beta*euc)
      end
    end
    diffs
  end
end

function astro_synapse_formation(locd)
  return function synapse_formation(dist::Float64, acell::CellInputs, bcell::CellInputs)
    (acell.ctype == 1) && ((bcell.ctype == 1) && (acell.cparams[1] == bcell.cparams[1]+1)
                          || (bcell.ctype == 2) && dist < locd)
    # TODO: also try the following:
    # (acell.ctype == 1) && (acell.cparams[1] == bcell.cparams[1]+1)
  end
end

function astro_synapse_weight(scalef)
  return function synapse_weight(a_morphs::Vector{Float64}, b_morphs::Vector{Float64},
                                acell::CellInputs, bcell::CellInputs)
    weight = 0.0
    if (acell.ctype == 1) && (bcell.ctype == 1)
      # correlation
      corr = (a_morphs[1] + a_morphs[2]) - (b_morphs[1] + b_morphs[2])
      # anti-correlation
      acorr = (b_morphs[1] + b_morphs[2]) - (a_morphs[1] + a_morphs[2])
      if acorr != 0.0
        weight = scalef * (corr / acorr) * ((1+b_morphs[4]) / (1+b_morphs[3]))
      end
    end
    weight
  end
end

function astro_nt_output(vt::Float64=2.0)
  return function nt_output(synapse_input::Float64, cell::CellInputs)
    out = 0.0
    if cell.ctype == 1 # soma
      # ReLu
      if synapse_input > 0.0
        out = synapse_input
      end
      # LIF
      # if cell.ntcont + synapse_input > vt
      #   out = vt
      # end
    end # astrocytes have no output synaptic connections, for now
    out
  end
end

function astro_nt_update(vt::Float64=2.0, leak::Float64=0.01, rate::Float64=0.3)
  return function nt_update(synapse_input::Float64, synapse_output::Float64, cell::CellInputs)
    update = 0.0
    if cell.ctype == 1 # soma
      # ReLu, reset to most recent output
      update = synapse_output - cell.ntconc
      # LIF, update on fire
      # if synapse_output == vt
      #   update = -cell.ntconc
      # else
      #   update = synapse_input - leak*cell.ntconc
      # end
    else # astrocyte
      update = rate * (synapse_input - cell.ntconc)
    end
    update
  end
end

function astrocyte_controller(cthresh::Float64=1.0, beta::Float64=5.0, hstat::Float64=0.5, locd::Float64=0.01,
                              scalef::Float64=0.1, vt::Float64=2.0, leak::Float64=0.01, rate::Float64=0.3)
  division = (morphogens::Vector{Float64}, cell::CellInputs)->false
  child_type = (morphogens::Vector{Float64}, cell::CellInputs)->1
  child_params = (morphogens::Vector{Float64}, ccell_type::Int64, pcell::CellInputs)->ones(N_PARAMS)
  child_position = (morphogens::Vector{Float64}, cell::CellInputs)->zeros(N_D)
  apoptosis = (morphogens::Vector{Float64}, cell::CellInputs)->false
  morphogen_diff = astro_morphogen_diff(cthresh, beta, hstat, locd)
  cell_movement = (morphogens::Vector{Float64}, gradients::Array{Float64}, cell::CellInputs)->zeros(N_D)
  synapse_formation = astro_synapse_formation(locd)
  astro_synapse_weight = astro_synapse_weight(scalef)
  nt_output = astro_nt_output(vt)
  nt_update = astro_nt_update(vt, leak, rate)

  Controller(division, child_type, child_params, child_position, apoptosis, morphogen_diff, cell_movement,
             synapse_formation, synapse_weight, nt_output, nt_update)
end
