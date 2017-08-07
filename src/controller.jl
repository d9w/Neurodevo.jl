export Controller

mutable struct Controller
    cell_division::Function
    child_parameters::Function
    child_state::Function
    cell_update::Function
    cell_death::Function
    synapse_formation::Function
    synapse_parameters::Function
    synapse_state::Function
    synapse_update::Function
    synapse_death::Function
    input::Function
    output::Function
    reward::Function
end

function Controller(cell_state::Int64, cell_params::Int64, syn_state::Int64,
                    syn_params::Int64)
    cell_division = (x...)->false
    child_parameters = (x...)->zeros(cell_params)
    child_state = (x...)->zeros(cell_state)
    cell_update = (x...)->zeros(cell_state)
    cell_death = (x...)->false
    synapse_formation = (x...)->false
    synapse_parameters = (x...)->zeros(syn_params)
    synapse_state = (x...)->zeros(syn_state)
    synapse_update = (x...)->zeros(syn_state)
    synapse_death = (x...)->false
    input = (x...)->zeros(3)
    output = (x...)->zeros(3)
    reward = (x...)->zeros(3)
    Controller(cell_division, child_parameters, child_state, cell_update,
               cell_death, synapse_formation, synapse_parameters, synapse_state,
               synapse_update, synapse_death, input, output, reward)
end
