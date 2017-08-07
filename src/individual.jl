mutable struct Controller
    cell_division::Function
    child_parameters::Function
    child_state::Function
    synpase_formation::Function
    synapse_parameters::Function
    synapse_state::Function
    cell_update::Function
    cell_death::Function
    synapse_update::Function
    synapse_death::Function
    input::Function
    output::Function
    reward::Function
end

mutable struct Individual
    cell_params::Array{Float64}
    cell_states::Array{Float64}
    synapse_params::Array{Float64}
    synapse_state::Array{Float64}
    synapse_index::Array{Float64}
    fitness_history::Array{Int64}
    fitness_values::Vector{Float64}
    controller::Controller
end

