export Controller, Eval, Individual, Environment

# type Cell
#     ID::Int64
#     params::Vector{Float64}
#     state::Vector{Float64}
# end

# type Synpase
#     params::Vector{Float64}
#     state::Vector{Float64}
#     pre::Int64
#     post::Int64
# end

type Controller
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

type Individual
    cell_params::Array{Float64}
    cell_states::Array{Float64}
    synapse_params::Array{Float64}
    synapse_state::Array{Float64}
    synapse_index::Array{Float64}
    fitness_history::Array{Int64}
    fitness_values::Vector{Float64}
    controller::Controller
end

type Environment
    individuals::Array{Individual}
    next_task::Function
    mate::Function
    compete::Function
    dominates::Function
    new_state::Array{Function}
    reward::Array{Function}
    mutation::Float64
    crossover::Float64
    param_length::Int64
    state_length::Int64
end
