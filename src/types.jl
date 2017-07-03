export Cell, Synapse, Controller, Eval, Individual, Environment

type Cell
    state::Vector{Float64}
    params::Vector{Float64}
    ID::Int64
end

type Synpase
    state::Vector{Float64}
    params::Vector{Float64}
    C1::Int64
    C2::Int64
end

type Controller
    cell_division::Function
    child_parameters::Function
    child_state::Function
    synpase_formation::Function
    synapse_parameters::Function
    synapse_state::Function
    self_update::Function
    cell_update::Function
    cell_death::Function
    synapse_update::Function
    synapse_death::Function
    global_input::Function
    global_output::Function
    global_reward::Function
end

type Eval
    age::Int64
    fitness::Float64
    task_id::Int64
end

type Individual
    cells::Array{Cell}
    synapses::Array{Synapse}
    controller::Controller
    fitness::Array{Eval}
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
end
