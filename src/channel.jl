export Branch, Channel

struct Branch
    state::Array{Float64}
    params::Array{Float64}
    connections::Array{Int64}
end

mutable struct Channel
    input::Float64
    branches::Array{Branch}
    neuron_id::Int64
    
end

function Branch(cfg::Dict, params::Array{Float64})
    Branch(zeros(cfg["n_branch_state"]),
           params,
           zeros(Int64, cfg["n_connections"]))
end

