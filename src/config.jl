export Config

struct Config
    seed::Int64
    n_cell_state::Int64
    n_cell_params::Int64
    n_branch_state::Int64
    n_branch_params::Int64
    memory_max::Int64 # 1024^3 (1 GiB), in bytes
    time_max::Float64 # 10.0 s, in seconds
end

function Config()
    seed = 0
    n_cell_state = 5
    n_cell_params = 5
    n_branch_state = 5
    n_branch_params = 5
    memory_max = 1024^3
    time_max = 10.0

    Config(seed, n_cell_state, n_cell_params, n_branch_state, n_branch_params,
           memory_max, time_max)
end

# function Config(genes::Array{Float64})
# end
