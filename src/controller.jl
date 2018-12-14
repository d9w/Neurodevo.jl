export Controller

struct Controller
    cell_division::Function # cell params -> bool
    new_cell_params::Function # cell params -> child params
    cell_death::Function # cell params -> bool
    cell_state_update::Function # cell params, state, input -> state, output
    cell_param_update::Function # cell params, state -> params

    connect::Function # cell 1 params, cc 2 params -> bool
    new_conn_params::Function # cell 1 params, cc 2 params -> conn params
    disconnect::Function # c1 params, c2 params, conn params, state -> bool
    conn_state_update::Function # conn params, state, cell input, conn input -> state, output
    conn_param_update::Function # c1 params, c2 params, conn params, state -> params
end

function Controller(cfg::Dict)
    nouts = output_lengths(cfg)
    cell_division(x::Array{Float64}) = false
    new_cell_params(x::Array{Float64}) = zeros(nouts[2])
    cell_death(x::Array{Float64}) = false
    cell_state_update(x::Array{Float64}) = zeros(nouts[4])
    cell_param_update(x::Array{Float64}) = zeros(nouts[5])
    connect(x::Array{Float64}) = false
    new_conn_params(x::Array{Float64}) = zeros(nouts[7])
    disconnect(x::Array{Float64}) = false
    conn_state_update(x::Array{Float64}) = zeros(nouts[9])
    conn_param_update(x::Array{Float64}) = zeros(nouts[10])

    Controller(cell_division, new_cell_params, cell_death,
               cell_state_update, cell_param_update,
               connect, new_conn_params, disconnect,
               conn_state_update, conn_param_update)
end
