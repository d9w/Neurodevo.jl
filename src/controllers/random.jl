export rand_controller

function rand_controller(cfg::Dict)
    Random.seed!(cfg["seed"])
    nouts = output_lengths(cfg)
    cell_division(x::Array{Float64}) = rand(Bool)
    new_cell_params(x::Array{Float64}) = rand(nouts[2])
    cell_death(x::Array{Float64}) = rand(Bool)
    cell_state_update(x::Array{Float64}) = rand(nouts[4])
    cell_param_update(x::Array{Float64}) = rand(nouts[5])
    connect(x::Array{Float64}) = rand(Bool)
    new_conn_params(x::Array{Float64}) = rand(nouts[7])
    disconnect(x::Array{Float64}) = rand(Bool)
    conn_state_update(x::Array{Float64}) = rand(nouts[9])
    conn_param_update(x::Array{Float64}) = rand(nouts[10])

    Controller(cell_division, new_cell_params, cell_death,
               cell_state_update, cell_param_update,
               connect, new_conn_params, disconnect,
               conn_state_update, conn_param_update)
end
